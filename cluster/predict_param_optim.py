# /usr/bin/env python3.10
# -*- coding: utf-8 -*-
# Problem definition and a pipeline for sensitivity analysis
import argparse
import itertools
import logging
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pygmo as pg
import scipy.linalg
import statsmodels.api as sm
from pykmd.kmd import StreamingKmdModel
from pykmd.libs.datclass import Procedure, RankMode, RefineRitz, Scaling, WeightsScheme
from pykmd.libs.lift import blocked_hankel_lift
from sklearn.preprocessing import MinMaxScaler

# IMPORTANT: add grip-emg-optimize to PYTHONPATH
from optim.problem import SignalProcessingParams


def predict_grip_pykmd(
    grip_approx: np.ndarray,
    num_delays: int,
    predict_horizon: int,
    thinning_step: int = 1,
    threshold: list = [-10, 600],
    force_rank: int = 2,
) -> np.ndarray:
    """
    Predict grip force using PyKMD library
    """
    # Add transformations
    grip_approx_obs = np.vstack((grip_approx,))
    # Number of observables
    obs_num = grip_approx_obs.shape[0]
    # Add time delays on observables
    grip_approx_lift = block_hankel_lift(
        grip_approx_obs, delay=num_delays, thin_hankel=thinning_step
    )
    # Add interactions between all time delays
    # Preallocate grip_approx_lift_td_interact for storing interactions between time delays
    grip_approx_lift_td_interact = np.empty(
        (
            grip_approx_lift.shape[0] + math.comb(grip_approx_lift.shape[0], 2),
            grip_approx_lift.shape[1],
        )
    )
    # Store grip_approx_lift at the end of grip_approx_lift_td_interact
    grip_approx_lift_td_interact[-grip_approx_lift.shape[0] :, :] = grip_approx_lift
    for i, (mul1, mul2) in enumerate(itertools.combinations(grip_approx_lift, 2)):
        # Populate rows in b with pairwise multiplication
        grip_approx_lift_td_interact[i, :] = np.log(mul1 + 10) * np.log(mul2 + 10)

    model = StreamingKmdModel(observables=grip_approx_lift_td_interact)
    model.initial_hankel_lifting(hankel_delays=0)
    model.initial_qr_compress()
    # print(f'Lifted observables shape: {model.lifted_observables.shape}')

    # Perform QR-DMD
    model.compute_model(
        procedure=Procedure.RRR,
        scaling_type=Scaling.COLUMNS,
        centering_type=None,
        # balancing=Balancing.BALANCE,
        refine_ritz=RefineRitz.MINIMIZE_RESIDUALS(100),
        rank_mode=RankMode.FORCE(force_rank),
        # rank_mode=RankMode.THRESH_NONZERO_FIRST,
    )
    # Print rank
    # print(f'Rank: {model.ritz_values_lmbl.shape[0]}')

    # Compute coefficients/amplitudes using all snapshots
    # eps = np.finfo(float).eps
    model.solve_coeffs_alpha(
        tolerance_level=4e10,
        observable_weights_scheme=WeightsScheme.UNIFORM,
    )

    # Extract row from predictions on lifted observables corresponding to last time delay
    predictions = model.predict_lifted(predict_horizon=predict_horizon).real[
        -obs_num, :
    ]
    # Threshold too low values to minimal grip_approx and too high values to maximum grip_approx
    predictions[predictions < threshold[0]] = threshold[0]
    predictions[predictions > threshold[1]] = threshold[1]

    return predictions.ravel()


def construct_grid_obs(
    time_delay_lifted: np.ndarray,
    grid_div: np.ndarray,
    keep_indices: list = None,
) -> tuple[np.ndarray, list]:
    # Function for creating gridded identity observables from time delay lifted data
    grid_obs = []
    num_delays = time_delay_lifted.shape[0]
    sparsify_indices = []
    for ind, bounds in enumerate(
        itertools.product(
            np.lib.stride_tricks.sliding_window_view(grid_div, 2), repeat=num_delays
        )
    ):
        # Skip if index not kept
        if keep_indices and (ind not in keep_indices):
            continue
        # Check all time delays against grid bounds, collapse using all and convert to float (0, 1)
        temp_obs = np.logical_and(
            #
            time_delay_lifted.T > np.array(bounds)[:, 0],
            time_delay_lifted.T <= np.array(bounds)[:, 1],
        )
        temp_obs = temp_obs.all(axis=1).astype(float)
        # Check if observable density is at least 0.1 % to keep it
        if (temp_obs.sum() / temp_obs.shape[0]) > 0.001:
            # if (temp_obs.sum() / temp_obs.shape[0]) > 0.000:
            grid_obs.append(temp_obs)
            sparsify_indices.append(ind)
        elif keep_indices and (ind in keep_indices):
            # If index is kept, but all zeros, add zeros to grid_obs
            grid_obs.append(temp_obs)
            sparsify_indices.append(ind)
    return np.array(grid_obs), sparsify_indices


def block_hankel_lift(
    observables: np.ndarray, delay: int, thin_hankel: int = 1
) -> np.ndarray:
    """Block hankel lift for observables - time delay embedding."""
    # HAs to be 2D array, convert if not
    if observables.ndim == 1:
        observables = observables.reshape(1, -1)
    obs_num = observables.shape[0]
    # Create hankel for each observable (list of arrays)
    # Function accepts first column and last row
    hankel_concat = [
        scipy.linalg.hankel(observables[i, :][: delay + 1], observables[i, :][delay:])
        for i in range(obs_num)
    ]
    # Create a block hankel from each hankel
    block_hankel = []
    # State + delays
    for row in range(delay + 1):
        # Observables
        for obs_n in range(obs_num):
            block_hankel.append(hankel_concat[obs_n][row])
    # Thin the block hankel but keep last snapshot
    block_hankel = np.array(block_hankel)[:, ::-1][:, ::thin_hankel][:, ::-1]
    return block_hankel


@dataclass(slots=True)
class StaticKoopmanInput:
    observables: np.ndarray
    lin_trans: np.ndarray


class ParameterOptimization:
    def __init__(
        self,
        measure_position: int,
        in_data_file: str,
        in_optimal_mask: str,
        output_folder: str,
        batch_wnd_coeff: float,
    ):
        self.measure_position = measure_position
        self.in_data_file = in_data_file
        self.in_optimal_mask = in_optimal_mask
        self.output_folder = Path(output_folder)
        self.batch_wnd_coeff: float = batch_wnd_coeff
        # Create folder if it doesn't exist
        self.output_folder.mkdir(parents=True, exist_ok=True)

    def init_optimization_problem(self):
        # Create problem instance (load data and prepare for optimization)
        # print(f"Measuring position: {self.measure_position}")
        self.prob = pg.problem(
            SignalProcessingParams(
                measure_position=self.measure_position,
                in_data_file=self.in_data_file,
            )
        )
        logging.error(self.prob)
        # Extract problem
        self.prob_extract = self.prob.extract(SignalProcessingParams)

        # Load and apply optimal mask to all measurements on single position
        opt_mask = pd.read_csv(self.in_optimal_mask)
        opt_mask_vals = opt_mask["mask"].to_numpy()
        self.proc_data = self.prob_extract.process_store_emg_corrs(opt_mask_vals)

    def run_grid_search(self):
        # Run the pipeline - sample, evaluate, analyze
        logging.error("Running grid search!")
        # print("Running grid search!")
        logging.error(f"Prediction window size: {self.batch_wnd_coeff}")
        # print(f"Prediction window size: {self.batch_wnd_coeff}")

        #### Grid search for optimal parameters for PyKMD prediction
        # Specify downsampling step - estimation and prediction
        dat_step = 8
        # Number of time delays - estimation
        num_delays = 60
        # Batch size as parameter
        batch_size = 496 // dat_step
        # Selection of time delays for grid observables
        grid_td_select = [1, num_delays // 2, num_delays]
        # Slice_to_train Koopman - all the data
        slice_train = np.s_[:]
        # Grid division for gridded identity observables
        grid_div = np.linspace(0, 1, 22) ** 1.8

        # Grid search for optimal parameters
        self.optim_results = {
            "position": [],
            "filename": [],
            "batch_wnd_coeff": [],
            "batch_smooth_coeff": [],
            "thin_step": [],
            "predict_horizon": [],
            "num_delays_predict": [],
            "force_rank": [],
            "wmape_pred": [],
        }

        # Iterate over all keys, and obtain wmape approximations and parameters
        for id_ind in range(len(self.proc_data["time_t"])):
            if id_ind not in [0, 1]:
                continue
            logging.error(f"Processing file: {self.proc_data['file_names'][id_ind]}")

            # Load data
            emg_proc = self.proc_data["emg_processed"][id_ind]
            grip = self.proc_data["grip"][id_ind]

            # Scale emg data to [0, 1] and grip to min-max range of the downsampled data
            transformer_emg = MinMaxScaler(feature_range=(0, 1))
            transformer_grip = MinMaxScaler(
                feature_range=(grip[::dat_step].min(), grip[::dat_step].max())
            )

            emg_proc_scaled = transformer_emg.fit_transform(
                emg_proc[::dat_step][slice_train].reshape(-1, 1),
                grip[::dat_step][slice_train].reshape(-1, 1),
            ).ravel()
            grip_trans = transformer_grip.fit_transform(
                grip[::dat_step][slice_train].reshape(-1, 1)
            ).ravel()

            emg_td_embed = blocked_hankel_lift(
                observables=emg_proc_scaled.reshape(1, -1), num_delays=num_delays
            )
            # Use only 2 - 3 state/time delays since number of combinations explodes
            grid_obs, grid_obs_ind = construct_grid_obs(
                emg_td_embed[grid_td_select, :], grid_div
            )

            grip_emg_lift = StaticKoopmanInput(
                observables=np.vstack((emg_td_embed, grid_obs)),
                lin_trans=np.vstack(
                    (
                        blocked_hankel_lift(grip_trans.reshape(1, -1), num_delays),
                        np.zeros_like(grid_obs),
                    )
                ),
            )
            # Static Koopman operator
            koopman_operator = grip_emg_lift.lin_trans @ scipy.linalg.pinv(
                grip_emg_lift.observables
            )

            for batch_smooth_coeff in np.linspace(1.2, 1.9, 8):
                if batch_smooth_coeff != 1.2:
                    continue
                for thin_step in range(3, 9):
                    # if thin_step != 3:
                    #     continue
                    predict_horizon = batch_size // thin_step
                    for num_delays_predict in range(4, 11):
                        # if num_delays_predict != 4:
                        #     continue
                        for force_rank in range(3, num_delays_predict + 1):
                            # if force_rank != 4:
                            #     continue
                            # Limit force_rank to 7
                            if force_rank > 7:
                                continue
                            # Position doesn't change
                            self.optim_results["position"].append(self.measure_position)
                            # Append batch window coefficient
                            self.optim_results["batch_wnd_coeff"].append(
                                self.batch_wnd_coeff
                            )
                            # Store filename
                            self.optim_results["filename"].append(
                                self.proc_data["file_names"][id_ind]
                            )
                            self.optim_results["batch_smooth_coeff"].append(
                                batch_smooth_coeff
                            )
                            self.optim_results["thin_step"].append(thin_step)
                            self.optim_results["predict_horizon"].append(
                                predict_horizon
                            )
                            self.optim_results["num_delays_predict"].append(
                                num_delays_predict
                            )
                            self.optim_results["force_rank"].append(force_rank)

                            # Compute grip approximation and prediction on bulk data
                            grip_blk_orig = []
                            grip_blk_approx = []
                            grip_blk_approx_smooth = []
                            grip_blk_approx_predict = []
                            try:
                                for em_blk, grp_blk in np.nditer(
                                    [emg_proc, grip],
                                    flags=["external_loop", "buffered"],
                                    buffersize=496,
                                ):
                                    if em_blk.size < 496 or grp_blk.size < 496:
                                        continue
                                    # Downsample EMG signal
                                    em_blk_dwns = em_blk[::dat_step]
                                    # Append downsampled grip
                                    grip_blk_orig.append(
                                        transformer_grip.transform(
                                            grp_blk[::dat_step].reshape(-1, 1)
                                        ).ravel()
                                    )
                                    emg_blk_trans = transformer_emg.transform(
                                        em_blk_dwns.reshape(-1, 1)
                                    ).ravel()

                                    emg_td_embed_blk = blocked_hankel_lift(
                                        emg_blk_trans.reshape(1, -1),
                                        num_delays=num_delays,
                                    )
                                    # Constuct grid observables, but keep indices from learned model on entire dataset
                                    grid_obs, _ = construct_grid_obs(
                                        emg_td_embed_blk[grid_td_select, :],
                                        grid_div,
                                        keep_indices=grid_obs_ind,
                                    )
                                    grip_blk_approx.append(
                                        koopman_operator
                                        @ np.vstack((emg_td_embed_blk, grid_obs))
                                    )
                                    grip_blk_approx[-1] = np.hstack(
                                        (
                                            grip_blk_approx[-1][0, :],
                                            grip_blk_approx[-1][1 : num_delays + 1, -1],
                                        )
                                    ).real
                                    # Threshold too low values to -1 N if below
                                    grip_blk_approx[-1][grip_blk_approx[-1] < -1] = -1
                                    # Stacked history - batch_smooth_coeff is always <= then 2, so last two batches will be enough
                                    grip_approx_hist = np.hstack(grip_blk_approx[-2:])
                                    # Smooth approximation
                                    # Fraction of data for smoothing
                                    if len(grip_approx_hist) > 1 * batch_size:
                                        smooth_frac = (
                                            batch_size * batch_smooth_coeff
                                        ) / len(grip_approx_hist)
                                    else:
                                        # Smooth using entire batch only for the first batch
                                        smooth_frac = 1
                                    # Smooth data using LOWESS
                                    lowess_smooth = sm.nonparametric.lowess(
                                        grip_approx_hist,
                                        range(len(grip_approx_hist)),
                                        frac=smooth_frac,
                                        it=0,
                                        is_sorted=True,
                                        missing="none",
                                        return_sorted=False,
                                    )
                                    # Append only last batch size to smoothed data
                                    grip_blk_approx_smooth.append(
                                        lowess_smooth[-batch_size:]
                                    )
                                    # Threshold too low values to -1 N if below
                                    grip_blk_approx_smooth[-1][
                                        grip_blk_approx_smooth[-1] < -1
                                    ] = -1
                                    grip_blk_approx_predict.append(
                                        predict_grip_pykmd(
                                            # Take batch size modified by batch_wnd_coeff for prediction
                                            grip_approx=np.hstack(
                                                grip_blk_approx_smooth[-2:]
                                            )[
                                                -int(
                                                    batch_size * self.batch_wnd_coeff
                                                ) :
                                            ].reshape(
                                                1, -1
                                            ),
                                            num_delays=num_delays_predict,
                                            predict_horizon=predict_horizon,
                                            thinning_step=thin_step,
                                            # Threshold for predictions, min and max from initial experiment
                                            threshold=[
                                                grip[::dat_step].min(),
                                                grip[::dat_step].max(),
                                            ],
                                            force_rank=force_rank,
                                        )
                                    )
                            except ValueError:
                                # If there is an error, append nan to wmape prediction
                                self.optim_results["wmape_pred"].append(np.nan)
                                continue

                            num_batches = len(grip_blk_orig)
                            grip_blk_orig = np.hstack(grip_blk_orig)
                            grip_blk_approx_predict = np.hstack(grip_blk_approx_predict)

                            plot_indices = np.arange(len(grip_blk_orig))
                            # Plot indices for predictions
                            plot_predict_indices = (
                                np.arange(1, num_batches + 1) * batch_size - 1
                            )
                            plot_predict_indices = (
                                plot_predict_indices.reshape(-1, 1)
                                + thin_step * np.arange(1, predict_horizon + 1)
                            ).ravel()
                            wmape_pred = (
                                np.abs(
                                    grip_blk_orig[
                                        plot_predict_indices[
                                            plot_predict_indices <= plot_indices.max()
                                        ]
                                    ]
                                    - grip_blk_approx_predict[
                                        plot_predict_indices <= plot_indices.max()
                                    ]
                                ).sum()
                                / np.abs(
                                    grip_blk_orig[
                                        plot_predict_indices[
                                            plot_predict_indices <= plot_indices.max()
                                        ]
                                    ]
                                ).sum()
                                * 100
                            )
                            # Store wmape prediction
                            self.optim_results["wmape_pred"].append(wmape_pred)

    def store_results(self):
        # Convert results to pandas
        optim_results_df = pd.DataFrame.from_dict(self.optim_results)

        # Export results to csv
        logging.error("Storing grid search results as a csv!")
        file_path = (
            self.output_folder
            / f"optim_grid_pred_wnd{self.batch_wnd_coeff}_pos{self.measure_position}.csv"
        )
        optim_results_df.to_csv(file_path, index_label="index")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.WARNING,
        format="%(asctime)s|%(name)s|%(levelname)s|%(message)s",
    )
    # Add measure_position: int, in_data_file: str, in_optimal_mask: str,
    # output_folder: str, to argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-pos",
        "--measure-position",
        type=int,
        help="Measurement position",
        required=True,
    )
    parser.add_argument(
        "-batch",
        "--batch-wnd-coeff",
        type=float,
        help="Window size modifying coefficient for prediction",
        required=True,
    )
    parser.add_argument(
        "-in",
        "--in-data-file",
        type=str,
        help="Input data file with measurements for sensitivity analysis",
        required=True,
    )
    parser.add_argument(
        "-mask",
        "--in-optimal-mask",
        type=str,
        help="Input optimal mask csv file",
        required=True,
    )
    parser.add_argument(
        "-out",
        "--output-folder",
        type=str,
        help="Output folder for sensitivity analysis results",
        required=True,
    )
    args = parser.parse_args()
    parameter_optimization = ParameterOptimization(
        measure_position=args.measure_position,
        batch_wnd_coeff=args.batch_wnd_coeff,
        in_data_file=args.in_data_file,
        in_optimal_mask=args.in_optimal_mask,
        output_folder=args.output_folder,
    )
    parameter_optimization.init_optimization_problem()
    parameter_optimization.run_grid_search()
    parameter_optimization.store_results()
