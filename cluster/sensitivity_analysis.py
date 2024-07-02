# /usr/bin/env python3.10
# -*- coding: utf-8 -*-
# Problem definition and a pipeline for sensitivity analysis
import argparse
import logging
from pathlib import Path

import numpy as np
import pygmo as pg
from SALib import ProblemSpec

# IMPORTANT: add grip-emg-optimize to PYTHONPATH
from optim.problem import SignalProcessingParams


class SensitivityAnalysis:
    def __init__(
        self,
        measure_position: int,
        groups: bool,
        in_data_file: str,
        output_folder: str,
        sample_size: int,
        num_resamples: int,
        nprocs: int,
    ):
        self.measure_position = measure_position
        self.groups = groups
        self.in_data_file = in_data_file
        self.output_folder = Path(output_folder)
        self.sample_size = sample_size
        self.num_resamples = num_resamples
        self.nprocs = nprocs
        # Create folder if it doesn't exist
        self.output_folder.mkdir(parents=True, exist_ok=True)

    def init_optimization_problem(self):
        # Create problem instance (load data and prepare for optimization/sensitivity)
        # TODO: problem parameters: bounds set inside the problem class
        self.prob = pg.problem(
            SignalProcessingParams(
                measure_position=self.measure_position,
                in_data_file=self.in_data_file,
            )
        )
        print(self.prob)
        # Extract problem
        self.prob_extract = self.prob.extract(SignalProcessingParams)

        # Decision vector parameter names using FFT frequencies + decay factor and
        # smoothing window size
        self.param_names = list(
            np.round(self.prob_extract.emg_fft_freq[1:], decimals=2)
        )
        self.param_names = ["f" + str(i) for i in self.param_names] + [
            "decay",
            "window",
        ]
        self.param_names = [i.replace(".", "_") for i in self.param_names]

    # Wrapper for evaluation function
    def corr_eval(self, X: np.ndarray) -> np.ndarray:
        # import numpy as np

        # Iterate over the sample rows (decision vectors) and evaluate the fitness function
        results = np.fromiter(
            iter=(self.prob.fitness(i) for i in X), dtype=np.float64, count=X.shape[0]
        )
        return 1.0 - results.ravel()

    def define_sensitivity_problem(self):
        # Define the sensitivity problem model
        # TODO hdmr variables: maxiter: int, m: int, K: int, R: int, lambdax: float
        problem_spec_dict = {
            "names": self.param_names,
            # Get bounds from the problem instance, just reframe in [low, high] format
            "bounds": [[low, up] for low, up in zip(*self.prob_extract.get_bounds())],
            "outputs": ["corr"],
        }

        if self.groups:
            # Add group information
            problem_spec_dict["groups"] = ["freqs"] * (len(self.param_names) - 2) + [
                "decay",
                "window",
            ]

        self.sensitivity_problem = ProblemSpec(problem_spec_dict)
        print(repr(self.sensitivity_problem))

    def run_sensitivity_analysis(self):
        # Run the pipeline - sample, evaluate, analyze
        logging.info("Running sensitivity analysis!")
        logging.info(f"Sample size: {self.sample_size}")
        logging.info(f"Number of parallel processors: {self.nprocs}")
        logging.info(f"Number of resamples for CI: {self.num_resamples}")
        (
            self.sensitivity_problem.sample_latin(N=self.sample_size)
            # self.sensitivity_problem.sample_saltelli(
            #     N=self.sample_size, calc_second_order=False)
            .evaluate(self.corr_eval, nprocs=self.nprocs)
            # .analyze_hdmr(
            #     nprocs=self.nprocs, maxorder=2, maxiter=200, m=2, K=20, print_to_console=True
            # )
            .analyze_rbd_fast(
                nprocs=self.nprocs,
                num_resamples=self.num_resamples,
                M=10,
                conf_level=0.95,
                print_to_console=True,
            )
            # .analyze_sobol(
            #     nprocs=self.nprocs,
            #     num_resamples=self.num_resamples,
            #     calc_second_order=False,
            #     conf_level=0.95,
            #     print_to_console=True,
            # )
        )

    def store_results(self):
        # Convert results to pandas
        sensitivity_list = self.sensitivity_problem.to_df()

        # For RBD FAST, list is a df
        sensitivity_df = sensitivity_list

        # For Sobol, join all dataframes
        # sensitivity_df = pd.DataFrame(sensitivity_list[0])
        # Join the dataframes by index passing a list of df
        # sensitivity_df = sensitivity_df.join(sensitivity_list[1:])

        # Export results to csv
        logging.info("Storing sensitivity analysis results as a csv!")
        file_path = (
            self.output_folder
            / f"sensitivity_groups{self.groups}_rbd_fast_samples{self.sample_size}_resamples{self.num_resamples}_pos{self.measure_position}_narrow_bounds15.csv"
        )
        sensitivity_df.to_csv(file_path, index_label="vars")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s|%(name)s|%(levelname)s|%(message)s",
    )
    # Add measure_position: int, groups: bool, in_data_file: str, output_folder: str,
    # sample_size: int, 4 num_resamples: int, nprocs: int to argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-pos",
        "--measure-position",
        type=int,
        help="Measurement position",
        required=True,
    )
    parser.add_argument(
        "-gr",
        "--groups",
        help="Group parameter in sensitivity analysis",
        action="store_true",
    )
    parser.add_argument(
        "-in",
        "--in-data-file",
        type=str,
        help="Input data file with measurements for sensitivity analysis",
        required=True,
    )
    parser.add_argument(
        "-out",
        "--output-folder",
        type=str,
        help="Output folder for sensitivity analysis results",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--sample-size",
        type=int,
        help="Sample size to generate for sensitivity analysis",
        required=True,
    )
    parser.add_argument(
        "-res",
        "--num-resamples",
        type=int,
        help="Number of resamples for confidence interval",
        required=True,
    )
    parser.add_argument(
        "-npr",
        "--nprocs",
        type=int,
        help="Number of parallel processors to use",
        required=True,
    )
    args = parser.parse_args()
    sensitivity_analysis = SensitivityAnalysis(
        measure_position=args.measure_position,
        groups=args.groups,
        in_data_file=args.in_data_file,
        output_folder=args.output_folder,
        sample_size=args.sample_size,
        num_resamples=args.num_resamples,
        nprocs=args.nprocs,
    )
    sensitivity_analysis.init_optimization_problem()
    sensitivity_analysis.define_sensitivity_problem()
    sensitivity_analysis.run_sensitivity_analysis()
    sensitivity_analysis.store_results()
