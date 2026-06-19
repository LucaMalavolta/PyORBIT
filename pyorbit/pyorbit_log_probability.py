from __future__ import print_function

import argparse
import contextlib
import io
import json
import os
import sys

import pyorbit
from pyorbit.subroutines.common import np
from pyorbit.subroutines.input_parser import yaml_parser
import pyorbit.subroutines.results_analysis as results_analysis
from pyorbit.samplers.pyorbit_log_probability import (
    prepare_log_probability_model,
    pyorbit_log_probability,
)

__all__ = ["pyorbit_log_probability_cli"]


def _print_banner():
    print()
    print(":::---::: PyORBIT v{0} :::---:::".format(pyorbit.__version__))
    print()
    print("Python version in use:")
    print(sys.version)


def _read_mapping_value(data, key):
    if key is None:
        if isinstance(data, dict) and "theta" in data:
            return data["theta"]
        return data

    value = data
    for key_part in key.split("."):
        value = value[key_part]
    return value


def _load_theta_file(theta_file, key=None):
    extension = os.path.splitext(theta_file)[1].lower()

    if extension == ".npy":
        return np.load(theta_file)

    if extension == ".npz":
        theta_data = np.load(theta_file)
        if key is None:
            if "theta" in theta_data.files:
                key = "theta"
            elif len(theta_data.files) == 1:
                key = theta_data.files[0]
            else:
                raise ValueError(
                    "NPZ theta files with multiple arrays require --theta-key"
                )
        return theta_data[key]

    if extension in [".json"]:
        with open(theta_file, "r") as theta_handle:
            return _read_mapping_value(json.load(theta_handle), key)

    if extension in [".yaml", ".yml"]:
        try:
            import yaml
        except ImportError:
            raise ImportError("pyyaml is required to read YAML theta files")

        with open(theta_file, "r") as theta_handle:
            return _read_mapping_value(yaml.safe_load(theta_handle), key)

    return np.loadtxt(theta_file)


def _print_parameter_order(theta_dictionary):
    print()
    print("Theta parameters, in input order:")
    for theta_name, theta_i in sorted(theta_dictionary.items(), key=lambda item: item[1]):
        print("  {0:4d}  {1}".format(theta_i, theta_name))


def _print_theta_yaml_template(theta_dictionary):
    for theta_name, theta_i in sorted(theta_dictionary.items(), key=lambda item: item[1]):
        print("{0}: ".format(json.dumps(str(theta_name))))


def _print_result(result):
    print()
    print("Model source: ", result["model_source"])
    print("Dimensions = ", result["theta"].shape[-1])
    print()

    if np.ndim(result["log_probability"]) == 0:
        print("log_prior       = {0:.16e}".format(result["log_prior"]))
        print("log_likelihood  = {0:.16e}".format(result["log_likelihood"]))
        print("log_probability = {0:.16e}".format(result["log_probability"]))
        return

    print("index log_prior log_likelihood log_probability")
    for ii in range(0, len(result["log_probability"])):
        print(
            "{0:d} {1:.16e} {2:.16e} {3:.16e}".format(
                ii,
                result["log_prior"][ii],
                result["log_likelihood"][ii],
                result["log_probability"][ii],
            )
        )


def pyorbit_log_probability_cli():
    parser = argparse.ArgumentParser(
        prog="PyORBIT_LogProbability.py",
        description="Compute log_probability = log_prior + log_likelihood for theta",
    )
    parser.add_argument("config_file", type=str, help="config file")
    parser.add_argument(
        "theta",
        type=float,
        nargs="*",
        help="theta values, in sampler order",
    )
    parser.add_argument(
        "-t",
        "--theta-file",
        type=str,
        default=None,
        help="file containing theta values",
    )
    parser.add_argument(
        "--theta-key",
        type=str,
        default=None,
        help="key for JSON/YAML/NPZ theta files; nested keys use dot notation",
    )
    parser.add_argument(
        "--reload-from",
        type=str,
        choices=["auto", "none", "emcee", "pyde"],
        default="auto",
        help="saved state used to recover theta ordering and bounds",
    )
    parser.add_argument(
        "--print-theta-dictionary",
        action="store_true",
        help="print a YAML theta template in sampler order and exit",
    )
    parser.add_argument(
        "--print-parameters",
        action="store_true",
        help="print theta parameter names in the order they must be provided and exit",
    )

    args = parser.parse_args()

    config_in = yaml_parser(args.config_file)

    if args.print_parameters and (args.theta_file is not None or len(args.theta) > 0):
        parser.error("--print-parameters does not accept theta values")

    if args.print_theta_dictionary and (args.theta_file is not None or len(args.theta) > 0):
        parser.error("--print-theta-dictionary does not accept theta values")

    if args.print_parameters or args.print_theta_dictionary:
        if args.print_theta_dictionary:
            stdout_stream = io.StringIO()
            stderr_stream = io.StringIO()
        else:
            stdout_stream = sys.stdout
            stderr_stream = sys.stderr

        with contextlib.redirect_stdout(stdout_stream), contextlib.redirect_stderr(stderr_stream):
            mc = prepare_log_probability_model(
                config_in,
                reload_from=args.reload_from,
            )
        theta_dictionary = results_analysis.get_theta_dictionary(mc)
        if args.print_theta_dictionary:
            _print_theta_yaml_template(theta_dictionary)
            return

        print()
        print("Dimensions = ", mc.ndim)
        _print_parameter_order(theta_dictionary)
        return

    if args.theta_file is not None and len(args.theta) > 0:
        parser.error("provide theta either as values or with --theta-file, not both")

    if args.theta_file is None and len(args.theta) == 0:
        parser.error("provide theta values or --theta-file")

    theta = (
        _load_theta_file(args.theta_file, args.theta_key)
        if args.theta_file is not None
        else args.theta
    )

    result = pyorbit_log_probability(
        config_in,
        theta,
        reload_from=args.reload_from,
        return_output=False,
    )

    _print_banner()
    _print_result(result)


if __name__ == "__main__":
    pyorbit_log_probability_cli()
