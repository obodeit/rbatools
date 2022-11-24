#!/usr/bin/env python3
"""Maximise growth rate and store results as sbtab"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import argparse
import sys

# package imports
from rbatools.rba_session import SessionRBA


def main():
    parser = argparse.ArgumentParser(description='Maximise growth rate and store results as sbtab')
    parser.add_argument('model_dir', metavar='model-dir', type=str,
                        help='Directory of RBA-model')
    parser.add_argument('--lp-solver', type=str, default="glpk",
                        help=(
                            'LP solver ("cplex", "glpk";'
                            'default: "glpk"'
                        ))
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Directory to save the SBtab file')

    args = parser.parse_args()

    Session=SessionRBA(xml_dir=args.model_dir,lp_solver=args.lp_solver)
    mumax=Session.find_max_growth_rate()
    print("Optimal growth-rate: {}".format(mumax))

    Session.record_results(run_name="run")
    Session.record_parameters(run_name="run")
    Session.write_results(session_name="rbatools_session")

    if args.output_dir is not None:
        output_dir=args.output_dir+"/"
    else:
        output_dir=""
    Session.SimulationData.export_sbtab(filename=output_dir+"SimulationResults_SBtab")
    print("Results written to: {}".format(output_dir+"SimulationResults_SBtab.tsv"))
    Session.SimulationParameters.export_sbtab(filename=output_dir+"SimulationParameters_SBtab")
    print("Paramters written to: {}".format(output_dir+"SimulationParameters_SBtab.tsv"))

if __name__ == '__main__':
    main()
