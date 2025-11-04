#!/usr/bin/env python3

# Usage:
#   ./multiRuns.py --test_dir DummyTest --input_file Spray-Conv.inp --pele_exec PeleC-Spray-Conv --run_cmd "mpiexec -np 4"

# Input:
#   * input_file: name of the input file, assumes it ends in "inp" if not given
#   * test_dir: name of test directory, test directory must contain the executable, input file, and gridfiles
#   * run_cmd: command to run test cases with, like 'mpiexec -np 4' or just './'
#   * pele_exec: PeleC executable, assumed to be in test_dir
#   * max_error: Maximum allowable error when using FCompare
#   * fcomp_exec: FCompare executable

# "Internal" user input
#   * bsizes: box sizes to test
#   * pltfiles: names of plot files to do compare on
#   * gridfiles: names of gridfiles used


import sys
import os
import shutil
import argparse

USAGE = """
    A script for doing CI with PeleC coupled with sprays
"""


def multiRun(args):
    print(" Scripted runs for spray CI ")
    # User data
    bsizes = [32, 64]
    iter_1 = 5
    iter_2 = 9

    # Get the PeleC exec
    run_cmd = args.run_cmd
    test_dir = args.test_dir
    executable = args.pele_exec
    if executable == "None":
        for f in os.listdir(test_dir):
            if f.startswith("PeleC") and f.endswith(".ex"):
                executable = f
    if not os.path.exists(executable):
        errorStatement = "Pele executable not found"
        raise ValueError(errorStatement)

    # Check for the input file: first input.* is default
    if args.input_file == "None":
        for f in os.listdir(test_dir):
            if f.endswith("inp"):
                args.input_file = f
                break
    args.input_file = test_dir + "/" + args.input_file
    if not os.path.exists(args.input_file):
        errorStatement = args.input_file + " file not found"
        raise ValueError(errorStatement)
    # Delete previous test directory if it exists
    test_name = "testfiles"
    for f in os.listdir(test_dir):
        if f == test_name:
            os.system("rm -r {}".format(test_dir + "/" + test_name))
    test_name = test_dir + "/" + test_name
    os.system("mkdir -p {}".format(test_name))
    pltfiles = []
    for box in bsizes:
        print(" Running {}x{} box sizes".format(box, box))
        runtime_params = get_runtime_params(box, args)
        first_params = "amr.check_file={}/{}_chk ".format(test_name, box)
        first_params += "max_step = {} ".format(iter_2)
        first_params += "amr.check_int = {} ".format(iter_1)
        first_params += "amr.checkpoint_files_output = 1 "
        first_params += "amr.file_name_digits = 5 "
        first_params += "amr.plot_file = {}/{}_1_plt ".format(test_name, box)
        first_params += "amrex.abort_on_unused_inputs=1 "
        # Run first test to iter_2 and save a checkpoint at iter_1
        assert(os.system(
            "{} {} {} {} {}".format(
                run_cmd, executable, args.input_file, runtime_params, first_params
            )
        ) == 0)
        second_params = "max_step = {} ".format(iter_2)
        second_params += "amr.restart = {}/{}_chk0000{} ".format(test_name, box, iter_1)
        second_params += "amr.checkpoint_files_output = 0 "
        second_params += "amr.plot_file = {}/{}_2_plt ".format(test_name, box)
        # Restart a case from checkpoint at iter_1
        assert(os.system(
            "{} {} {} {} {}".format(
                run_cmd, executable, args.input_file, runtime_params, second_params
            )
        ) == 0)
        # Now delete the checkpoint files and starting plot files
        print("Deleting unused plot and checkpoint files")
        os.system("rm -r {}/*plt00000".format(test_name, box))
        os.system("rm -r {}/{}_chk*".format(test_name, box))
        pltname1 = "{}_1_plt0000{}".format(box, iter_2)
        pltname2 = "{}_2_plt0000{}".format(box, iter_2)
        pltfiles.append(pltname1)
        pltfiles.append(pltname2)
    return pltfiles


def get_runtime_params(box, args):
    runtime_params = "amr.plot_int=100 "
    runtime_params += "amr.initial_grid_file={}/gridfile_{}_1.dat ".format(
        args.test_dir, box
    )
    runtime_params += "amr.regrid_file={}/gridfile_{}_2.dat ".format(args.test_dir, box)
    runtime_params += "amr.max_grid_size={} ".format(box)
    runtime_params += "amr.blocking_factor={} ".format(box)
    return runtime_params


def compare(pltfiles, args):

    # User data
    vars = ["density", "rho_NC10H22", "xmom", "Temp", "rho_E", "spray_vol"]

    # Get a local copy of post-processing executable
    test_dir = args.test_dir
    fcomp_exec = args.fcomp_exec
    # Check if fcompare executable is in current directory
    if not os.path.exists(fcomp_exec):
        errorStatement = fcomp_exec + " executable not found"
        raise ValueError(errorStatement)
    test_name = test_dir + "/testfiles"
    # We have 3 comparisons to make,
    # 32_1 to 64_1, 32_1 to 64_2, and 64_1 to 64_2
    numcomps = 3
    comp1 = [0, 0, 2]
    comp2 = [2, 3, 3]
    maxerror = args.max_error
    errfail = 0
    for comp in range(numcomps):
        indx1 = comp1[comp]
        indx2 = comp2[comp]
        plt1 = test_name + "/" + pltfiles[indx1]
        plt2 = test_name + "/" + pltfiles[indx2]
        print("Comparing the L2 norm of the error between " + plt1 + " and " + plt2)
        outfile = "{}/error_{}.analysis.out".format(test_name, comp)
        if not os.path.exists(plt1):
            errorStatement = plt1 + " file not found"
            raise ValueError(errorStatement)
        if not os.path.exists(plt2):
            errorStatement = plt2 + " file not found"
            raise ValueError(errorStatement)
        os.system("{} -n 2 -a {} {} > {}".format(fcomp_exec, plt1, plt2, outfile))
        # Extract errors on each variable
        curlev = 0
        varcount = 0
        with open(outfile) as fp:
            for i, line in enumerate(fp):
                if i >= 2:
                    var = line.split()[0]
                    if var == "level":
                        curlev = int(line.split()[2])
                    for v in range(len(vars)):
                        if var == vars[v]:
                            if curlev == 0:
                                varcount += 1
                            cer = float(line.split()[2])
                            if cer > maxerror:
                                errfail += 1
                                perr = "Error in {} on level {}: {} > {}".format(
                                    var, curlev, cer, maxerror
                                )
                                print(perr)
        if varcount != len(vars):
            error = "Not all variables were found in plot file"
            raise ValueError(error)
    if errfail > 0:
        errorStatement = "Failure detected in spray test"
        raise ValueError(errorStatement)
    else:
        print("Spray tests passed")


def parse_args(arg_string=None):
    parser = argparse.ArgumentParser(description=USAGE)

    parser.add_argument(
        "--test_dir",
        type=str,
        default=".",
        help="directory where executable, gridfiles, and input files.",
    )

    parser.add_argument(
        "--input_file",
        type=str,
        default="None",
        metavar="input.2d",
        help="input file name. Default = first inputs.* in current directory.",
    )

    parser.add_argument(
        "--pele_exec", type=str, default="None", help="PeleC executable."
    )

    parser.add_argument(
        "--run_cmd", type=str, default="", help="MPI or serial run command."
    )

    parser.add_argument(
        "--max_error", type=float, default=1.0e-12, help="max error for tests."
    )

    parser.add_argument(
        "--fcomp_exec", type=str, default="None", help="name fcompare executable."
    )

    if not arg_string is None:
        args, unknown = parser.parse_known_args(arg_string)
    else:
        args, unknown = parser.parse_known_args()

    return args


if __name__ == "__main__":
    args = parse_args(arg_string=sys.argv[1:])
    pltfiles = multiRun(args)
    compare(pltfiles, args)
