#!/usr/bin/env python3

import subprocess
import multiprocessing as mp
import os
import argparse
import numpy as np

ap = argparse.ArgumentParser(
description="Run isotherm test for validation of GCMC code")
ap.add_argument("fastmc_exe", metavar="x", type=str,
                help="Path to executable for fastmc.")
ap.add_argument("num_cores", metavar="n", type=int,
                help="Number of cores to use for long test.")
ap.add_argument("test_type", metavar="t", type=str,
                help="Type of test, acceptable inputs are 'long' or 'short'." +
                     "Long test runs full isotherm (parallelizable up to 10 cores," +
                     "8 mins per point)." +
                     "Short test runs a single isotherm point (single core ~8min).")
args = ap.parse_args()
N = args.num_cores
exe = args.fastmc_exe

# Values to compare to for isotherm test
# (these are avg number of guests in the following
# format: P: [CO2, N2]
reference_dict = {
        0.01: [4.58, 0.474],
        0.049997: [19.82, 1.015],
        0.099988: [32.10, 1.283],
        0.149971: [40.67, 1.435],
        0.199949: [46.63, 1.536],
        0.399798: [58.63, 1.714],
        0.599547: [65.17, 1.767],
        0.799195: [68.70, 1.850],
        0.998743: [72.70, 1.773],
        1.198191: [74.04, 1.908],
        }

reference_energy = -11190.606165

script_path, script_name = os.path.split(os.path.realpath(__file__))
test_path = '{}/test'.format(script_path)
dirs = [d for d in os.listdir(test_path)
        if os.path.isdir('{}/{}'.format(test_path, d))]

def run_bash(cmd):
    p = subprocess.Popen([cmd],
                            shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    out, err = p.communicate()
    return out.decode('ascii')

def main(dir):
    data_dict = {}
    P = 0
    avg_guests = []
    std_dev = []
    full_dir_path = '{}/{}'.format(test_path, dir)
    os.chdir(full_dir_path)
    run_bash(exe)

    with open('{}/OUTPUT'.format(full_dir_path), 'r') as f:
        for line in f:
            if 'Electrostatic energy' in line:
                E = float(line.strip().split()[3])
            if 'average number of guests' in line:
                avg_guests.append(float(line.strip().split()[4]))
            if 'standard error' in line:
                std_dev.append(float(line.strip().split()[2]))
        f.close()

    std_dev = [std_dev[0], std_dev[2]]

    with open('{}/CONTROL'.format(full_dir_path), 'r') as f:
        for line in f:
            if 'temperature' in line:
                T = float(line.strip().split()[1])
            if 'pressure' in line:
                P += float(line.strip().split()[2])
        f.close()

    min_uptake_CO2 = (np.asarray(avg_guests) - np.asarray(std_dev))[0]
    min_uptake_N2 = (np.asarray(avg_guests) - np.asarray(std_dev))[1]
    max_uptake_N2 = (np.asarray(avg_guests) + np.asarray(std_dev))[1]
    max_uptake_CO2 = (np.asarray(avg_guests) + np.asarray(std_dev))[0]

    return [E, min_uptake_CO2, max_uptake_CO2, min_uptake_N2, max_uptake_N2, P]


if __name__ == '__main__':

    if args.test_type == 'long':
        overall_dict = {}

        pool = mp.Pool(processes=N)
        for results in pool.imap_unordered(main, dirs):
            overall_dict[results.pop()] = results
        overall_dict = dict(sorted(overall_dict.items()))

        # Give a little extra tolerance (small # steps) 
        for key in overall_dict.keys():
            min_uptake_CO2 = overall_dict[key][1] * 0.8
            max_uptake_CO2 = overall_dict[key][2] * 1.2
            min_uptake_N2 = overall_dict[key][3] * 0.8
            max_uptake_N2 = overall_dict[key][4] * 1.2
            ref_uptake_CO2 = reference_dict[key][0]
            ref_uptake_N2 = reference_dict[key][1]
            energy = overall_dict[key][0]

            if (min_uptake_CO2 <= ref_uptake_CO2) <= max_uptake_CO2 and (min_uptake_N2 <= ref_uptake_N2 <= max_uptake_N2) and (energy == reference_energy):
                pass_test = True
            else:
                pass_test = False

            if pass_test:
                print("Pressure: {}     PASS".format(key))
            else:
                print("Pressure: {}     FAIL".format(key))

    elif args.test_type == 'short':

        # Run for single isotherm point (P = 0.20 bar)

        result = main(dirs[4])
        P = result.pop()

	# Give a little extra tolerance (small # steps)
        min_uptake_CO2 = result[1] * 0.9
        max_uptake_CO2 = result[2] * 1.1
        min_uptake_N2 = result[3] * 0.9
        max_uptake_N2 = result[4] * 1.1
        ref_uptake_CO2 = reference_dict[P][0]
        ref_uptake_N2 = reference_dict[P][1]
        energy = result[0]
        
        if (min_uptake_CO2 <= ref_uptake_CO2) <= max_uptake_CO2 and (min_uptake_N2 <= ref_uptake_N2 <= max_uptake_N2) and (energy == reference_energy):
            pass_test = True
        else:
            pass_test = False

        if pass_test:
            print("Pressure: {}     PASS".format(P))
        else:
            print("Pressure: {}     FAIL".format(P))

    else:
        print("Unknown input for test type... please specify either 'long' or 'short'")
