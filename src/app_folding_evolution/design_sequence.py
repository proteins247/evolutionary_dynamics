#!/usr/bin/env python
"""
design_sequence.py

Design sequence using Z-score method

author : Victor Zhao

"""
import sys
import os
import argparse
import pickle as pkl

import numpy as np
import pandas as pd

sys.path.append("/n/shakfs1/users/vzhao/project__cotranslational_fold/simulations/trajectory_analysis")
import lattice_protein as lp


CONFORMATION_FILE = os.environ["FOLDEVO_SHARE"] + "/LPforms/LP10000_latpack-form.txt"


def read_conformations():
    conformations = []
    with open(CONFORMATION_FILE) as f:
        for line in f:
            conformations.append(line.strip())
    return conformations


def design_sequence(input_sequence, points, potential, target=-45,
                    max_iterations=50000, temperature=0.1,
                    mutation_range=None):
    sequence = list(input_sequence)  # make a copy
    energy, contact_count = lp.energy(points, sequence, potential)
    z_score = lp.z_score(points, sequence, potential, energy, contact_count)
    n_iterations = 0
    while z_score > target and n_iterations < max_iterations:
        new_sequence = lp.mutate_number_sequence(
            sequence, length_range=mutation_range)
        new_energy = lp.energy(points, new_sequence, potential)[0]
        new_z_score = lp.z_score(points, new_sequence, potential,
                                 new_energy, contact_count)

        if new_z_score < z_score:
            sequence = new_sequence
            z_score = new_z_score
            energy = new_energy
        elif np.random.rand() < np.exp(
                -(new_z_score - z_score) / temperature):
            sequence = new_sequence
            z_score = new_z_score
            energy = new_energy
        n_iterations += 1
    if n_iterations >= max_iterations:
        print("--> Reached iteration limit", max_iterations)
    return sequence, energy, z_score


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("conformation", type=int)
    parser.add_argument("--temperature", type=float, default=0.1)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--target-score", type=float, default=-45)
    parser.add_argument("--max-iterations", type=int, default=30000)
    return parser.parse_args()


def main(args):
    conformations = read_conformations()
    conformation = conformations[args.conformation]
    potential = lp.read_energy_file(
        "/n/home00/vzhao/code/evolutionary_dynamics/share/MJ96/energy__MJ85-tableVI.dat")
    points = lp.conformation_to_points(conformation)

    print("Designing conformation %d: %s" % (args.conformation, conformation))

    contact_order = lp.contact_order(conformation)
    print("Contact order:", contact_order)

    if args.seed is not None:
        np.random.seed(args.seed)

    # Random sequence
    length = len(conformation) + 1
    sequence = lp.random_sequence(length)
    energy, contact_count = lp.energy(points, sequence, potential)
    z_score = lp.z_score(points, sequence, potential, energy, contact_count)

    print("Starting sequence:", lp.number_sequence_to_letters(sequence))
    print("Starting energy:", energy)
    print("Starting z-score:", z_score)

    sequence, energy, z_score = design_sequence(
        sequence, points, potential, temperature=args.temperature,
        target=args.target_score, max_iterations=args.max_iterations)

    print("Ending sequence:", lp.number_sequence_to_letters(sequence))
    print("Ending energy:", energy)
    print("Ending z-score:", z_score)

    outfile = "%d.sourceme" % (args.conformation)
    print("Writing to", outfile)
    print()
    with open(outfile, 'w') as f:
        f.write("CONFORMATIONNO=%d\n" % args.conformation)
        f.write("CONTACTORDER=%f\n" % contact_order)
        f.write("DESIGNTEMP=%.2f\n" % args.temperature)
        f.write("AASEQ=%s\n" % lp.number_sequence_to_letters(sequence))
        f.write("NATIVEENERGY=%.2f\n" % energy)
        f.write("CONFORMATION=%s\n" % conformation)
        f.write("ZSCORE=%.2f\n" % z_score)
    if z_score > args.target_score:
        return 1
    else:
        return 0


if __name__ == "__main__":
    sys.exit(main(parse_args()))
