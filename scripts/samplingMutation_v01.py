import fnmatch
from multiprocessing.pool import Pool
import os
import argparse
import logging
from functools import lru_cache
import pandas as pd
import numpy as np
from collections import defaultdict

# script to generate randomised mutations
# written by Jordi Deu-Pons

# RANDOMIZATIONS
RANDOMIZATION = 1000

# Signature column names
REF = 'Signature_reference'
ALT = 'Signature_alternate'
PROBABILITY = 'Probability'

# Data column names
POS = 'pos'
TRI = 'triID'
COUNT = 'count'

@lru_cache(maxsize=1)
def tri_dictionary():
    tri_id = 1
    tri_dict = {}
    atgc = 'ATGC'
    for a in atgc:
        for b in atgc:
            for c in atgc:
                tri_dict[tri_id] = a + b + c
                tri_id += 1
    return tri_dict


def simulation(data_file, signature):

    # Load data
    #data = pd.read_csv(data_file, sep='\t', compression='gzip', header=None, names=["chr", "start", "end", POS, TRI, COUNT])
    print("reading the data in simulation step")
    data = pd.read_csv(data_file, sep='\t', header=None, names=["chr", "start", "end", POS, TRI, COUNT])
    mutation_column = COUNT
    print("finished reading the data in simulation step")
    tri_idmap = tri_dictionary()

    # Count observed mutations
    mutation_count = data[mutation_column].sum()
    mutation_real_pos = (data[data[mutation_column] > 0])[POS]

    if mutation_count == 0:
        return mutation_count, [], [], None

    # Simulate random mutations
    pos_vector = []
    prb_vector = []
    for ix, row in data.iterrows():
        pos_vector += [row[POS]]*3
        prb_vector += (signature.loc(axis=0)[tri_idmap[row[TRI]], ])[PROBABILITY].tolist()
        print("mid simulation step")

    # Normalize probabilities
    prb_vector = np.array(prb_vector)
    prb_vector = prb_vector / prb_vector.sum()

    mutation_rand_pos = np.random.choice(pos_vector, size=mutation_count*RANDOMIZATION, replace=True, p=prb_vector)

    # Full probabilities array
    probs_dict = defaultdict(list)
    for i, p in enumerate(pos_vector):
        probs_dict[p].append(prb_vector[i])
    probs = np.array([sum(probs_dict[p]) for p in range(-1000, 1001)])
    print("simulation step is done")
    return mutation_count, mutation_real_pos.tolist(), mutation_rand_pos.tolist(), probs


def run_tf(signature, data_folder, output_file, pattern, pool):

    # If output files exists skip
    if os.path.exists(output_file + 'rand.bin') and os.path.exists(output_file + 'real.bin'):
        logging.info("%s [skip]", output_file)
        return

    to_run = []
    for file in os.listdir(data_folder):
        if fnmatch.fnmatch(file, pattern):
            to_run.append((os.path.join(data_folder, file), signature))

    total_mutations = 0
    total_real_pos = []
    total_rand_pos = []
    total_prob = None
    for count, real_pos, rand_pos, prob in pool.starmap(simulation, to_run):
        total_mutations += count
        total_real_pos += real_pos
        total_rand_pos += rand_pos
        if total_prob is None:
            total_prob = prob
        elif prob is not None:
            total_prob = (total_prob + prob) / 2

    # Store output
    np.array(total_rand_pos).astype(np.int32).tofile(output_file + "rand.bin")
    #np.array(total_real_pos).astype(np.int32).tofile(output_file + "real.bin")
    #np.array(total_prob).astype(np.float32).tofile(output_file + "prob.bin")


def cmdline():

    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--signature', dest='signature_file', help='Signature file with the probabilities')
    parser.add_argument('-c', '--column', dest='signature_column', help='Signature probability column')
    parser.add_argument('-f', '--folder', dest='tf_folder', help='Folder with transcription factor files')
    parser.add_argument('-p', '--pattern', dest='pattern', help='Transcription factor files pattern')
    parser.add_argument('-o', '--output', dest='output_file', help='Output file')
    parser.add_argument('--cores', dest='cores', type=int, default=os.cpu_count(), help='Maximum number of CPU cores to use')
    parser.add_argument('--debug', dest='debug', action='store_true')
    parser.add_argument('--multiple', dest='multiple_tf', action='store_true', help='Execute all subfolders as a different TF')
    args = parser.parse_args()

    # Configure the logging
    level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=level)
    logging.debug(args)

    # Load the signature
    signature = pd.read_csv(args.signature_file, sep='\t')
    signature = signature[[REF, ALT, args.signature_column]]
    signature.columns = [REF, ALT, PROBABILITY]
    signature.set_index([REF, ALT], inplace=True)

    pool = Pool(processes=args.cores)
    if not args.multiple_tf:
        logging.debug("Start")
        run_tf(signature, args.tf_folder, args.output_file, args.pattern, pool)
        logging.debug("Done")
    else:
        subfolders = [o for o in os.listdir(args.tf_folder) if os.path.isdir(os.path.join(args.tf_folder, o))]
        for s in subfolders:
            logging.debug("Start %s", s)
            run_tf(signature, os.path.join(args.tf_folder, s), args.output_file + s + "_", args.pattern, pool)
            logging.debug("Done %s", s)

if __name__ == "__main__":
    cmdline()