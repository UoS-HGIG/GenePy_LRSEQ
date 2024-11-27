import numpy as np
import pandas as pd
from numba import njit, cuda, prange
import math
import pyarrow.csv as pa_csv
from typing import Tuple
import argparse


def read_meta_file(filepath: str) -> pd.DataFrame:
    parse_options = pa_csv.ParseOptions(delimiter="\t")
    read_options = pa_csv.ReadOptions(block_size=1e9)
    data = pa_csv.read_csv(
        filepath, parse_options=parse_options, read_options=read_options
    )
    data = data.to_pandas()
    return data


def get_scores(data: np.array) -> np.array:
    ##cadd score; reformat CADD range to 0-1
    ## the UKBB 200k cohort raw variant has CADD16 ranges from -18.793437 to 19.100986 including the prescore and permutated score (GRCh38_v1.6;VEP 100_GRCh38)
    scores = data[:, 16:25]
    scores = scores.astype("float")
    scores = (scores - (-19.811548)) / (25.028523 - (-19.811548))
    scores[np.isnan(scores)] = 0
    s0 = scores[:, 0]
    scores = np.insert(scores, 0, s0, axis=1)
    scores[np.isnan(scores)] = 0

    return scores


def get_allele_freq(data: np.array) -> np.array:
    ##allele frequency as it is in the UKBB cohort; this is currently based on the raw pVCF data.
    af = data[:, 6:15]
    af = af.astype("float")
    ##the ref allele frequency
    af0 = 1 - np.nansum(af, axis=1)
    af0 = af0[:, np.newaxis]
    af[np.isnan(af)] = 1
    af = np.hstack([af0, af])

    return af


def format_data(data: pd.DataFrame) -> Tuple[np.array, np.array, np.array, np.array]:
    header = data.columns.values
    data = np.array(data)
    scores = get_scores(data=data)
    af = get_allele_freq(data=data)
    samples_header = header[26:]
    samples = data[:, 26:]
    samples[samples == "0"] = "0/0"
    samples = samples.astype("str")

    return scores, af, samples, samples_header


def parse_command_line_args():
    parser = argparse.ArgumentParser(description="GenePy2 - Make score matrix")
    parser.add_argument("--gene", type=str, help="Gene name", required=True)
    parser.add_argument("--cadd", type=str, help="CADD score", required=True)
    parser.add_argument("--gpu", action="store_true", help="Use GPU")
    parser.add_argument(
        "--gpu-threads",
        type=int,
        help="Number of threads in each GPU block",
        nargs="?",
        const=256,
        default=256,
    )
    args = parser.parse_args()
    return args


@njit()
def get_score(S: np.array, af: np.array, db1: np.array, s_int: np.array):
    for i in range(db1.shape[0]):
        for j in prange(db1.shape[1]):
            if s_int[i, j, 0] == 0:
                db1[i][j] = S[i, s_int[i][j][1]] * (
                    -math.log10(af[i, s_int[i][j][0]] * af[i, s_int[i][j][1]])
                )

            elif s_int[i, j, 1] == 0:
                db1[i][j] = S[i, s_int[i][j][0]] * (
                    -math.log10(af[i, s_int[i][j][0]] * af[i, s_int[i][j][1]])
                )
            else:
                db1[i][j] = (
                    (S[i, s_int[i][j][0]] + S[i, s_int[i][j][1]])
                    * 0.5
                    * (-math.log10(af[i, s_int[i][j][0]] * af[i, s_int[i][j][1]]))
                )
    return db1


@cuda.jit
def get_score_kernel(S_d, db1_d, af_d, s_int_d):
    start = cuda.grid(1)
    stride = cuda.gridsize(1)

    for i in range(db1_d.shape[0]):
        for j in range(start, db1_d.shape[1], stride):

            if s_int_d[i, j, 0] == 0:
                db1_d[i][j] = S_d[i, s_int_d[i][j][1]] * (
                    -math.log10(af_d[i, s_int_d[i][j][0]] * af_d[i, s_int_d[i][j][1]])
                )

            elif s_int_d[i, j, 1] == 0:
                db1_d[i][j] = S_d[i, s_int_d[i][j][0]] * (
                    -math.log10(af_d[i, s_int_d[i][j][0]] * af_d[i, s_int_d[i][j][1]])
                )
            else:
                db1_d[i][j] = (
                    (S_d[i, s_int_d[i][j][0]] + S_d[i, s_int_d[i][j][1]])
                    * 0.5
                    * (-math.log10(af_d[i, s_int_d[i][j][0]] * af_d[i, s_int_d[i][j][1]]))
                )


def get_scores_cuda(
    scores: np.array,
    af: np.array,
    db1: np.array,
    samples_int: np.array,
    threads_per_block: int,
):

    num_cols = db1.shape[1]
    S_d = cuda.to_device(np.ascontiguousarray(scores))
    db1_d = cuda.to_device(np.ascontiguousarray(db1))
    af_d = cuda.to_device(np.ascontiguousarray(af))
    s_int_d = cuda.to_device(np.ascontiguousarray(samples_int))

    blockspergrid = (num_cols + threads_per_block - 1) // threads_per_block
    get_score_kernel[blockspergrid, threads_per_block](S_d, db1_d, af_d, s_int_d)

    S_d = None
    af_d = None
    s_int_d = None
    db1 = db1_d.copy_to_host()
    db1_d = None

    return db1


def index(data: np.array) -> np.array:
    int_data = data.ravel().view("uint8")[::4]
    int_data = int_data.reshape(data.shape[0], -1, 3) - 48
    int_data[(int_data[:, :, 2] > 253) | (int_data[:, :, 0] > 253)] = [0, 0, 0]
    int_data = int_data[:, :, [0, 2]].reshape(int_data.shape[0], int_data.shape[1], 2)
    return int_data.astype("uint8")


def nan_if(arr, value):
    return np.where(arr == value, np.nan, arr)


def score_db(
    samples: np.array,
    scores: np.array,
    af: np.array,
    gpu: bool,
    gpu_threads_per_block: int,
):
    samples_int = index(samples)
    db1 = np.zeros_like(samples, dtype=float)
    if gpu:
        db1 = get_scores_cuda(
            scores, af, db1, samples_int, threads_per_block=gpu_threads_per_block
        )
    else:
        db1 = get_score(scores, af, db1, samples_int)
    db1[np.all(samples_int == np.asarray([0.0, 0.0], dtype="uint8"), axis=-1)] = 0.0
    out1 = np.nansum(nan_if(db1, "0.0"), axis=0)
    gg = np.array([gene] * len(samples_header))
    U = np.vstack((samples_header, out1, gg)).T
    return U


if __name__ == "__main__":
    args = parse_command_line_args()
    meta_file = args.gene + ".meta"
    gene = args.gene
    cadd = args.cadd
    gpu = args.gpu
    gpu_threads_per_block = args.gpu_threads
    data = read_meta_file(filepath=meta_file)
    scores, af, data, samples_header = format_data(data=data)

    if (np.isnan(scores).sum()) < (
        scores.shape[0]
    ):  # compute metascores if at least 1 variant
        U = score_db(
            samples=data,
            scores=scores,
            af=af,
            gpu=gpu,
            gpu_threads_per_block=gpu_threads_per_block,
        )
        np.savetxt(
            "./" + args.gene + "_" + args.cadd + "_matrix", U, fmt="%s", delimiter="\t"
        )
