#!/usr/bin/env python
"""Run a SVM on collected alignment data
"""
import sys
from svm_analysis import run_svm_on_motif
from argparse import ArgumentParser
from multiprocessing import Process, current_process, Manager


def parse_args():
    parser = ArgumentParser(description=__doc__)

    # query files
    parser.add_argument('--C_files', '-c', action='store',
                        dest='c_files', required=True, type=str, default=None,
                        help="directory with C files")
    parser.add_argument('--mC_files', '-mc', action='store',
                        dest='mc_files', required=True, type=str, default=None,
                        help="directory with mC files")
    parser.add_argument('--hmC_files', '-hmc', action='store',
                        dest='hmc_files', required=True, type=str, default=None,
                        help="directory with hmC files")
    parser.add_argument('-nb_files', '-nb', action='store', dest='nb_files', required=False,
                        default=50, type=int, help="maximum number of reads to align")
    parser.add_argument('--jobs', '-j', action='store', dest='jobs', required=False,
                        default=4, type=int, help="number of jobs to run concurrently")
    parser.add_argument('--iter', '-i', action='store', dest='iter', required=False,
                        default=4, type=int, help="number of iterations to do")
    parser.add_argument('--train_test', '-s', action='store', dest='split', required=False,
                        default=0.5, type=float, help="train/test split")
    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put results")
    parser.add_argument('--weighted', '-w', action='store', dest='weighted', default=False,
                        help="use weighted samples?")
    parser.add_argument('--kernel', '-k', action='store', dest='kernel', default='linear',
                        required=False, type=str)
    parser.add_argument('--C', '-C', action='store', dest='C', default=1.0,
                        required=False, type=float)
    parser.add_argument('--forward', '-f', action='store', dest='forward', default=True,
                        help='forward mapped reads?')
    args = parser.parse_args()
    return args


def run_svm(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            run_svm_on_motif(**f)
    except Exception:
        done_queue.put("%s failed" % current_process().name)


def main(args):
    args = parse_args()

    start_message = """
    Starting SVM analysis.
    Looking at {nbFiles} files.
    Forward mapped strand: {forward}.
    Using weights: {weights}.
    Iterations: {iter}.
    Train/test split: {train_test}.
    Kernel: {kernel}
    Output to: {out}""".format(nbFiles=args.nb_files, forward=args.forward, weights=args.weighted,
                               iter=args.iter, train_test=args.split, kernel=args.kernel, out=args.out)

    print >> sys.stderr, start_message

    motifs = [747, 354, 148, 796, 289, 363, 755, 626, 813, 653, 525, 80, 874]

    workers = args.jobs
    work_queue = Manager().Queue()
    done_queue = Manager().Queue()
    jobs = []

    for motif in motifs:
        svm_args = {
            "c_files": args.c_files,
            "mc_files": args.mc_files,
            "hmc_files": args.hmc_files,
            "weighted": args.weighted,
            "forward": args.forward,
            "ref_start": motif,
            "train_test_split": args.split,
            "iterations": args.iter,
            "out_path": args.out,
            "kernel": args.kernel,
            "max_samples": args.nb_files,
            "C": args.C
        }
        work_queue.put(svm_args)

    for w in xrange(workers):
        p = Process(target=run_svm, args=(work_queue, done_queue))
        p.start()
        jobs.append(p)
        work_queue.put('STOP')

    for p in jobs:
        p.join()

    done_queue.put('STOP')


if __name__ == "__main__":
    sys.exit(main(sys.argv))







