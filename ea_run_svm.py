#!/usr/bin/env python
"""Run a SVM on collected alignment data
"""
import sys
from ea_svm_analysis import eventAlign_run_svm
from argparse import ArgumentParser
from multiprocessing import Process, current_process, Manager


def parse_args():
    parser = ArgumentParser(description=__doc__)

    # query files
    parser.add_argument('--C_tsv', '-c', action='store',
                        dest='c_tsv', required=True, type=str, default=None,
                        help="file with C alignments")
    parser.add_argument('--mC_tsv', '-mc', action='store',
                        dest='mc_tsv', required=True, type=str, default=None,
                        help="file with mC alignments")
    parser.add_argument('--hmC_tsv', '-hmc', action='store',
                        dest='hmc_tsv', required=True, type=str, default=None,
                        help="file with hmC alignments")
    parser.add_argument('--kernel', '-k', action='store', dest='kernel', default='linear',
                        required=False, type=str)
    parser.add_argument('--C', '-C', action='store', dest='C', default=1.0,
                        required=False, type=float)
    parser.add_argument('--backward', '-bw', action='store_false', dest='forward',
                        default=True, help='forward mapped reads?')
    parser.add_argument('-nb_files', '-nb', action='store', dest='nb_files', required=False,
                        default=50, type=int, help="maximum number of reads to use")
    parser.add_argument('--jobs', '-j', action='store', dest='jobs', required=False,
                        default=4, type=int, help="number of jobs to run concurrently")
    parser.add_argument('--iter', '-i', action='store', dest='iter', required=False,
                        default=4, type=int, help="number of iterations to do")
    parser.add_argument('--train_test', '-s', action='store', dest='split', required=False,
                        default=0.5, type=float, help="train/test split")
    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put results")
    args = parser.parse_args()
    return args


def run_svm(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            eventAlign_run_svm(**f)
    except Exception:
        done_queue.put("%s failed" % current_process().name)


def main(args):
    args = parse_args()

    start_message = """
    Running SVM on SignalAlign Data!
    Command line: {cmd}
    Starting SVM analysis.
    Looking at {nbFiles} files.
    Forward mapped strand: {forward}.
    Iterations: {iter}.
    Train/test split: {train_test}.
    Kernel: {kernel}
    Output to: {out}""".format(nbFiles=args.nb_files, forward=args.forward,
                               iter=args.iter, train_test=args.split, kernel=args.kernel, out=args.out,
                               cmd=" ".join(sys.argv[:]))

    print >> sys.stderr, start_message

    motifs = [747, 354, 148, 796, 289, 363, 755, 626, 813, 653, 525, 80, 874]

    workers = args.jobs
    work_queue = Manager().Queue()
    done_queue = Manager().Queue()
    jobs = []

    for motif in motifs:
        svm_args = {
            "c_tsv": args.c_tsv,
            "mc_tsv": args.mc_tsv,
            "hmc_tsv": args.hmc_tsv,
            "forward": args.forward,
            "motif_start": motif,
            "split": args.split,
            "iterations": args.iter,
            "out_path": args.out,
            "kernel": args.kernel,
            "max_reads": args.nb_files,
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

    print >> sys.stderr, "\n\tFinished SVM"


if __name__ == "__main__":
    sys.exit(main(sys.argv))







