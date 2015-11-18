#!/usr/bin/env python
"""Run SVM on eventAlign data
"""
from __future__ import print_function
import numpy as np
from random import sample
from sklearn import svm, mixture, linear_model

c_ea = "../cPecan/tests/marginTest_10_30_15/eventalign/testC_NP.fasta.sorted.bam.tsv"
mc_files = "../cPecan/tests/testAlignments/mc100/tempFiles_alignment/"
hmc_files = "../cPecan/tests/testAlignments/hmc100/tempFiles_alignment/"

labels = []


def get_working_reads_(eventAlignTable, nb_reads_to_use):
    read_indexes = set([x[3] for x in eventAlignTable])
    working_reads = sample(read_indexes, nb_reads_to_use)
    working_reads = [int(x) for x in working_reads]
    return working_reads


def get_eventalign_motif_features(tsv, nb_reads_to_get, label, motif_start, split):
    # figure out how many reads are in the tsv
    eventalign_table = np.loadtxt(tsv, dtype=str, skiprows=1)

    # sample which reads to use, returns a list
    working_read_idxs = get_working_reads_(eventalign_table, nb_reads_to_get)

    # extract the lines for just the working reads
    working_data = [line for line in eventalign_table if int(line[3]) in working_read_idxs]

    # collect the motifs from those reads
    feature_vectors = np.zeros([len(working_read_idxs), 12])

    kmer_length = 6
    motif_range = range(motif_start, motif_start + kmer_length)
    labels = []

    for i, read in enumerate(working_read_idxs):
        feature_vector = np.zeros(12)
        read_data = [line for line in working_data if int(line[3]) == read]
        for entry in read_data:
            if entry[4] == "t" and int(entry[1]) in motif_range:
                e_index = motif_range.index(int(entry[1]))
                delta_mean = float(entry[6]) - float(entry[10])
                feature_vector[e_index] = delta_mean

            if entry[4] == "c" and int(entry[1]) in motif_range:
                e_index = motif_range.index(int(entry[1])) + 6
                delta_mean = float(entry[6]) - float(entry[10])
                feature_vector[e_index] = delta_mean
        feature_vectors[i:i + 1] = feature_vector
        labels.append(label)

    # split into training and test vectors
    split_point = int(split * len(feature_vectors))
    training_vectors = feature_vectors[:split_point]
    training_labels = labels[:split_point]
    test_vectors = feature_vectors[split_point:]
    test_labels = labels[split_point:]

    return training_vectors, training_labels, test_vectors, test_labels


def eventAlign_run_svm(c_tsv, mc_tsv, hmc_tsv, forward, motif_start, split, iterations, out_path, kernel,
                       max_reads, C):

    if forward:
        direction_label = ".forward."
    else:
        direction_label = ".backward."

    out_file = open(out_path + str(motif_start) + direction_label + str(iterations) + ".tsv", 'wa')

    # bin to hold accuracies for each iteration
    scores = []

    for i in xrange(iterations):
        labels = []
        c_train, c_train_labels, c_test, c_test_labels = get_eventalign_motif_features(tsv=c_tsv,
                                                                                       nb_reads_to_get=max_reads,
                                                                                       label=0,
                                                                                       motif_start=motif_start,
                                                                                       split=split)
        mc_train, mc_train_labels, mc_test, mc_test_labels = get_eventalign_motif_features(tsv=mc_tsv,
                                                                                           nb_reads_to_get=max_reads,
                                                                                           label=1,
                                                                                           motif_start=motif_start,
                                                                                           split=split)
        hmc_train, hmc_train_labels, hmc_test, hmc_test_labels = get_eventalign_motif_features(
                                                                                        tsv=hmc_tsv,
                                                                                        nb_reads_to_get=max_reads,
                                                                                        label=2,
                                                                                        motif_start=motif_start,
                                                                                        split=split)
        training_data = np.vstack((c_train, mc_train, hmc_train))
        training_labels = np.concatenate((c_train_labels, mc_train_labels, hmc_train_labels))

        clf = svm.SVC(kernel=kernel, C=C)

        clf.fit(training_data, training_labels)

        test_data = np.vstack((c_test, mc_test, hmc_test))
        test_labels = np.concatenate((c_test_labels, mc_test_labels, hmc_test_labels))

        score = clf.score(test_data, test_labels)

        scores.append(score)

    print(">{motif}\t{mean}".format(motif=motif_start, mean=np.mean(scores)), end='\n', file=out_file)






