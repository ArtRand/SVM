#!/usr/bin/env python
"""Run a SVM on collected alignment data
"""
from __future__ import print_function, division
import numpy as np
import os
from sklearn import svm
from random import shuffle
import matplotlib.pyplot as plt


def get_motif_ranges(ref_start, forward):
    kmer_length = 6
    if forward:
        template_motif_range = range(ref_start, ref_start + kmer_length)
        complement_start = 891 - (ref_start + (kmer_length - 1))
        complement_motif_range = range(complement_start, complement_start + 6)
        return template_motif_range, complement_motif_range
    if not forward:
        template_start = 891 - (ref_start + (kmer_length - 1))
        template_motif_range = range(template_start, template_start + 6)
        complement_motif_range = range(ref_start, ref_start + kmer_length)
        return template_motif_range, complement_motif_range


def cull_motif_features(start, tsv, forward):
    # load the tsv
    data = np.loadtxt(tsv, dtype=str)
    template_motif_range, complement_motif_range = get_motif_ranges(start, forward)

    # buld a feature vector that has the first 6 elements as the template features and the second
    # six elements as the complement features, the features are selected as the ones with the maximum
    # posterior probability
    feature_vector = np.zeros(12)
    feature_posteriors = np.zeros(12)

    for line in data:
        if line[4] == "t" and int(line[0]) in template_motif_range:
            # determine which event in the motif this is
            e_index = template_motif_range.index(int(line[0]))
            delta_mean = float(line[5]) - float(line[9])
            delta_noise = float(line[6]) - float(line[10])
            posterior = line[8]
            # if the posterior for this event is higher than the one we have previously seen,
            if posterior > feature_posteriors[e_index]:
                feature_vector[e_index] = delta_mean
                feature_posteriors[e_index] = posterior
        if line[4] == "c" and int(line[0]) in complement_motif_range:
            e_index = complement_motif_range.index(int(line[0])) + 6
            delta_mean = float(line[5]) - float(line[9])
            delta_noise = float(line[6]) - float(line[10])
            posterior = line[8]
            if posterior > feature_posteriors[e_index]:
                feature_vector[e_index] = delta_mean
                feature_posteriors[e_index] = posterior

    return np.mean(feature_posteriors), feature_vector


def collect_data_vectors(path, forward, labels, label, portion, motif_start, max_samples):
    """collects the training data
    """
    # collect the files
    if forward:
        tsvs = [x for x in os.listdir(path) if x.endswith(".forward.tsv")]
    else:
        tsvs = [x for x in os.listdir(path) if x.endswith(".backward.tsv")]

    # shuffle
    shuffle(tsvs)

    if max_samples < len(tsvs):
        tsvs = tsvs[:max_samples]

    # get the number of files we're going to use
    split_index = int(portion * len(tsvs))

    # container for training and test data
    train_data = np.zeros([split_index, 12])
    test_data = np.zeros([len(tsvs) - split_index, 12])

    # container for weights
    weights = np.zeros(split_index)

    for i, f in enumerate(tsvs[:split_index]):
        weight, vector = cull_motif_features(motif_start, path + f, forward)
        train_data[i:i + 1] = vector
        weights[i] = weight
        labels.append(label)

    for i, f in enumerate(tsvs[split_index:]):
        weight, vector = cull_motif_features(motif_start, path + f, forward)
        test_data[i:i+1] = vector

    return train_data, weights, labels, test_data


def run_svm_on_motif(c_files, mc_files, hmc_files, weighted, forward, ref_start, train_test_split, iterations,
                     out_path, kernel, max_samples, C):

    if forward:
        direction_label = ".forward."
    else:
        direction_label = ".backward."

    out_file = open(out_path + str(ref_start) + direction_label + str(iterations) + ".tsv", 'wa')

    # bin to hold accuracies for each iteration
    calls = []

    for i in xrange(iterations):
        labels = []
        c_train, c_weights, labels, c_test = collect_data_vectors(path=c_files, forward=forward, labels=labels,
                                                                  label=0, portion=train_test_split,
                                                                  motif_start=ref_start, max_samples=max_samples)

        mc_train, mc_weights, labels, mc_test = collect_data_vectors(path=mc_files, forward=forward, labels=labels,
                                                                     label=1, portion=train_test_split,
                                                                     motif_start=ref_start, max_samples=max_samples)

        hmc_train, hmc_weights, labels, hmc_test = collect_data_vectors(path=hmc_files, forward=forward, labels=labels,
                                                                        label=2, portion=train_test_split,
                                                                        motif_start=ref_start, max_samples=max_samples)
        training_data = np.vstack((c_train, mc_train, hmc_train))
        weights = np.concatenate((c_weights, mc_weights, hmc_weights))
        clf = svm.SVC(kernel=kernel, C=C)
        if weighted:
            clf.fit(training_data, labels, sample_weight=weights)
        else:
            clf.fit(training_data, labels)


        c_predictions = clf.predict(c_test)
        nb_c_correct = list(c_predictions).count(0)
        c_accuracy = list(c_predictions).count(0) / float(len(c_predictions))
        print(c_accuracy, end='\t', file=out_file)

        mc_predictions = clf.predict(mc_test)
        nb_mc_correct = list(mc_predictions).count(1)
        mc_accuracy = list(mc_predictions).count(1) / float(len(mc_predictions))
        print(mc_accuracy, end='\t', file=out_file)

        hmc_predictions = clf.predict(hmc_test)
        nb_hmc_correct = list(hmc_predictions).count(2)
        hmc_accuracy = list(hmc_predictions).count(2) / float(len(hmc_predictions))
        print(hmc_accuracy, end='\n', file=out_file)
        accuracy = (nb_c_correct + nb_mc_correct + nb_hmc_correct) / (len(c_predictions) + len(mc_predictions) +
                                                                      len(hmc_predictions))
        calls.append(accuracy)

    print(">{motif}\t{mean}".format(motif=ref_start, mean=np.mean(calls)), end='\n', file=out_file)

    return

