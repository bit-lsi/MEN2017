# -*- coding: utf-8 -*-

import numpy as np
import scipy.stats as stats


def get_set_of_significant_predictions(predictions):
    """  Gets the set of positive and negative predictions, the combination of the sets Sh+ and Sh-

    :param predictions:
    :return:
    """
    counter = 0
    num_predictions = predictions.shape[0]
    significant_predictions = np.zeros((num_predictions, 1))
    for i in range(0, num_predictions):
        if int(predictions[i, 1] != 0):
            significant_predictions[counter, 0] = predictions[i, 0]
            counter = counter + 1

    significant_predictions = significant_predictions[1:counter]

    return significant_predictions


def get_set_of_differentially_expressed_genes(results):
    """ Gets the set of differentially expressed genes in the results, G+ as defined by in Causal reasoning on biological networks

    :param results:
    :return:
    """
    counter = 0
    num_results = results.shape[0]
    differentially_expressed_genes = np.zeros((num_results, 1))
    for i in range(0, num_results):
        if int(results[i, 1] != 0):
            differentially_expressed_genes[counter, 0] = results[i, 0]
            counter = counter + 1

    differentially_expressed_genes = differentially_expressed_genes[0:counter]

    return differentially_expressed_genes


def calculate_enrichment_p_value(predictions, results):
    """ Calculate a enrichment p-value for a given hypothesis by comparing the corresponding predicted and observed
     gene changes

    :param predictions:
    :param numpy.array results:
    :return:
    """

    set_of_significant_predictions = get_set_of_significant_predictions(predictions)
    set_of_differentially_expressed_genes = get_set_of_differentially_expressed_genes(results)

    set_of_non_differentially_expressed_genes = np.setdiff1d(results[:, 0], set_of_differentially_expressed_genes)

    # n_pp + n_pm + n_mp + n_mm
    num_significant_predictions_that_are_responsive = len(
        np.intersect1d(set_of_significant_predictions, set_of_differentially_expressed_genes))
    # n+0 + n-0
    num_significant_predictions_that_are_unresponsive = len(
        np.intersect1d(set_of_significant_predictions, set_of_non_differentially_expressed_genes))
    # n0+ + n0-
    num_zeroPredictions_that_are_responsive = len(
        set_of_differentially_expressed_genes) - num_significant_predictions_that_are_responsive
    # n00
    num_zeroPredictions_that_are_unreponsive = len(
        set_of_non_differentially_expressed_genes) - num_significant_predictions_that_are_unresponsive

    contingency_table = np.array(
        [[num_significant_predictions_that_are_responsive, num_significant_predictions_that_are_unresponsive],
         [num_zeroPredictions_that_are_responsive,
          num_zeroPredictions_that_are_unreponsive]])

    enrichment_p_value = stats.fisher_exact(contingency_table, alternative="greater")

    return enrichment_p_value


if __name__ == '__main__':
    predictions = np.array([[1, 2, 3, 1, 1, -1], [1, 2, 3, 1, 1, -1]])
    results = np.array([[1, 2, 3, 4, 1, 1, -1, 1], [1, 2, 3, 4, 1, 1, -1, 1]])
    calculate_enrichment_p_value(predictions, results)
