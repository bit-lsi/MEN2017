# -*- coding: utf-8 -*-

import math
import numpy as np


def get_all_rounding_combinations(contingency_table):
    """ Returns all possible rounding combinations of a 2x2 table.

    :param contingency_table:
    :return:
    """

    floor_pos_pos = math.floor(contingency_table[0])
    floor_pos_neg = math.floor(contingency_table[1])
    floor_neg_pos = math.floor(contingency_table[2])
    floor_neg_neg = math.floor(contingency_table[3])

    ceeling_pos_pos = math.ceil(contingency_table[0])
    ceeling_pos_neg = math.ceil(contingency_table[1])
    ceeling_neg_pos = math.ceil(contingency_table[2])
    ceeling_neg_neg = math.ceil(contingency_table[3])

    round_combinations = np.array([
        [floor_pos_pos, floor_pos_neg, floor_neg_pos, floor_neg_neg],
        [floor_pos_pos, floor_pos_neg, floor_neg_pos, ceeling_neg_neg],
        [floor_pos_pos, floor_pos_neg, ceeling_neg_pos, floor_neg_neg],
        [floor_pos_pos, floor_pos_neg, ceeling_neg_pos, ceeling_neg_neg],
        [floor_pos_pos, ceeling_pos_neg, floor_neg_pos, floor_neg_neg],
        [floor_pos_pos, ceeling_pos_neg, floor_neg_pos, ceeling_neg_neg],
        [floor_pos_pos, ceeling_pos_neg, ceeling_neg_pos, floor_neg_neg],
        [floor_pos_pos, ceeling_pos_neg, ceeling_neg_pos, ceeling_neg_neg],
        [ceeling_pos_pos, floor_pos_neg, floor_neg_pos, floor_neg_neg],
        [ceeling_pos_pos, floor_pos_neg, floor_neg_pos, ceeling_neg_neg],
        [ceeling_pos_pos, floor_pos_neg, ceeling_neg_pos, floor_neg_neg],
        [ceeling_pos_pos, floor_pos_neg, ceeling_neg_pos, ceeling_neg_neg],
        [ceeling_pos_pos, ceeling_pos_neg, floor_neg_pos, floor_neg_neg],
        [ceeling_pos_pos, ceeling_pos_neg, floor_neg_pos, ceeling_neg_neg],
        [ceeling_pos_pos, ceeling_pos_neg, ceeling_neg_pos, floor_neg_neg],
        [ceeling_pos_pos, ceeling_pos_neg, ceeling_neg_pos, ceeling_neg_neg]]
    )

    keep = [True] * 16

    if floor_pos_pos == ceeling_pos_pos:
        keep[8:15] = False

    if floor_pos_neg == ceeling_pos_neg:
        keep[6:9] = False
        keep[12:15] = False

    if floor_neg_pos == ceeling_neg_pos:
        keep[2:3] = False
        keep[6:7] = False
        keep[10:11] = False
        keep[14:15] = False

    if (floor_neg_neg == ceeling_neg_neg):
        for i in range(0, 15, 2):
            keep[i] = False

    return round_combinations[keep,]


def find_approximate_values_that_will_maximise_D_value(predictionListStats, experimentalDataStats):
    """ Finds an approximate table values to maximise D.

    :param list predictionListStats: a list containing the values q+, q- and q0 which are numbers of positive, negative and non-significant/contradictory predictions
    :param list experimentalDataStats: a list containing the values n+, n- and n0 which are numbers of positive, negative and non-significant/contradictory predictions
    :rtype list
    :return twoByTwoContingencyTable: a list  which is a 2x2 contingency table which approximately maximises D
    """

    q_p = predictionListStats[0]
    q_m = predictionListStats[1]
    n_p = experimentalDataStats[0]
    n_m = experimentalDataStats[1]
    Tval = sum(predictionListStats)

    # The values of n++, n+-, n-+ and n-- that give the maximum D-value are given by the formula within the paper -
    # Assessing statistical significance in causal graphs, page 6.  The formula is n_ab is approximately equal to q_a*n_b/T, where T = q+ + q- + q0 = n+ + n- + n0, and a,b are either + or -.

    n_pp = q_p * n_p / Tval
    n_pm = q_p * n_m / Tval
    n_mp = q_m * n_p / Tval
    n_mm = q_m * n_m / Tval
    twoByTwoContingencyTable = [n_pp, n_pm, n_mp, n_mm]

    return twoByTwoContingencyTable


def compute_final_distribution(result_matrix):
    """ Computes a final reference distribution of the score used to compute the final p-value.

    :param numpy.array result_matrix: a numpy matrix
    :rtype numpy.array
    :return distributionMatrix
    """

    if result_matrix.size == 2:
        result_matrix = result_matrix.transpose()

    maxScore = max(result_matrix[:, 0])
    minScore = min(result_matrix[:, 0])

    # Pre-allocate the size of storage array
    distribution_matrix = np.zeros(shape=((maxScore - minScore + 1), 2))
    numRows = result_matrix.shape[0]
    counter = 0
    for score in range(minScore, maxScore + 1):
        probability = 0
        for i in range(numRows):
            if score == result_matrix[i, 0]:
                probability = probability + result_matrix[i, 1]
        if probability > 0:
            distribution_matrix[counter,] = [score, probability]
            counter += 1

    # Remove the rows that were not populated in the for loops above
    distribution_matrix = distribution_matrix[0:counter, :]

    return distribution_matrix
