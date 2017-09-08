# -*- coding: utf-8 -*-

import math
import numpy as np


def get_all_rounding_combinations(contingency_table):
    """ Returns all possible rounding combinations of a 2x2 table.

    :param contingency_table:
    :return:
    """

    contingency_table_floor = [
        math.floor(element)
        for element in contingency_table
    ]

    contingency_table_ceiling = [
        math.ceil(element)
        for element in contingency_table
    ]

    floor_pos_pos, floor_pos_neg, floor_neg_pos, floor_neg_neg = contingency_table_floor[0:4]

    ceeling_pos_pos, ceeling_pos_neg, ceeling_neg_pos, ceeling_neg_neg = contingency_table_ceiling[0:4]

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
