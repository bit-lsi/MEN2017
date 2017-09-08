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
