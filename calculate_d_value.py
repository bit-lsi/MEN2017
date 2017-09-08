# -*- coding: utf-8 -*-

# ## Implementation of CalculateSignificanceUsingCubicAlgorithm1b.r
# 
# ### Adapted from CausalR

import numpy as np
import math

from miscellaneous import get_all_rounding_combinations
from weights import calculate_weight_given_values_in_three_by_three_contingency_table

def find_approximate_values_maximise_d_value(prediction_stats, experimental_data_stats):
    '''Helper function for findMaximumDValue: Finds an approximate table values to maximise D.
    Given the values of q+, q-, q0, n+, n- and n0 this function will produce the approximate values
    of n++, n+-, n-+ and n-- that will maximise the D value. See Assessing statistical significance of 
    casual graphs, page 6. The values are approximate since they need to be rounded, although the direction
    of rounding is not clear at this stage.

    @param predictionListStats a vector containing the values q+, q- and q0: 
        numbers of positive, negative and non-significant/contradictory predictions
    @param experimentalDataStats a vector containing the values n+, n- and n0:
        numbers of positive, negative and non-significant/contradictory predictions
    @return a 2x2 contingency table which approximately maximises D'''

    q_p, q_m = prediction_stats[0:2]

    n_p, n_m = experimental_data_stats[0:2]

    Tval = sum(prediction_stats)
    # The values of n++, n+-, n-+ and n-- that give the maximum D-value are given by the formula within the paper - Assessing statistical significance
    # in causal graphs, page 6.  The formula is n_ab is approximately equal to q_a*n_b/T, where T = q+ + q- + q0 = n+ + n- + n0, and a,b are either + or
    # -.

    n_pp = q_p * n_p / Tval
    n_pm = q_p * n_m / Tval
    n_mp = q_m * n_p / Tval
    n_mm = q_m * n_m / Tval

    contingency_table = [n_pp, n_pm, n_mp, n_mm]

    return contingency_table


def get_maximum_d_value_from_2_x_2_contingency_table(contingency_table, prediction_stats, experimental_data_stats,
                                                     log_of_factorial_prediction_stats,
                                                     return_log):
    '''Helper function for findMaximumDValue: Get maximum D value from two-by-two contingency table 
        - computes the maximum D value (or weight) given approximate values of n++, n+-, n-+ and n--.
        These values are approximate and in general are non-integer values; they are found by using an 
        approximation detailed in the paper Assessing statistical significance in causal graphs on page 6 i.e. n_ab is approximately
        equal to q_a*n_b/t where a and b are either +, - or 0. The value is an approximation since the direction in which the 
        number should be rounded is not clear and hence this function runs through all possible combinations of rounding before
        concluding the maximum D-value.

    @param contingency_table approximate values of n++, n+-, n-+ and n--, these values arecalculated to optimise the D-value
    @param prediction_stats a vector containing the values q+, q- and q0 the number of positive/negative/non-significant (or contradictory) predictions)
    @param experimental_data_stats a vector containing the values n+, n- and n0 (the number of positive/negative/non-significant (or contradictory) transcripts in the results)
    @param log_of_factorial_prediction_stats a vector containing the log of the factorial value for each entry in predictionListStats
    @param return_log whether or not he value should be returned as a log (TRUE) or not (FALSE)
    @return the maximal D-value'''

    # We need to find all the combinations (maximum of 16) when applying the rounding to the values in twoByTwoContingencyTable.  The exception is if
    # all four values are integers.  combinations_rounding will be a nx4 matrix containing the possible values of n_pp, n_pm, n_mp, n_mm
    if np.round(contingency_table) == contingency_table:
        # All four values are integers --> twoByTwoContingencyTable = [n_pp, n_pm, n_mp, n_mm] (c(...) in R)
        combinations_rounding = np.matrix(contingency_table)  # , ncol = 4)
    else:
        combinations_rounding = get_all_rounding_combinations(contingency_table)

    maximumDValue = 0

    for i in range(combinations_rounding.shape[0]):  # nrwo in R (loop over number of rows of the matrix)
        numbers_correctand_incorrect_predictions = combinations_rounding[
            i]  # get all elements of the i-th row ([i,] in R)
        three_x_three_contingency_table = populateTheThreeByThreeContingencyTable(
            numbers_correctand_incorrect_predictions[0, 0], numbers_correctand_incorrect_predictions[0, 1],
            numbers_correctand_incorrect_predictions[0, 2], numbers_correctand_incorrect_predictions[0, 3],
            prediction_stats, experimental_data_stats)
        # Some of the combinationations produce an invalid contingency table - this we can check by seeing if any of the values in the table are negative.
        if three_x_three_contingency_table.min() >= 0:  # if 3x3ContTable is matrix or nparray! otherwise change to min(x)
            weight = calculate_weight_given_values_in_three_by_three_contingency_table(three_x_three_contingency_table,
                                                                              log_of_factorial_prediction_stats,
                                                                              return_log)
            # what type is weight? If not integer/float, max needs to be changed (e.g. for list max(max(weight), maximumDValue))
            maximumDValue = max(maximumDValue, weight)

    return maximumDValue


def find_maximum_d_value(prediction_stats, experimental_data_stats, log_factorial_prediction_stats, returnlog=False):
    '''Helper function for cubic1b: Find maximum D value - computes the maximum possible D-value for given values q+, q-, q0 and n+, n-, n0.
    @param prediction_stats a vector containing the predicted values q+, q- and q0:
    numbers of positive, negative and non-significant/contradictory predictions
    @param experimental_data_stats A vector containing the observed values n+, n- and n0:
    numbers of positive, negative and non-significant/contradictory observations
    @param log_factorial_prediction_stats a vector containing the log of the factorial value for
    each entry in prediction_stats
    @param returnlog should the result be returned as a log; default FALSE
    @return the maximum possible D value'''

    contingency_table = find_approximate_values_maximise_d_value(prediction_stats, experimental_data_stats)

    maximum_d_value = get_maximum_d_value_from_2_x_2_contingency_table(contingency_table, prediction_stats,
                                                                     experimental_data_stats,
                                                                     log_factorial_prediction_stats,
                                                                     returnlog)

    return maximum_d_value


def get_weights_above_hypothesis_score_for_3_2_table(weights, r_p, r_m, r_z, n_p, n_m, predictionListStats,
                                                     experimentalDataStats,
                                                     logOfFactorialOfPredictionListStats, hypothesisScore, logepsDMax,
                                                     logDMax):
    '''Helper function for cubic1b: '''
    return [0, 0]


