# -*- coding: utf-8 -*-

import math

from calculate_d_value import find_maximum_d_value, get_weights_above_hypothesis_score_for_3_2_table


def cubic1b(hypothesis_score, prediction_stats, experimental_data_stats, epsilon):
    '''Calculate the p-value of a score given the hypothesis score and the
    distribution table (calculated using the cubic algorithm  1b in
    Assessing statistical significance in causal graphs - Chindelevitch et al)
    '''
    # Calculate the log of the factorial for all values in predictionListStats,
    # this is to reduce the number of computations done later on.
    log_factorial_prediction_stats = [math.log(math.factorial(n)) for n in prediction_stats]
    logDMax = find_maximum_d_value(prediction_stats, experimental_data_stats, log_factorial_prediction_stats, True)

    q_p, q_m = prediction_stats[0:2]
    n_p, n_m = experimental_data_stats[0:2]

    # Array for the D values - first element is the total D values which have a score equal to or better than the cut-off, second element is all D
    # values (i.e. the ratio is the P value).
    logepsDMax = math.log(epsilon) + logDMax  # or [n+math.log(epsilon) for n in logDMax] #
    weights = [0, 0]

    # In algorithm a superfamily is defined by the left 3x2 submatrix of the contingency table. Now we only need to interate over two quantities: r+ and
    # r-.
    for r_p in range(0, q_p):
        # for 0 to the limit of positive predictions
        for r_m in range(0, q_m):
            # for 0 to the limit of number of negative results Calculate r0
            r_z = n_p + n_m - (r_p + r_m)
            # Need to check that n00 is postive, to ensure it is a valid superfamily
            if r_z <= prediction_stats[2] and r_z >= 0:
                # The five values in the array below define the superfamily
                weights = get_weights_above_hypothesis_score_for_3_2_table(weights, r_p, r_m, r_z, n_p, n_m,
                                                                           prediction_stats, experimental_data_stats,
                                                                           log_factorial_prediction_stats,
                                                                           hypothesis_score, logepsDMax, logDMax)

    p_value = weights[0] / weights[1]

    return p_value
