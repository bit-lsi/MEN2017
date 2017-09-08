# -*- coding: utf-8 -*-

import math


def get_approximate_maximum_d_value_from_three_by_two_contingency_table(three_by_two_contingency_table,
                                                                        prediction_list_stats,
                                                                        log_of_factorial_of_prediction_list_stats,
                                                                        returnlog=False):
    n_pz = prediction_list_stats[1] - (three_by_two_contingency_table[1] + three_by_two_contingency_table[2])
    n_mz = prediction_list_stats[2] - (three_by_two_contingency_table[3] + three_by_two_contingency_table[4])
    n_zz = prediction_list_stats[3] - (three_by_two_contingency_table[5] + three_by_two_contingency_table[6])
    three_by_three_contingency_table = [three_by_two_contingency_table, n_pz, n_mz, n_zz]
    weight = calculate_weight_given_values_in_three_by_three_contingency_table(three_by_three_contingency_table,
                                                                               log_of_factorial_of_prediction_list_stats,
                                                                               returnlog)
    return weight


def get_max_D_value_for_a_three_by_two_family(r_p, r_m, r_z, n_p, n_m, predictionListStats,
                                              logOfFactorialOfPredictionListStats, returnlog=False):
    total = n_p + n_m

    # Compute the values of n++, n+-, n-+, n--, n0+, and n0- that maximise the D-value; these correspond to the following formula:
    # n_ab is approximately equal to q_a*n_b/T, where T = n+ + n-, and a,b are either + or -.
    # See Assessing statistical significance in causal graphs, page 7 - the formula is not stated explicitly but follows from the logic of algorithm 1a).

    if total > 0:
        n_pp = math.ceil(r_p * n_p / total)
        n_pm = r_p - n_pp
        n_mp = math.ceil(r_m * n_p / total)

    # Check this rounding produces a valid combination i.e. n++ + n-+ <= n+
    if (n_mp + n_pp) > n_p:
        n_mp = math.floor(r_m * n_p / total)
        n_mm = r_m - n_mp
        n_zp = n_p - (n_pp + n_mp)
        n_zm = n_m - (n_pm + n_mm)
    else:
        n_pp = 0
        n_pm = 0
        n_mp = 0
        n_mm = 0
        n_zp = 0
        n_zm = 0

    threeByTwoContingencyTable = [n_pp, n_pm, n_mp, n_mm, n_zp, n_zm]

    maximumDFamValue = get_approximate_maximum_d_value_from_three_by_two_contingency_table(threeByTwoContingencyTable,
                                                                                           predictionListStats,
                                                                                           logOfFactorialOfPredictionListStats,
                                                                                           returnlog)

    return maximumDFamValue


def calculate_weight_given_values_in_three_by_three_contingency_table(threeByThreeContingencyTable,
                                                                      logOfFactorialOfPredictionListStats,
                                                                      returnlog=False):
    """Given the values in the three by three contingency table and the values of the number of positive/negative/non-significant predictions (q+, q-, q0) this function calculates the D-value (or weight).
    
    :param threeByThreeContingencyTable a 3x3 contingency table
    :param logOfFactorialOfPredictionListStats log of Factorial of prediction statistics
    :param returnlog should the result be returned as a log value. Default is FALSE.
    
    :return a D-value (or weight)
    """

    # The D value is defined as: (q+! / (n++! * n+-! * n+0!)) * (q-! / (n-+! * n--! * n-0!)) * (q0! / (n0+! * n0-! * n00!))
    if returnlog:
        return sum(logOfFactorialOfPredictionListStats) - sum(math.log(math.factorial(threeByThreeContingencyTable)))
    else:
        return (
            math.exp(
                sum(logOfFactorialOfPredictionListStats) - sum(math.log(math.factorial(threeByThreeContingencyTable)))))


def get_weights_above_hypothesis_score_for_a_three_by_two_table(
        weights, r_p, r_m, r_z, n_p, n_m, prediction_list_stats, experimentalDataStats,
        log_of_factorial_of_prediction_list_stats,
        hypothesisScore, logepsDMax, logDMax):
    """ Finds the D-Values (weights) from any 3x2 contingency tables that have a score above and including the
    hypothesis score. It also calculates the total weight, and returns a 2x1 vector of the two values.
    The ratio of these values is the p-value

    :param weights:
    :param r_p:
    :param r_m:
    :param r_z:
    :param n_p:
    :param n_m:
    :param prediction_list_stats:
    :param experimentalDataStats:
    :param log_of_factorial_of_prediction_list_stats:
    :param hypothesisScore:
    :param logepsDMax:
    :param logDMax:
    :return:
    """

    log_approx_d_family_max = get_max_D_value_for_a_three_by_two_family(r_p, r_m, r_z, n_p, n_m, prediction_list_stats,
                                                                        log_of_factorial_of_prediction_list_stats, True)

    # To reduce the overall runtime of the algorithm, the approximate maxDValue for the family is calculated, rather than the actual.
    if log_approx_d_family_max > logepsDMax:
        # Calculate the values of n++ and n-+ that were used to calculate approxDFamilyMax

        # t = n_p + n_m ( the sum of all the values in the 3x2 submatrix)

        total = n_p + n_m

        if total > 0:
            n_pp = math.ceil(r_p * n_p / total)
        else:
            n_pp = 0

        # Iterate over possible values of n++, starting from the one calculated above.  Since the values decrease monotonically either side of this value,
        # we start at this value and iterate in both directions until the D value is less than a given threshold. This will give us the non-neglible parts
        # of the distribution.

        # An upper bound for possible values of n++ is the minimum of r+ and n+
        for m_pp in range(n_pp, min(r_p, n_p)):

            # We still need to define one more value in the 3x2 matrix to determine all the rest n+- is bounded by the minumum of r- and n+
            for m_mp in range(0, min(r_m, n_p)):

                # n+-
                m_pm = r_p - m_pp
                # n--
                m_mm = r_m - m_mp
                # n0+
                m_zp = n_p - (m_pp + m_mp)
                # n0-
                m_zm = n_m - (m_pm + m_mm)
                contingencyTableValues = [m_pp, m_pm, m_mp, m_mm, m_zp, m_zm]

                m_pz = prediction_list_stats[1] - r_p
                m_mz = prediction_list_stats[2] - r_m
                m_zz = prediction_list_stats[3] - r_z

                # None of these values can be less than zero
                if min(contingencyTableValues) < 0:
                    continue

                threeByThreeContingencyTable = [m_pp, m_pm, m_pz, m_mp, m_mm, m_mz, m_zp, m_zm, m_zz]
                logDValue = sum(log_of_factorial_of_prediction_list_stats) - sum(
                    math.log(math.factorial(threeByThreeContingencyTable)))
                # Only need to compute the score if the D-value is non-neglible
                if logDValue <= logepsDMax:
                    continue

                weight_contrib = math.exp(logDValue - logDMax)
                weights[2] = weights[2] + weight_contrib

                # Rather than just calculate the DValue by taking the exponential of logDvalue, we first subtract a constant.  This means when we take the
                # exponential we get the DValue divided by a different constant, but this cancels out when we take ratios later when computing the p-value
                if m_pp + m_mm - (m_pm + m_mp) >= hypothesisScore:
                    weights[1] = weights[1] + weight_contrib

            # The loop below does exactly the same thing as the one above; it covers the cases that are not done in the loop above. A lower bound for possible
            # values of n++ is 0 To avoid n_pp become negative or double counting when n_pp = 0, we add the following if statement
            if (n_pp > 0):
                for m_pp in range(n_pp - 1, 0):

                    # Iterate over all values of m_mp
                    for m_mp in range(0, min(r_m, n_p)):

                        # n+-
                        m_pm = r_p - m_pp
                        # n--
                        m_mm = r_m - m_mp
                        # n0+
                        m_zp = n_p - (m_pp + m_mp)
                        # n0-
                        m_zm = n_m - (m_pm + m_mm)
                        contingencyTableValues = [m_pp, m_pm, m_mp, m_mm, m_zp, m_zm]

                        m_pz = prediction_list_stats[1] - r_p
                        m_mz = prediction_list_stats[2] - r_m
                        m_zz = prediction_list_stats[3] - r_z

                        # None of these values can be less than zero
                        if min(contingencyTableValues) < 0:
                            continue

                        threeByThreeContingencyTable = [m_pp, m_pm, m_pz, m_mp, m_mm, m_mz, m_zp, m_zm, m_zz]
                        logDValue = math.log(math.factorial(log_of_factorial_of_prediction_list_stats)) - math.log(
                            math.factorial(threeByThreeContingencyTable))

                        # Only need to compute the score if the D-value is non-neglible
                        if logDValue <= logepsDMax:
                            continue

                        # Rather than just calculate the DValue by taking the exponential of logDvalue, we first subtract a constant.  This means when we take the
                        # exponential we get the DValue divided by a different constant, but this cancels out when we take ratios later when computing the p-value
                        weight_contrib = math.exp(logDValue - logDMax)
                        weights[2] = weights[2] + weight_contrib

                        if (m_pp + m_mm - (m_pm + m_mp) >= hypothesisScore):
                            weights[1] = weights[1] + weight_contrib

    return weights
