# ' @title returns approximate maximum D value or weight for a 3x2 superfamily
# ' @description
# ' Computes an approximate maximum D value (or weight) for a superfamily (3x2 table).
# ' The result is only approximate as only the first valid D value that is return. This has been
# ' done to speed up the overall algorithm.
# '
# ' @param threeByTwoContingencyTable  approximate values of n++, n+-, n-+, n--, n0+ and n0-, these values are
# ' calculated to optimise the D-value (see page 6 of Assessing
# ' statistical significance of causal graphs)
# ' @param predictionListStats a vector containing the values q+, q- and q0
# ' (the number of positive/negative/non-significant
# ' (or contradictory) predictions)
# ' @param logOfFactorialOfPredictionListStats a vector containing the log of the factorial value for
# ' each entry in predictionListStats
# ' @param returnlog return the result as a log, default is FALSE
# ' @return an approximate maximum D value or weight


def get_approximate_maximum_d_value_from_three_by_two_contingency_table(three_by_two_contingency_table, prediction_list_stats,
                                                                        log_of_factorial_of_prediction_list_stats,
                                                                        returnlog=False):


# *** Summary *** Description: A function that will compute an approximate maximum D value (or weight) for a superfamily (3x2 table).  The reason it
# is described as approiximate is that this function returns the first valid D value that is found. This has been done to speed up the overall
# algorithm.


    n_pz = prediction_list_stats[1] - (three_by_two_contingency_table[1] + three_by_two_contingency_table[2])
    n_mz = prediction_list_stats[2] - (three_by_two_contingency_table[3] + three_by_two_contingency_table[4])
    n_zz = prediction_list_stats[3] - (three_by_two_contingency_table[5] + three_by_two_contingency_table[6])
    three_by_three_contingency_table = [three_by_two_contingency_table, n_pz, n_mz, n_zz]
    weight = calculate_weight_given_values_in_three_by_three_contingency_table(three_by_three_contingency_table,
                                                                      log_of_factorial_of_prediction_list_stats, returnlog)
    return (weight)
