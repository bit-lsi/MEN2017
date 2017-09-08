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

