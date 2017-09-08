def findMaximumDValue(predictionListStats, experimentalDataStats, logOfFactorialOfPredictionListStats, returnlog = FALSE):
    '''Helper function for cubic1b: Find maximum D value - computes the maximum possible D-value for given values q+, q-, q0 and n+, n-, n0.
    @param predictionListStats a vector containing the predicted values q+, q- and q0:
    numbers of positive, negative and non-significant/contradictory predictions
    @param experimentalDataStats A vector containing the observed values n+, n- and n0:
    numbers of positive, negative and non-significant/contradictory observations
    @param logOfFactorialOfPredictionListStats a vector containing the log of the factorial value for
    each entry in predictionListStats
    @param returnlog should the result be returned as a log; default FALSE
    @return the maximum possible D value'''

    twoByTwoContingencyTable = findApproximateValuesThatWillMaximiseDValue(predictionListStats, experimentalDataStats)
    
    maximumDValue = getMaximumDValueFromTwoByTwoContingencyTable(twoByTwoContingencyTable, predictionListStats, experimentalDataStats, logOfFactorialOfPredictionListStats, 
        returnlog)
    
    return maximumDValue