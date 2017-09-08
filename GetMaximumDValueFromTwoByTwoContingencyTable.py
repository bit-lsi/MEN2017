# -*- coding: utf-8 -*-

import numpy as np

def getMaximumDValueFromTwoByTwoContingencyTable(twoByTwoContingencyTable, predictionListStats, experimentalDataStats, logOfFactorialOfPredictionListStats,
        returnlog):
    '''Helper function for findMaximumDValue: Get maximum D value from two-by-two contingency table 
        - computes the maximum D value (or weight) given approximate values of n++, n+-, n-+ and n--.
        These values are approximate and in general are non-integer values; they are found by using an 
        approximation detailed in the paper Assessing statistical significance in causal graphs on page 6 i.e. n_ab is approximately
        equal to q_a*n_b/t where a and b are either +, - or 0. The value is an approximation since the direction in which the 
        number should be rounded is not clear and hence this function runs through all possible combinations of rounding before
        concluding the maximum D-value.

    @param twoByTwoContingencyTable approximate values of n++, n+-, n-+ and n--, these values arecalculated to optimise the D-value
    @param predictionListStats a vector containing the values q+, q- and q0 the number of positive/negative/non-significant (or contradictory) predictions)
    @param experimentalDataStats a vector containing the values n+, n- and n0 (the number of positive/negative/non-significant (or contradictory) transcripts in the results)
    @param logOfFactorialOfPredictionListStats a vector containing the log of the factorial value for each entry in predictionListStats
    @param returnlog whether or not he value should be returned as a log (TRUE) or not (FALSE)
    @return the maximal D-value'''
    
    # We need to find all the combinations (maximum of 16) when applying the rounding to the values in twoByTwoContingencyTable.  The exception is if
    # all four values are integers.  combinationsOfRounding will be a nx4 matrix containing the possible values of n_pp, n_pm, n_mp, n_mm
    if np.round(twoByTwoContingencyTable) == twoByTwoContingencyTable:
        # All four values are integers --> twoByTwoContingencyTable = [n_pp, n_pm, n_mp, n_mm] (c(...) in R)
        combinationsOfRounding = np.matrix(twoByTwoContingencyTable)#, ncol = 4)
    else:
        combinationsOfRounding = getAllPossibleRoundingCombinations(twoByTwoContingencyTable)

    maximumDValue = 0
    
    for i in range(combinationsOfRounding.shape[0]): #nrwo in R (loop over number of rows of the matrix)
        numbersOfCorrectandIncorrectPredictions = combinationsOfRounding[i] #get all elements of the i-th row ([i,] in R)
        threeByThreeContingencyTable = populateTheThreeByThreeContingencyTable(numbersOfCorrectandIncorrectPredictions[0, 0], numbersOfCorrectandIncorrectPredictions[0, 1], 
            numbersOfCorrectandIncorrectPredictions[0, 2], numbersOfCorrectandIncorrectPredictions[0,3], predictionListStats, experimentalDataStats)
        # Some of the combinationations produce an invalid contingency table - this we can check by seeing if any of the values in the table are negative.
        if threeByThreeContingencyTable.min() >= 0: #if 3x3ContTable is matrix or nparray! otherwise change to min(x)
            weight = calculateWeightGivenValuesInThreeByThreeContingencyTable(threeByThreeContingencyTable, logOfFactorialOfPredictionListStats, returnlog)
            #what type is weight? If not integer/float, max needs to be changed (e.g. for list max(max(weight), maximumDValue))
            maximumDValue = max(maximumDValue, weight)
    return maximumDValue