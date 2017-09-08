
# coding: utf-8

# ## Implementation of CalculateSignificanceUsingCubicAlgorithm1b.r
# 
# ### Adapted from CausalR

import numpy as np
import math

def calculateWeightGivenValuesInThreeByThreeContingencyTable(threeByThreeContingencyTable, logOfFactorialOfPredictionListStats, returnlog):
    '''Helper Function for getMaximumDValueFromTwoByTwoContingencyTable'''

def populateTheThreeByThreeContingencyTable(numbersOfCorrectandIncorrectPredictions[1], numbersOfCorrectandIncorrectPredictions[2], 
            numbersOfCorrectandIncorrectPredictions[3], numbersOfCorrectandIncorrectPredictions[4], predictionListStats, experimentalDataStats):
    '''Helper Function for getMaximumDValueFromTwoByTwoContingencyTable'''

    pass



def findApproximateValuesThatWillMaximiseDValue(predictionListStats, experimentalDataStats):
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
 
    q_p = predictionListStats[0]
    q_m = predictionListStats[1]
    
    n_p = experimentalDataStats[0]
    n_m = experimentalDataStats[1]
    
    Tval = sum(predictionListStats)
    # The values of n++, n+-, n-+ and n-- that give the maximum D-value are given by the formula within the paper - Assessing statistical significance
    # in causal graphs, page 6.  The formula is n_ab is approximately equal to q_a*n_b/T, where T = q+ + q- + q0 = n+ + n- + n0, and a,b are either + or
    # -.
    n_pp = q_p * n_p/Tval
    n_pm = q_p * n_m/Tval
    n_mp = q_m * n_p/Tval
    n_mm = q_m * n_m/Tval
    
    twoByTwoContingencyTable = [n_pp, n_pm, n_mp, n_mm]
    
    return twoByTwoContingencyTable

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


def getWeightsAboveHypothesisScoreForAThreeByTwoTable(weights, r_p, r_m, r_z, n_p, n_m, predictionListStats, experimentalDataStats, 
                logOfFactorialOfPredictionListStats, hypothesisScore, logepsDMax, logDMax):
    '''Helper function for cubic1b: '''
    return [0, 0]

def cubic1b(hypScore, predictionListStats, expDataStats, epsilon):
    '''Calculate the p-value of a score given the hypothesis score and the 
    distribution table (calculated using the cubic algorithm  1b in 
    Assessing statistical significance in causal graphs - Chindelevitch et al) 
    '''
    # Calculate the log of the factorial for all values in predictionListStats,
    #this is to reduce the number of computations done later on.
    logOfFactorialOfPredictionListStats = [math.log(math.factorial(n)) for n in predictionListStats]
    logDMax = findMaximumDValue(predictionListStats, experimentalDataStats, logOfFactorialOfPredictionListStats, TRUE)
    
    # Number of positive predictions from the network. This will be an upper bound for r+ (r_p)
    q_p = predictionListStats[0] #c(up, down, ambig)
    # Number of negative predictions from the network. This will be an upper bound for r- (r_m)
    q_m = predictionListStats[1]
    
    # Number of positive results. This will be an upper bound for c+ (c_p)
    n_p = experimentalDataStats[0] #c(up, down, nochange)
    # Number of negative results. This will be an upper bound for c- (c_m)
    n_m = experimentalDataStats[1]
    
    # Array for the D values - first element is the total D values which have a score equal to or better than the cut-off, second element is all D
    # values (i.e. the ratio is the P value).
    logepsDMax = math.log(epsilon) + logDMax # or [n+math.log(epsilon) for n in logDMax] #
    weights = [0, 0]
    
    # In algorithm a superfamily is defined by the left 3x2 submatrix of the contingency table. Now we only need to interate over two quantities: r+ and
    # r-.
    for r_p in range(0, q_p):
        # for 0 to the limit of positive predictions
        for r_m in range(0, q_m):
            # for 0 to the limit of number of negative results Calculate r0
            r_z = n_p + n_m - (r_p + r_m)
            # Need to check that n00 is postive, to ensure it is a valid superfamily
            if r_z <= predictionListStats[2] and r_z >= 0: 
                # The five values in the array below define the superfamily
                weights = getWeightsAboveHypothesisScoreForAThreeByTwoTable(weights, r_p, r_m, r_z, n_p, n_m, predictionListStats, experimentalDataStats, 
                logOfFactorialOfPredictionListStats, hypothesisScore, logepsDMax, logDMax)
    
    pValue = weights[0]/weights[1]
    
    return(pValue)

    
    

    



