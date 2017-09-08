# -*- coding: utf-8 -*-

import math

def calculateSignificanceUsingCubicAlgorithm1b(hypScore, predictionListStats, expDataStats, epsilon):
    '''Calculate the p-value of a score given the hypothesis score and the 
    distribution table (calculated using the cubic algorithm  1b in 
    Assessing statistical significance in causal graphs - Chindelevitch et al) 
    '''
    # Calculate the log of the factorial for all values in predictionListStats,
    #this is to reduce the number of computations done later on.
    logOfFactorialOfPredictionListStats = [math.log(math.factorial(n)) for n in predictionListStats]
    logDMax = findMaximumDValue(predictionListStats, experimentalDataStats, logOfFactorialOfPredictionListStats, TRUE)

    # Number of positive predictions from the network. This will be an upper bound for r+ (r_p)
    q_p, q_m = predictionListStats[0:2]

    # Number of negative predictions from the network. This will be an upper bound for r- (r_m)
    q_m = predictionListStats[1]

    # Number of positive results. This will be an upper bound for c+ (c_p)
    # Number of negative results. This will be an upper bound for c- (c_m)
    n_p, n_m = experimentalDataStats[0:2]


    
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
    
    return pValue