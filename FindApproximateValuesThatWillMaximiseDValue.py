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