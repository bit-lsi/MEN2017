def find_approximate_values_that_will_maximise_D_value(predictionListStats, experimentalDataStats):
    """ Finds an approximate table values to maximise D.
    
    :param list predictionListStats: a list containing the values q+, q- and q0 which are numbers of positive, negative and non-significant/contradictory predictions
    :param list experimentalDataStats: a list containing the values n+, n- and n0 which are numbers of positive, negative and non-significant/contradictory predictions
    :rtype list
    :return twoByTwoContingencyTable: a list  which is a 2x2 contingency table which approximately maximises D
    """
    
    q_p = predictionListStats[0]
    q_m = predictionListStats[1]
    n_p = experimentalDataStats[0]
    n_m = experimentalDataStats[1]
    Tval = sum(predictionListStats)
    
    #The values of n++, n+-, n-+ and n-- that give the maximum D-value are given by the formula within the paper - 
    #Assessing statistical significance in causal graphs, page 6.  The formula is n_ab is approximately equal to q_a*n_b/T, where T = q+ + q- + q0 = n+ + n- + n0, and a,b are either + or -.
       
    n_pp = q_p * n_p/Tval
    n_pm = q_p * n_m/Tval
    n_mp = q_m * n_p/Tval
    n_mm = q_m * n_m/Tval
    twoByTwoContingencyTable = [n_pp, n_pm, n_mp, n_mm]
    
    return twoByTwoContingencyTable
    