import math

def get_max_D_value_for_a_three_by_two_family(r_p, r_m, r_z, n_p, n_m, predictionListStats, logOfFactorialOfPredictionListStats, returnlog = False):
    
    total = n_p + n_m
    
    #Compute the values of n++, n+-, n-+, n--, n0+, and n0- that maximise the D-value; these correspond to the following formula: 
    #n_ab is approximately equal to q_a*n_b/T, where T = n+ + n-, and a,b are either + or -.  
    #See Assessing statistical significance in causal graphs, page 7 - the formula is not stated explicitly but follows from the logic of algorithm 1a).
    
    if total > 0:
        n_pp = math.ceil(r_p * n_p/total)
        n_pm = r_p - n_pp
        n_mp = math.ceil(r_m * n_p/total)
    
    # Check this rounding produces a valid combination i.e. n++ + n-+ <= n+
    if (n_mp + n_pp) > n_p:
        n_mp = floor(r_m * n_p/total)
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
    
    maximumDFamValue = GetApproximateMaximumDValueFromThreeByTwoContingencyTable(threeByTwoContingencyTable, predictionListStats, logOfFactorialOfPredictionListStats, returnlog)
    
    return maximumDFamValue