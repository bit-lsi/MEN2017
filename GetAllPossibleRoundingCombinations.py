def getAllPossibleRoundingCombinations(twoByTwoContingencyTable):
    '''Helper function for getMaximumDValueFromTwoByTwoContingencyTable: get score for numbers of correct and incorrect predictions
    Returns all possible rounding combinations of a 2x2 table.  Given the values of n++, n+-, n-+ and n-- (stored in twoByTwoContingencyTable) this function will
    compute all possibilities of rounding each value up or down. 

    @param twoByTwoContingencyTable    Approximate values of n++, n+-, n-+ and n--, these values are
        calculated to optimise the D-value (see page 6 of Assessing statistical significance of causal graphs)
    @return a matrix of rounding combinations'''
    from math import floor, ceil
    import numpy as np
    
    f_n_pp <- floor(twoByTwoContingencyTable[0]) #(floor in R/python --> largest integer value less than or equal to x)
    f_n_pm <- floor(twoByTwoContingencyTable[1])
    f_n_mp <- floor(twoByTwoContingencyTable[2])
    f_n_mm <- floor(twoByTwoContingencyTable[3])
    
    c_n_pp <- ceil(twoByTwoContingencyTable[0]) #(ceil(ing) in R/python --> smallest integer value bigger than or equal to x)
    c_n_pm <- ceil(twoByTwoContingencyTable[1])
    c_n_mp <- ceil(twoByTwoContingencyTable[2])
    c_n_mm <- ceil(twoByTwoContingencyTable[3])
    
    roundingCombinations <- np.matrix([[f_n_pp, f_n_pm, f_n_mp, f_n_mm], [f_n_pp, f_n_pm, f_n_mp, c_n_mm], [f_n_pp, f_n_pm, c_n_mp, f_n_mm], [f_n_pp, f_n_pm, 
        c_n_mp, c_n_mm], [f_n_pp, c_n_pm, f_n_mp, f_n_mm], [f_n_pp, c_n_pm, f_n_mp, c_n_mm], [f_n_pp, c_n_pm, c_n_mp, f_n_mm], [f_n_pp, c_n_pm, c_n_mp, c_n_mm], 
        [c_n_pp, f_n_pm, f_n_mp, f_n_mm], [c_n_pp, f_n_pm, f_n_mp, c_n_mm], [c_n_pp, f_n_pm, c_n_mp, f_n_mm], [c_n_pp, f_n_pm, c_n_mp, c_n_mm], [c_n_pp, c_n_pm, 
        f_n_mp, f_n_mm], [c_n_pp, c_n_pm, f_n_mp, c_n_mm], [c_n_pp, c_n_pm, c_n_mp, f_n_mm], [c_n_pp, c_n_pm, c_n_mp, c_n_mm]] # in R, ncol = 4, byrow = TRUE
    
    
    # If any of the input values are integers then there will be duplicates in the list Remove the duplicates
    keep = [True for _ in range(16)]
    if f_n_pp == c_n_pp: 
        keep[8:16] = [False for _ in range(8)] #keep[9:16] <- FALSE
    if f_n_pm == c_n_pm:
        keep[4:8] = [False for _ in range(5)]
        keep[12:16] = [False for _ in range(5)] #keep[c(5:8, 13:16)] <- FALSE
    if f_n_mp == c_n_mp:
        keep[2:4] = [False for _ in range(3)]
        keep[6:8] = [False for _ in range(3)]
        keep[10:12] = [False for _ in range(3)]
        keep[14:16] = [False for _ in range(3)] #keep[c(3:4, 7:8, 11:12, 15:16)] <- FALSE
    if f_n_mm == c_n_mm:
        keep[1] = False
        keep[3] = False
        keep[5] = False
        keep[7] = False
        keep[9] = False
        keep[11] = False
        keep[3] = False
        keep[15] = False #keep[c(2, 4, 6, 8, 10, 12, 14, 16)] <- FALSE
    
    return roundingCombinations[keep, ]