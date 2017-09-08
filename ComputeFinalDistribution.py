# -*- coding: utf-8 -*-

import numpy

def compute_final_distribution(result_matrix):
    """ Computes a final reference distribution of the score used to compute the final p-value. 
    
    :param numpy.array result_matrix: a numpy matrix
    :rtype numpy.array
    :return distributionMatrix
    """
    
    if result_matrix.size == 2: 
        result_matrix = result_matrix.transpose()
        
    maxScore = max(result_matrix[:,0])
    minScore = min(result_matrix[:,0])
    
    #Pre-allocate the size of storage array
    distribution_matrix = zeros(shape=((maxScore - minScore + 1),2))
    numRows = result_matrix.shape[0]
    counter = 0
    for score in range(minScore, maxScore+1):
        probability = 0
        for i in range(numRows):
            if score == result_matrix[i,0]:
                probability = probability + result_matrix[i,1]
        if probability > 0:
            distribution_matrix[counter, ] = [score, probability]
            counter += 1
    
    #Remove the rows that were not populated in the for loops above
    distribution_matrix = distribution_matrix[0:counter,:]
    
    return distribution_matrix
