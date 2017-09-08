import numpy

def computeFinalDistribution(resultsMatrix):
    """ Computes a final reference distribution of the score used to compute the final p-value. 
    
    :param numpy.array resultsMatrix: a numpy matrix
    :rtype numpy.array
    :return distributionMatrix
    """
    
    if resultsMatrix.size == 2: 
        resultsMatrix = resultsMatrix.transpose()
        
    maxScore = max(resultsMatrix[:,0])
    minScore = min(resultsMatrix[:,0])
    
    #Pre-allocate the size of storage array
    distributionMatrix = zeros(shape=((maxScore - minScore + 1),2))
    numRows = resultsMatrix.shape[0]
    counter = 0
    for score in range(minScore, maxScore+1):
        probability = 0
        for i in range(numRows):
            if score == resultsMatrix[i,0]:
                probability = probability + resultsMatrix[i,1]
        if probability > 0:
            distributionMatrix[counter, ] = [score, probability]
            counter += 1
    
    #Remove the rows that were not populated in the for loops above
    distributionMatrix = distributionMatrix[0:counter,:]
    
    return distributionMatrix