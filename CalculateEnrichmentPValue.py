# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import scipy.stats as stats


def data_frame_to_set(data_frame):
    return set([
        tuple(line)
        for line in data_frame.values.tolist()
    ])


# ' @title get set of significant predictions
# ' @description
# ' Gets the set of positive and negative predictions, the combination of the sets Sh+ and Sh- in Causal reasoning on biological networks: Interpreting transcriptional changes,  L Chindelevitch et al.
# ' @param predictions a table of predictions
# ' @return a matrix of positive and negative predictions

# ' @references
# ' L Chindelevitch et al.
# ' Causal reasoning on biological networks: Interpreting transcriptional changes.
# ' Bioinformatics, 28(8):1114-21, 2012.
def GetSetOfSignificantPredictions(predictions):
    counter = 1
    numPredictions = predictions.shape[0]
    significantPredictions = pd.DataFrame(0, index=numPredictions, columns=1)
    for i in range(1, numPredictions):
        if int(predictions.ix[i, 2] != 0):
            significantPredictions.ix[counter, 1] = predictions.ix[i, 1]
            counter = counter + 1

    significantPredictions = significantPredictions.loc[1:counter - 1]

    return significantPredictions


# ' @title get set of differientially expressed genes
# ' @description
# ' Gets the set of differentially expressed genes in the results, G+ as defined by in Causal reasoning on biological networks: Interpreting transcriptional changes,  L Chindelevitch et al.
# ' @param results a table of results
# ' @return a matrix of differentially expressed genes

# ' @references
# ' L Chindelevitch et al.
# ' Causal reasoning on biological networks: Interpreting transcriptional changes.
# ' Bioinformatics, 28(8):1114-21, 2012.
def GetSetOfDifferentiallyExpressedGenes(results):
    counter = 1
    numResults = results.shape[0]
    differentiallyExpressedGenes = pd.DataFrame(0, index=numResults, columns=1)
    for i in range(1, numResults):
        if int(results.ix[i, 2] != 0):
            differentiallyExpressedGenes.ix[counter, 1] = results[i, 1]
            counter = counter + 1

    differentiallyExpressedGenes = differentiallyExpressedGenes.loc[1:counter - 1]

    return differentiallyExpressedGenes


# ' calculates an enrichment p-value
# ' @description
# ' Calculate a enrichment p-value for a given hypothesis by comparing the corresponding predicted and observed gene changes

# ' @export
# ' @concept CausalR
# ' @param  predictions  predictions of changes from the CCG for a particular hypothesis
# ' @param  results  gene changes observed in the experimental data
# ' @return  an enrichment p-value
# ' @examples
# ' predictions <- matrix(c(1,2,3,1,1,-1), ncol = 2)
# ' results<- matrix(c(1,2,3,4,1,1,-1,1), ncol = 2)
def CalculateEnrichmentPValue(predictions, results):
    """

    :param predictions:
    :param pd.DataFrame results:
    :return:
    """

    setOfSignificantPredictions = data_frame_to_set(GetSetOfSignificantPredictions(predictions))
    setOfDifferentiallyExpressedGenes = data_frame_to_set(GetSetOfDifferentiallyExpressedGenes(results))

    setOfNonDifferentiallyExpressedGenes = results.ix[:, 1].difference(setOfDifferentiallyExpressedGenes)

    # n_pp + n_pm + n_mp + n_mm
    numSignificantPredictionsThatAreResponsive = len(
        setOfSignificantPredictions.intersection(setOfDifferentiallyExpressedGenes))
    # n+0 + n-0
    numSignificantPredictionsThatAreUnresponsive = len(
        setOfSignificantPredictions.intersection(setOfNonDifferentiallyExpressedGenes))
    # n0+ + n0-
    numZeroPredictionsThatAreResponsive = len(
        setOfDifferentiallyExpressedGenes) - numSignificantPredictionsThatAreResponsive
    # n00
    numZeroPredictionsThatAreUnreponsive = len(
        setOfNonDifferentiallyExpressedGenes) - numSignificantPredictionsThatAreUnresponsive

    contingencyTable = pd.DataFrame(
        np.r_[numSignificantPredictionsThatAreResponsive, numSignificantPredictionsThatAreUnresponsive,
              numZeroPredictionsThatAreResponsive,
              numZeroPredictionsThatAreUnreponsive], index=2)


    enrichmentPValue = stats.fisher_exact(contingencyTable, alternative="greater")

    print(enrichmentPValue)

    return enrichmentPValue


if __name__ == '__main__':
    predictions = pd.DataFrame(np.r_[1,2,3,1,1,-1], columns = 2)
    results = pd.DataFrame(np.r_[1,2,3,4,1,1,-1,1], columns = 2)
    CalculateEnrichmentPValue(predictions, results)