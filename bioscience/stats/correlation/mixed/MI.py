from bioscience.base import *

import math
import sys
import os
import threading
import warnings
import numpy as np
import time

def mi(dataset, deviceCount = 0, mode = 1, debug = False):
    """
    Application of the Mutual Information (MI) correlation method.
    
    :param dataset: The dataset object to be binarized.
    :type dataset: :class:`bioscience.base.models.Dataset`
    
    param deviceCount: Number of GPU devices to execute
    :type deviceCount: int
    
    :param mode: Type of execution of the algorithm: `mode=1` for sequential execution, `mode=2` for parallel execution on CPUs and `mode=3` for execution on a multi-GPU architecture.
    :type mode: int
    
    :return: A CorrelationModel object that stores values generated by a correlation method.
    :rtype: :class:`bioscience.base.models.CorrelationModel`      
    """ 
    
    oModel = None
    if (dataset is not None):
        sMode = ""
        if mode == 2: # NUMBA: CPU Parallel mode
            # To be developed
            sMode = "NUMBA - CPU Parallel mode (to be developed)"
        elif mode == 3: # NUMBA: GPU Parallel mode
            # To be developed
            sMode = "NUMBA - GPU Parallel mode (to be developed)"
        else: # Sequential mode
            oModel = __miSequential(dataset, debug)
            deviceCount = 0
            sMode = "CPU Sequential"
    
    return oModel

def __miSequential(dataset, debug):
    iRows = dataset.data.shape[0]
    iCols = dataset.data.shape[1]
    fExecutionTime = None
    
    maxPairs = 0
    for i in range(iRows):
        for j in range(i + 1, iRows):
            maxPairs += 1
                    
    resultsCorrelation = np.zeros(maxPairs)
    
    if debug == True:
        start_time = time.time()
    
    pattern = 0
    while pattern < maxPairs:
        
        # Get R1 and R2 from index results vector
        r1 = 0
        r2 = -1
        auxPat = pattern - iRows + 1
        
        if auxPat < 0:
            r2 = auxPat + iRows

        j = iRows - 2
        while r2 == -1:
            auxPat -= j
            r1 += 1
            if auxPat < 0:
                r2 = (j + auxPat) + (r1 + 1)
            j -= 1
        
        if r1 < iRows and r2 < iRows:
            
            # Calculate the probability of R1
            probR1 = __calcProbability(dataset.data[r1])
            probR2 = __calcProbability(dataset.data[r2])
            
            # Calculate Joint probability
            combined = np.vstack((dataset.data[r1], dataset.data[r2])).T
            unique_pairs, counts = np.unique(combined, axis=0, return_counts=True)
            joint_probabilities = counts / len(dataset.data[r1])
            jointProb = dict(zip(map(tuple, unique_pairs), joint_probabilities))
            
            # Calculate mutual information
            mi = 0
            for (x, y), pR1R2 in jointProb.items():
                pR1 = probR1[x]
                pR2 = probR2[y]
                mi += pR1R2 * math.log(pR1R2 / (pR1 * pR2), 2)
            
            resultsCorrelation[pattern] = mi
                
        pattern += 1
    
    if debug == True:
        end_time = time.time()
        fExecutionTime = end_time - start_time
    
    oCorrelationResults = CorrelationModel(name=MI, results=resultsCorrelation, rows = iRows, executionTime=fExecutionTime)
    return oCorrelationResults

def __calcProbability(geneRow):
    
    unique, counts = np.unique(geneRow, return_counts=True)
    probabilities = counts / len(geneRow)
    return dict(zip(unique, probabilities))