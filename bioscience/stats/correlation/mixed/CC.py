from bioscience.base import *
from math import comb

import math
import sys
import os
import threading
import warnings
import numpy as np
import time

def cc(dataset, deviceCount = 0, mode = 1, debug = False):
    """
    Application of the Contingency Coefficient (CC) method.
    
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
            oModel = __ccSequential(dataset, debug)
            deviceCount = 0
            sMode = "CPU Sequential"
    
    return oModel

def __ccSequential(dataset, debug):
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
            
            # Calculate the contingency table
            r1Unique = np.unique(dataset.data[r1])
            r2Unique = np.unique(dataset.data[r2])
            contingencyTable = np.zeros((len(r1Unique), len(r2Unique)))
            
            print(contingencyTable)
            
            # Calculate chi-squared
            rowTotals = np.sum(contingencyTable, axis=1)
            colTotals = np.sum(contingencyTable, axis=0)
            total = np.sum(contingencyTable)
            
            print(rowTotals)
            print(colTotals)
            print(total)
            
            expected = np.outer(rowTotals, colTotals) / total
            chi2 = np.sum((contingencyTable - expected)**2 / expected)
            
            cc = np.sqrt(chi2 / (chi2 + total))
            
            resultsCorrelation[pattern] = cc
                
        pattern += 1
    
    if debug == True:
        end_time = time.time()
        fExecutionTime = end_time - start_time
    
    oCorrelationResults = CorrelationModel(name=CC, results=resultsCorrelation, rows = iRows, executionTime=fExecutionTime)
    return oCorrelationResults