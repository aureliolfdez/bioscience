from bioscience.base import *

import math
import sys
import os
import threading
import warnings
import numpy as np
import time

def quadrant(dataset, deviceCount = 0, mode = 1, debug = False):
    """
    Application of the quadrant (Q) metric.
    
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
            oModel = __quadrantSequential(dataset, debug)
            deviceCount = 0
            sMode = "CPU Sequential"
    
    return oModel

def __quadrantSequential(dataset, debug):
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
        
    for pattern, value in enumerate(resultsCorrelation):
        
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
            
            meanR1 = 0.0
            for value in dataset.data[r1]:
                meanR1 += value
            meanR1 /= iCols
            
            meanR2 = 0.0
            for value in dataset.data[r2]:
                meanR2 += value
            meanR2 /= iCols
            
            R1centered = dataset.data[r1] - meanR1
            R2centered = dataset.data[r2] - meanR2
            
            # Count points in each quadrant
            Q1 = np.sum((R1centered > 0) & (R2centered > 0))
            Q2 = np.sum((R1centered < 0) & (R2centered > 0))
            Q3 = np.sum((R1centered < 0) & (R2centered < 0))
            Q4 = np.sum((R1centered > 0) & (R2centered < 0))
            
            # Calculate Q measure
            qValue = (Q1 + Q3 - Q2 - Q4) / iCols
        
            resultsCorrelation[pattern] = qValue
    
    if debug == True:
        end_time = time.time()
        fExecutionTime = end_time - start_time
    
    oCorrelationResults = CorrelationModel(name=QUADRANT, results=resultsCorrelation, rows = iRows, executionTime=fExecutionTime)
    return oCorrelationResults