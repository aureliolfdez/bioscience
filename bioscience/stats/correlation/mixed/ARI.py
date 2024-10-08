from bioscience.base import *
from math import comb

import math
import sys
import os
import threading
import warnings
import numpy as np
import time

def ari(dataset, deviceCount = 0, mode = 1, debug = False):
    """
    Application of the Adjusted Rand Index (ARI) method.
    
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
            oModel = __ariSequential(dataset, debug)
            deviceCount = 0
            sMode = "CPU Sequential"
    
    return oModel

def __ariSequential(dataset, debug):
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

            for i, lR1 in enumerate(r1Unique):
                for j, lR2 in enumerate(r2Unique):
                    contingencyTable[i, j] = np.sum((dataset.data[r1] == lR1) & (dataset.data[r2] == lR2))
            
            # Calculate ARI value
            n = int(np.sum(contingencyTable))
            sum_rows = np.sum(contingencyTable, axis=1)
            sum_cols = np.sum(contingencyTable, axis=0)
            comb_table = comb(n, 2)
            comb_rows = np.sum([comb(int(row), 2) for row in sum_rows])
            comb_cols = np.sum([comb(int(col), 2) for col in sum_cols])
            comb_cells = np.sum([comb(int(cell), 2) for cell in contingencyTable.flatten()])
            expected_index = (comb_rows * comb_cols) / comb_table
            max_index = (comb_rows + comb_cols) / 2
            ari = (comb_cells - expected_index) / (max_index - expected_index)
            
            resultsCorrelation[pattern] = ari
                
        pattern += 1
    
    if debug == True:
        end_time = time.time()
        fExecutionTime = end_time - start_time
    
    oCorrelationResults = CorrelationModel(name=ARI, results=resultsCorrelation, rows = iRows, executionTime=fExecutionTime)
    return oCorrelationResults