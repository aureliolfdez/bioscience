from operator import index
import numpy as np
import pandas as pd
import os
import time
import ftplib
from urllib.parse import urlparse

from .models import *

def load(db, apiKey = None, separator = "\t", skipr = 0, naFilter = False, index_gene = -1, index_lengths = -1, head = None) -> pd.DataFrame:
    """
    Load any data from a txt or csv file. (Reuse function)
    
    :param db: The path where the file is stored or ID database
    :type db: str
    
    :param apiKey: API Key NCBI
    :type apiKey: str
        
    :param separator: An attribute indicating how the columns of the file are separated.
    :type separator: str,optional
         
    :param skipr: Number of rows the user wishes to omit from the file, defaults to 0.
    :type skipr: int, optional
    
    :param naFilter: Boolean to detect NA values in a file. NA values shall be replaced by 0's, defaults to False.
    :type naFilter: boolean, optional
    
    :param index_gene: Column position where the gene names are stored in the dataset, defaults to -1 (deactivated).
    :type index_gene: int, optional
    
    :param index_lengths: Column position where the gene lengths are store in the dataset, defaults to -1 (deactivated).
    :type index_lengths: int, optional
    
    :param head: Row number(s) containing column labels and marking the start of the data (zero-indexed), defaults to None.
    :type head: int, optional
        
    :return: A dataset object from the reading of a file
    :rtype: :class:`bioscience.base.models.Dataset`
    """
    dataset = None
    if db is not None:
        extensionsCsv = [".txt",".csv",".tsv"]    
        fileName, fileExtension = os.path.splitext(db)
        if fileExtension in extensionsCsv:
            return __loadSpecificFile(db, separator, skipr, naFilter, index_gene, index_lengths, head)
        else:
            if db is not None:
                __loadNcbiDb(db, key = apiKey)
            else:
                return None

def __loadSpecificFile(db, separator, skipr, naFilter, index_gene, index_lengths, head) -> pd.DataFrame:
    if naFilter is True:
        dfPandas = pd.read_csv(db, sep=separator, skiprows = skipr, na_filter=naFilter, header = head).fillna(0)
    else:
        dfPandas = pd.read_csv(db, sep=separator, skiprows = skipr, header = head)
            
    dataColumns = np.asarray(dfPandas.columns)
    dataset = np.asarray(dfPandas)
    # Dataset object
    geneNames = None
    lengths = None
            
    # Get gene names from dataset
    if index_gene >= 0:
        index_lengths = index_lengths -1 # Update gene lengths index due to remove the gene column of the dataset.
        geneNames = dataset[:,index_gene]
        dataset = np.delete(dataset, index_gene, 1)                
            
    # Get lengths from dataset
    if index_lengths >= 0:
        lengths = dataset[:,index_lengths]
        dataset = np.delete(dataset, index_lengths, 1)                
            
    dataColumns = np.delete(dataColumns, np.arange(0, 1 + index_gene))         
    return Dataset(dataset.astype(np.double), geneNames=geneNames, columnsNames=dataColumns, lengths=lengths)

def __loadNcbiDb(idGeo, key = None, fileType = "raw"):
    
    # 1) Database connection 
    client, idsByGeo = __connectDatabase(idGeo, key)
    if idsByGeo == None:
        print("Unable to connect to the database.")
        return None
        
    print("Database connected.")
    
    # 2) Get info database
    infoDataset = __getInfoDatabase(client, idGeo, idsByGeo)
    
    # 3) Download dataset
    __downloadDatabase(infoDataset)
    
    #xml_data = client.fetch_geo_data(geo_id)
    
    #if xml_data:
    #    client.parse_geo_data(xml_data)

def __connectDatabase(idGeo, key):
    client = NCBIClient(idDB = idGeo, apiKey = key)
    print("Connecting NCBI database...")
    idsByGeo = client.getIdsByGeo()
    time.sleep(1)
    return client, idsByGeo

def __getInfoDatabase(client, idGeo, idsByGeo):
    
    print("Getting information from the ",idGeo," database...")
    
    summaryGeo = None
    iCount = 0
    idValid = -1
    while (iCount < len(idsByGeo) and summaryGeo == None):
        summaryGeo = client.getSummaryById(idsByGeo[iCount])
        if summaryGeo == None:
            time.sleep(1)
        else:
            idValid = idsByGeo[iCount]
        iCount += 1
    
    infoDataset = None
    if summaryGeo is None:
        print("ERROR: Information not obtained.")
    else:
        print("Information obtained.")
        try:
            sAccession = summaryGeo['accession']
        except Exception:
            sAccession = None
            
        try:
            sTitle = summaryGeo['title']
        except Exception:
            sTitle = None
        
        try:
            sSummary = summaryGeo['summary']
        except Exception:
            sSummary = None
            
        try:
            sGPL = summaryGeo['gpl']
        except Exception:
            sGPL = None
        
        try:
            sGSE = summaryGeo['gse']
        except Exception:
            sGSE = None
        
        try:
            sTaxon = summaryGeo['taxon']
        except Exception:
            sTaxon = None
        
        try:
            sGdsType = summaryGeo['gdstype']
        except Exception:
            sGdsType = None
        
        try:
            sSuppFile = summaryGeo['suppfile']
        except Exception:
            sSuppFile = None
        
        try:
            sNSamples = summaryGeo['n_samples']
        except Exception:
            sNSamples = None
        
        try:
            sFTPLink = summaryGeo['ftplink']
        except Exception:
            sFTPLink = None
        
        try:
            sBioProject = summaryGeo['bioproject']
        except Exception:
            sBioProject = None
        
        try:
            aSamples = np.array([(item['accession'],item['title']) for item in summaryGeo['samples']])
        except Exception:
            aSamples = None
             
        infoDataset = NCBIDataset(accessionNumber = sAccession, title = sTitle, summary = sSummary, gpl = sGPL, gse = sGSE, taxonomy = sTaxon, gdstype = sGdsType, suppfile = sSuppFile, nSamples = sNSamples, link = sFTPLink, bioProject = sBioProject, samples = aSamples)
        return infoDataset

def __downloadDatabase(infoDataset):
    if infoDataset is not None:
        parsedUrl = urlparse(infoDataset.link)
        ftpHostname = parsedUrl.hostname
        ftpPath = parsedUrl.path
        
        ftp = ftplib.FTP(ftpHostname)
        ftp.login()
        ftp.cwd(ftpPath)
        
        # Files list
        filesList = ftp.nlst()
        for file in filesList:
            
            if file == 'matrix':
                __downloadMatrix(ftp,file)
            elif file == 'soft':
                __downloadSoft(ftp,file)
            elif file == 'suppl':
                __downloadSuppl(ftp,file)
        
        ftp.quit()

def __downloadMatrix(ftp,file):
    print("MATRIX")
    ftp.cwd(file)
    filesList = ftp.nlst()
    for file in filesList:
        if '.tar' in file or '.gz' in file or '.zip' in file:            
            print(file)
    ftp.cwd('..')
    
def __downloadSoft(ftp,file):
    print("SOFT")
    ftp.cwd(file)
    filesList = ftp.nlst()
    for file in filesList:
        if '.tar' in file or '.gz' in file or '.zip' in file:
            print(file)
    ftp.cwd('..')
    
def __downloadSuppl(ftp,file):
    print("SUPPL")
    ftp.cwd(file)
    filesList = ftp.nlst()
    for file in filesList:
        if '.tar' in file or '.gz' in file or '.zip' in file:
            print(file)
    ftp.cwd('..')

def __downloadGSE(infoDataset, fileType):
    
    parsedUrl = urlparse(infoDataset.link)
    ftpHostname = parsedUrl.hostname
    ftpPath = parsedUrl.path
        
    ftp = ftplib.FTP(ftpHostname)
    ftp.login()
    ftp.cwd(ftpPath)
        
    # Files list
    filesList = ftp.nlst()
    for file in filesList:
        print(file)
        
    ftp.quit()
    
def saveResultsIndex(path, models):    
    """
    Save the results index (rows and columns index of the dataset) of applying a data mining technique.
    
    :param path: The path where the file will be stored.
    :type path: str
        
    :param models: An attribute indicating how the columns of the file are separated.
    :type models: :class:`bioscience.dataMining.biclustering.BiclusteringModel`
         
    """
    if models is not None:
        iLevel = 1     
        for model in models:
            infoData = []
            for oBicluster in model.results:
                if oBicluster.rows is not None:
                    rows = ','.join(str(int(row)) for row in oBicluster.rows)
                else:
                    rows = ""
                
                if oBicluster.cols is not None:
                    cols = ','.join(str(int(col)) for col in oBicluster.cols)
                else:
                    cols = ""
                
                infoData.append(rows + ';' + cols)
            
            df = pd.DataFrame(infoData, columns=['Data'])
            df['Data'] = df['Data'].str.replace('"', '')
            df.to_csv(path+"index"+str(iLevel)+".csv", index=False, header=False)
            iLevel += 1
    
        print("Results index saved in: " + path)

def saveResults(path, models, data):
    """
    Save the results of applying a data mining technique.
    
    :param path: The path where the file will be stored.
    :type path: str
        
    :param models: The results of the data mining technique.
    :type models: :class:`bioscience.dataMining.biclustering.BiclusteringModel`
    
    :param data: The dataset object which stores the original dataset.
    :type data: :class:`bioscience.base.models.Dataset`
         
    """
    if models is not None:
        iLevel = 1
        for i, model in enumerate(models, start=1):
            infoModel = ""
            for j, oBicluster in enumerate(model.results, start=1):
                if isinstance(data, set):
                    geneNames = list(data)[iLevel-1].geneNames
                    colNames = list(data)[iLevel-1].columnsNames
                else:
                    geneNames = data.geneNames
                    colNames = data.columnsNames
                
                if oBicluster.rows is not None:
                    if geneNames is not None:
                        rows = ','.join(str(geneNames[int(row)]) for row in oBicluster.rows)
                    else:
                        rows = ','.join(str(row) for row in oBicluster.rows)
                else:
                    rows = ""
                
                if oBicluster.cols is not None:
                    if colNames is not None:
                        cols = ','.join(str(colNames[int(col)]) for col in oBicluster.cols)
                    else:
                        cols = ','.join(str(int(col)) for col in oBicluster.cols)
                else:
                    cols = ""
                
                infoBicluster = f"\nRESULT #{j} (ROWS: {rows}) - (COLS: {cols})\n"
                   
                if isinstance(data, set):
                    dataset = list(data)[i - 1].original
                else:
                    dataset = data.original
                for oRow in oBicluster.rows:
                    if oBicluster.cols is not None:
                        infoBicluster += ",".join(str(dataset[int(oRow)][int(oCol)]) for oCol in oBicluster.cols)
                        infoBicluster += "\n"

                infoModel += infoBicluster

            df = pd.DataFrame([infoModel], columns=['Data'])
            df['Data'] = df['Data'].str.replace('"', '')
            df.to_csv(f"{path}results{i}.csv", index=False, header=False)
            iLevel += 1

        print("Results saved in: " + path)
        
def saveGenes(path, models, data):
    """
    Save the gene names from the results of applying a data mining technique.
    
    :param path: The path where the file will be stored.
    :type path: str
        
    :param models: The results of the data mining technique.
    :type models: :class:`bioscience.dataMining.biclustering.BiclusteringModel`
    
    :param data: The dataset object which stores the original dataset.
    :type data: :class:`bioscience.base.models.Dataset`
         
    """
    if models is not None:
        iLevel = 1
        for model in models:
            if isinstance(data, set):
                geneNames = list(data)[iLevel-1].geneNames
            else:
                geneNames = data.geneNames
            
            infoData = []
            for oBicluster in model.results:
                if geneNames is not None:
                    rows = ','.join(str(geneNames[int(row)]) for row in oBicluster.rows)
                else:
                    rows = ','.join(str(row) for row in oBicluster.rows)
                infoData.append(rows)
            
            df = pd.DataFrame(infoData, columns=['Data'])
            df['Data'] = df['Data'].str.replace('"', '')
            df.to_csv(path+"genes"+str(iLevel)+".csv", index=False, header=False)
            iLevel += 1

    print("Genes saved in: " + path)

def saveBinaryDatasets(path, datasets):
    """
    If the dataset has been binarised, this function allows storing the binary dataset.
    
    :param path: The path where the file will be stored.
    :type path: str
        
    :param datasets: The dataset object which stores the binary dataset.
    :type datasets: :class:`bioscience.base.models.Dataset`
         
    """
    if datasets is not None:
        if isinstance(datasets, set):
            iLevel = 1
            for dataset in datasets:
                df = pd.DataFrame(dataset.data)
                df.to_csv(path+"dataset"+str(iLevel)+".csv", index=False, header=False)
                iLevel += 1
        else:
            df = pd.DataFrame(datasets.data)
            df.to_csv(path+"dataset.csv", index=False, header=False)
        
        print("Binary datasets saved in: " + path)


        
        
        
    
    