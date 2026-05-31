import numpy as np
import requests
from xml.etree import ElementTree

class GeneDataset:
    """
    This is a concept class representing a dataset.
    
    :param original: Here the original dataset is stored in a NumPy array.
    :type original: np.array
        
    :param data: This attribute is used to store the dataset once it has undergone the transformations desired by the user.
    :type data: np.array
         
    :param geneNames: Array with the name of the genes involved in the dataset. If the dataset does not have the name of the genes, it shall be replaced by a set of sequential numbers.
    :type geneNames: np.array, optional
    
    :param columnsNames: Array with the name of the columns involved in the dataset. If the dataset does not have the name of the columns, it shall be replaced by a set of sequential numbers.
    :type columnsNames: np.array, optional    
    
    :param annotations: Array that stores data from an annotation file for subsequent validation phases.
    :type annotations: np.array, optional
    
    :param cut: Cut-off parameter used in level binarisation.
    :type cut: float, optional
    
    """    
    
    def __init__(self, data, geneNames = None, columnsNames = None, annotations = None, cut = None):
        """
        Constructor method
        """
        self._original = data
        self._data = np.copy(self.original)
        self._geneNames = geneNames
        self._annotations = annotations
        self._cut = cut
        self._columnsNames = columnsNames
    
    @property
    def data(self):
        """
        Getter and setter methods of the data property.
        """
        return self._data
    
    @data.setter
    def data(self, data):
        self._data = data
    
    @property
    def geneNames(self):
        """
        Getter and setter methods of the geneNames property.
        """
        return self._geneNames
    
    @geneNames.setter
    def geneNames(self, geneNames):
        self._geneNames = geneNames
    
    @property
    def annotations(self):
        """
        Getter and setter methods of the annotations property.
        """
        return self._annotations
    
    @annotations.setter
    def annotations(self, annotations):
        self._annotations = annotations
    
    @property
    def cut(self):
        """
        Getter and setter methods of the cut property.
        """
        return self._cut
    
    @cut.setter
    def cut(self, cut):
        self._cut = cut
    
    @property
    def original(self):
        """
        Getter and setter methods of the original property.        
        """
        return self._original
    
    @original.setter
    def original(self, original):
        self._original = original
    
    @property
    def columnsNames(self):
        """
        Getter and setter methods of the columnsNames property.
        """
        return self._columnsNames
    
    @columnsNames.setter
    def columnsNames(self, columnsNames):
        self._columnsNames = columnsNames
    
    def __eq__(self, other):
        if isinstance(other, Dataset):
            return np.array_equal(self.data,other.data)
        else:
            return False
    
    def __hash__(self):
        return 1
    
    
class MicroarrayDataset(GeneDataset):
    """
    This is a concept class representing a microarray dataset.
    
    :param data: Matrix of expression values (intensities)
    :type data: np.array
    
    :param geneNames: Names of genes or probes
    :type geneNames: np.array, optional
    
    :param columnsNames: Names of samples
    :type columnsNames: np.array, optional
    
    :param annotations: Experimental annotations (e.g. control/treatment)
    :type annotations: np.array, optional
    
    :param cut: Cut-off for binarisation (optional, not commonly used)
    :type cut: float, optional
    
    :param probeNames: Names of probes (if different from geneNames)
    :type probeNames: np.array, optional
    """
    
    def __init__(self, data, geneNames=None, columnsNames=None, annotations=None, cut=None, probeNames=None):
        """
        Constructor method
        """
        # En Microarray NO usamos lengths, MIRAMOS GEN A GEN
        super().__init__(data, geneNames=geneNames, columnsNames=columnsNames, annotations=annotations, cut=cut)
        

class ScRNADataset(GeneDataset):
    """
    This is a concept class representing a single-cell RNA-Seq dataset.
    
    :param data: Matrix of gene expression counts
    :type data: np.array
    
    :param geneNames: Names of genes
    :type geneNames: np.array, optional
    
    :param columnsNames: Names of cells
    :type columnsNames: np.array, optional
    
    :param annotations: Cell annotations (cell type, cluster, etc.)
    :type annotations: np.array, optional
    
    :param cut: Cut-off for binarisation (not commonly used in scRNA-Seq)
    :type cut: float, optional
    """
    
    def __init__(self, data, geneNames=None, columnsNames=None, annotations=None, cut=None):
        """
        Constructor method
        """
        # En scRNA-Seq NO usamos lengths, MIRAMOS GEN A GEN
        super().__init__(data, geneNames=geneNames, columnsNames=columnsNames, lengths=None, annotations=annotations, cut=cut)


class BulkRNASeqDataset(GeneDataset):
    """
    This is a concept class representing a bulk RNA-Seq dataset.
    
    :param data: Matrix of gene expression counts
    :type data: np.array
    
    :param geneNames: Names of genes
    :type geneNames: np.array, optional
    
    :param columnsNames: Names of samples
    :type columnsNames: np.array, optional
    
    :param lengths: Array with gene length value (RNA-Seq)
    :type lengths: np.array, optional 
    
    :param lengths: Gene lengths (important for normalization)
    :type lengths: np.array, optional
    
    :param annotations: Sample annotations (e.g. control vs tumor)
    :type annotations: np.array, optional
    
    :param cut: Cut-off for binarisation
    :type cut: float, optional
    """
    
    def __init__(self, data, geneNames=None, columnsNames=None, lengths=None, annotations=None, cut=None):
        """
        Constructor method
        """
        # LENGTHS!!!!
        super().__init__(data, geneNames=geneNames, columnsNames=columnsNames, annotations=annotations, cut=cut)
        self._lengths = lengths

    @property
    def lengths(self):
        """
        Getter and setter methods of the lengths property.
        """
        return self._lengths
    
    @lengths.setter
    def lengths(self, lengths):
        self._lengths = lengths

class NetworkDataset(GeneDataset):
    """
    This is a concept class representing a network dataset.
    
    :param data: This attribute is used to store the dataset once it has undergone the transformations desired by the user.
    :type data: np.array
    
    :param geneNamesNodeA: Array with the name of the genes involved in the dataset. If the dataset does not have the name of the genes, it shall be replaced by a set of sequential numbers.
    :type geneNames: np.array, optional
    
    :param geneNamesNodeB: Array with the name of the genes involved in the dataset. If the dataset does not have the name of the genes, it shall be replaced by a set of sequential numbers.
    :type geneNames: np.array, optional
    
    :param columnsNames: Array with the name of the columns involved in the dataset. If the dataset does not have the name of the columns, it shall be replaced by a set of sequential numbers.
    :type columnsNames: np.array, optional
    
    :param extraInfo: Array that stores extra information about the network dataset.
    :type extraInfo: np.array, optional
    
    :param importantColumnsName: Name of the column that contains the most relevant information of the dataset.
    :type importantColumnsName: np.array, optional
    
    """    
    
    def __init__(self, data, geneNamesNodeA = None, geneNamesNodeB = None, columnsNames = None, extraInfo = None, importantColumnsName = None):
        """
        Constructor method
        """
        super().__init__(data, columnsNames = columnsNames)
        self._geneNamesNodeA = geneNamesNodeA
        self._geneNamesNodeB = geneNamesNodeB
        self._extraInfo = extraInfo
        self._importantColumnsName = importantColumnsName
    
    @property
    def geneNamesNodeA(self):
        """
        Getter and setter methods of the geneNamesNodeA property.
        """
        return self._geneNamesNodeA
    @geneNamesNodeA.setter
    def geneNamesNodeA(self, geneNamesNodeA):
        self._geneNamesNodeA = geneNamesNodeA

    @property
    def geneNamesNodeB(self):
        """
        Getter and setter methods of the geneNamesNodeB property.
        """
        return self._geneNamesNodeB
    
    @geneNamesNodeB.setter
    def geneNamesNodeB(self, geneNamesNodeB):
        self._geneNamesNodeB = geneNamesNodeB
        
    @property
    def extraInfo(self):
        """
        Getter and setter methods of the extraInfo property.
        """
        return self._extraInfo
    
    @extraInfo.setter
    def extraInfo(self, extraInfo):
        self._extraInfo = extraInfo

    @property
    def importantColumnsName(self):
        """
        Getter and setter methods of the importantColumnsName property.
        """
        return self._importantColumnsName
    
    @importantColumnsName.setter
    def importantColumnsName(self, importantColumnsName):
        self._importantColumnsName = importantColumnsName

class NCBIDataset:
    
    """
    This is a concept class representing a information of NCBI dataset.
    
    :param accessionNumber: Accession number of the dataset.
    :type accessionNumber: str
        
    :param data: This attribute is used to store the dataset once it has undergone the transformations desired by the user.
    :type data: np.array
         
    :param geneNames: Array with the name of the genes involved in the dataset. If the dataset does not have the name of the genes, it shall be replaced by a set of sequential numbers.
    :type geneNames: np.array, optional
    
    :param columnsNames: Array with the name of the columns involved in the dataset. If the dataset does not have the name of the columns, it shall be replaced by a set of sequential numbers.
    :type columnsNames: np.array, optional
    
    :param lengths: Array with gene length value (RNA-Seq)
    :type lengths: np.array, optional    
    
    :param annotations: Array that stores data from an annotation file for subsequent validation phases.
    :type annotations: np.array, optional
    
    :param cut: Cut-off parameter used in level binarisation.
    :type cut: float, optional
    
    """    
    
    def __init__(self, accessionNumber, title = None, summary = None, gpl = None, gse = None, taxonomy = None, gdstype = None, suppfile = None, nSamples = None, link = None, bioProject = None, samples = None):
        """
        Constructor method
        """
        self._accessionNumber = accessionNumber
        self._title = title
        self._summary = summary
        self._gpl = gpl
        self._gse = gse
        self._taxonomy = taxonomy
        self._gdstype = gdstype
        self._suppfile = suppfile
        self._nSamples = nSamples
        self._link = link
        self._bioProject = bioProject
        self._samples = samples
    
    @property
    def accessionNumber(self):
        """
        Getter and setter methods of the accessionNumber property.        
        """
        return self._accessionNumber
    
    @accessionNumber.setter
    def accessionNumber(self, accessionNumber):
        self._accessionNumber = accessionNumber
    
    @property
    def title(self):
        """
        Getter and setter methods of the title property.        
        """
        return self._title
    
    @title.setter
    def title(self, title):
        self._title = title
    
    @property
    def summary(self):
        """
        Getter and setter methods of the summary property.        
        """
        return self._summary
    
    @summary.setter
    def summary(self, summary):
        self._summary = summary
    
    @property
    def gpl(self):
        """
        Getter and setter methods of the gpl property.        
        """
        return self._gpl
    
    @gpl.setter
    def gpl(self, gpl):
        self._gpl = gpl
    
    @property
    def gse(self):
        """
        Getter and setter methods of the gse property.        
        """
        return self._gse
    
    @gse.setter
    def gse(self, gse):
        self._gse = gse
    
    @property
    def taxonomy(self):
        """
        Getter and setter methods of the taxonomy property.        
        """
        return self._taxonomy
    
    @taxonomy.setter
    def taxonomy(self, taxonomy):
        self._taxonomy = taxonomy
    
    @property
    def gdstype(self):
        """
        Getter and setter methods of the gdstype property.        
        """
        return self._gdstype
    
    @gdstype.setter
    def gdstype(self, gdstype):
        self._gdstype = gdstype
    
    @property
    def suppfile(self):
        """
        Getter and setter methods of the suppfile property.        
        """
        return self._suppfile
    
    @suppfile.setter
    def suppfile(self, suppfile):
        self._suppfile = suppfile
    
    @property
    def nSamples(self):
        """
        Getter and setter methods of the nSamples property.        
        """
        return self._nSamples
    
    @nSamples.setter
    def nSamples(self, nSamples):
        self._nSamples = nSamples
    
    @property
    def link(self):
        """
        Getter and setter methods of the link property.        
        """
        return self._link
    
    @link.setter
    def link(self, link):
        self._link = link
    
    @property
    def bioProject(self):
        """
        Getter and setter methods of the bioProject property.        
        """
        return self._bioProject
    
    @bioProject.setter
    def bioProject(self, bioProject):
        self._bioProject = bioProject
    
    @property
    def samples(self):
        """
        Getter and setter methods of the samples property.        
        """
        return self._samples
    
    @samples.setter
    def samples(self, samples):
        self._samples = samples
    
    def fullInfo(self):
        strReturn = ""
        strReturn += f"Accession number: {self.accessionNumber}"
        strReturn += f"Title: {self.title}"
        strReturn += f"GDS type: {self.gdstype}"        
        strReturn += f"Taxonomy: {self.taxonomy}"
        strReturn += f"Link: {self.link}"
        strReturn += f"BioProject: {self.bioProject}"        
        strReturn += f"Num. Samples: {self.nSamples}"
        strReturn += f"Samples:\n--------"
        
        
        
        strReturn += f"Summary:\n--------\n{self.summary}"
        return strReturn
        
    def __str__(self):
        strReturn = ""
        strReturn += f"Accession number: {self.accessionNumber}"
        strReturn += f"Title: {self.title}"
        strReturn += f"GDS type: {self.gdstype}"
        strReturn += f"Taxonomy: {self.taxonomy}"
        strReturn += f"Link: {self.link}"        
        strReturn += f"Num. Samples: {self.nSamples}"        
        return strReturn

class NCBIClient:    

    def __init__(self, idDB, apiKey = None):
        self._baseUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self._apiKey = apiKey
        self._idDB = idDB
        
    @property
    def apiKey(self):
        """
        Getter and setter methods of the apiKey property.
        """
        return self._apiKey
    
    @apiKey.setter
    def apiKey(self, apiKey):
        self._apiKey = apiKey
    
    @property
    def idDB(self):
        """
        Getter and setter methods of the idDB property.
        """
        return self._idDB
    
    @idDB.setter
    def idDB(self, idDB):
        self._idDB = idDB
    
    @property
    def baseUrl(self):
        """
        Getter and setter methods of the baseUrl property.
        """
        return self._baseUrl        

    def getIdsByGeo(self):
        
        if(self.apiKey is None):
                    
            params = {
                "db": "gds",
                "term": self.idDB,
                "retmode": "json",
                "retmax": 20000
            }
        
        else:
            
            params = {
                "db": "gds",
                "term": self.idDB,
                "retmode": "json",
                "retmax": 20000,
                "apiKey": self.apiKey
            }

        response = requests.get(self.baseUrl + "esearch.fcgi", params=params)

        matrixReturn = None
        if response.status_code == 200:
            data = response.json()      
            if data['esearchresult']['idlist']:
                matrixReturn = data['esearchresult']['idlist']
        
        return matrixReturn

    def getSummaryById(self, id):
        
        if(self.apiKey is None):
                    
            params = {
                "db": "gds",
                "id": id,
                "retmode": "json",
                "version": "2.0"
            }
        
        else:
            
            params = {
                "db": "gds",
                "id": id,
                "retmode": "json",
                "apiKey": self.apiKey,
                "version": "2.0"
            }

        response = requests.get(self._baseUrl + "esummary.fcgi", params=params)
        data = response.json()

        summaryGeo = None
        if data is not None:
                        
            try:
                accessionData = data['result'][id]['accession']
            except Exception:
                accessionData = None
            
            if accessionData is not None and accessionData == self.idDB:            
                summaryGeo = data['result'][id]
        
        return summaryGeo

    """
    def fetch_geo_data(self, geo_id):
        geo_record_id = self.search_geo(geo_id)

        if geo_record_id:
            params = {
                "db": "gds",
                "id": geo_record_id,
                "retmode": "xml",
                "email": self.EMAIL
            }

            response = requests.get(self.baseUrl + "efetch.fcgi", params=params)

            if response.status_code == 200:
                with open(f"{geo_id}_data.xml", "w") as file:
                    file.write(response.text)
                print(f"Datos guardados en {geo_id}_data.xml")

                return ElementTree.fromstring(response.content)
            else:
                print(f"Error al descargar los datos: {response.status_code}")
        else:
            print("No se pudo obtener el ID de GEO.")

    
    def fetch_geo_data_by_accession(self, accession_number):

        params = {
            "db": "gds",  
            "term": accession_number,  
            "retmode": "xml",  
            "email": self.EMAIL
        }

        response = requests.get(self.baseUrl + "esearch.fcgi", params=params)

        if response.status_code == 200:
            data = response.text
            if '<IdList>' in data:
                print(f"Accession Number {accession_number} encontrado.")
                
                start = data.find("<Id>")
                end = data.find("</Id>", start)
                geo_record_id = data[start + 4:end]
                
                print(f"ID de GEO: {geo_record_id}")

                return self.fetch_geo_data_by_id(geo_record_id)
            else:
                print(f"No se encontró información para el Accession Number: {accession_number}")
                return None
        else:
            print(f"Error en la solicitud: {response.status_code}")
            return None
    
    def fetch_geo_data_by_id(self, geo_record_id):
    
        params = {
            "db": "gds",
            "id": geo_record_id,
            "retmode": "xml",
            "email": self.EMAIL
        }

        response = requests.get(self.baseUrl + "efetch.fcgi", params=params)

        if response.status_code == 200:
            with open(f"{geo_record_id}_data.xml", "w") as file:
                file.write(response.text)
            print(f"Datos guardados en {geo_record_id}_data.xml")

            return ElementTree.fromstring(response.content)
        else:
            print(f"Error al descargar los datos: {response.status_code}")
            return None

    def parse_geo_data(self, xml_data):
        for docsum in xml_data.findall(".//DocSum"):
            title = docsum.find(".//Item[@Name='Title']").text
            print(f"Título: {title}")"""

class Validation:
    """
    This is a conceptual class representing a validation model after applying a data mining technique.
    
    :param measure: Name of the validation measure used.
    :type measure: str
        
    :param value: Value of the validation measure used.
    :type value: float
        
    """
    
    def __init__(self, measure, value):
        """
        Constructor method
        """
        self._measure = measure
        self._value = value
    
    @property
    def measure(self):
        """
        Getter and setter methods of the measure property.
        """
        return self._measure
    
    @measure.setter
    def measure(self, measure):
        self._measure = measure
    
    @property
    def value(self):
        """
        Getter and setter methods of the value property.
        """
        return self._value
    
    @value.setter
    def value(self, value):
        self._value = value
    
    def __eq__(self, other):
        if isinstance(other, Validation):
            return self.measure == other.measure
        else:
            return False
    
    def __hash__(self):
        return hash(self.measure)
    
    def __str__(self):
        return 'Measure (0): (1)'.format(self.measure, self.value)

class Bicluster:
    """
    This is a conceptual class representing a bicluster after applying a Biclustering technique.
    
    :param rows: Rows of the bicluster.
    :type rows: np.array
        
    :param cols: Columns of the bicluster.
    :type cols: np.array, optional
    
    :param data: Bicluster values according to the original dataset.
    :type data: np.array, optional
    
    :param validations: A set of instances from :class:`bioscience.base.models.Validation`.
    :type validations: np.array, optional
    
    """
    
    def __init__(self, rows, cols = None, data=None, validations=None):
        """
        Constructor method
        """
        self._rows = rows
        self._cols = cols
        self._data = data        
        self._validations = validations
    
    @property
    def rows(self):
        """
        Getter and setter methods of the rows property.
        """
        return self._rows
    
    @rows.setter
    def rows(self, rows):
        self._rows = rows
    
    @property
    def cols(self):
        """
        Getter and setter methods of the cols property.
        """
        return self._cols
    
    @cols.setter
    def cols(self, cols):
        self._cols = cols
    
    @property
    def data(self):
        """
        Getter and setter methods of the data property.
        """
        return self._data
    
    @data.setter
    def data(self, data):
        self._data = data
        
    @property
    def validations(self):
        """
        Getter and setter methods of the validations property.
        """
        return self._validations
    
    @validations.setter
    def validations(self, validations):
        self._validations = validations
        
    def sizeBicluster(self):
        """
        Number of total elements in the bicluster.
        """
        return len(self.rows) * len(self.cols)
    
    def sort(self):
        """ 
        Sort the column and row array by theirs indices.
        """
        self.rows.sort()
        self.cols.sort()
    
    def __str__(self):
        return 'Bicluster: rows{0} - cols{1}'.format(self.rows,self.cols)

class BiclusteringModel:
    
    """
    This is a conceptual class representing a set of biclusters generated after applying a Biclustering technique.
    
    :param results: Data structure (set) that stores all biclusters after running a Biclustering algorithm.
    :type results: set(:class:`bioscience.base.models.Bicluster`), optional
        
    :param executionTime: Time taken to execute the Biclustering method.
    :type executionTime: float, optional
    
    """
    
    def __init__(self, results = None):     
        """
        Constructor method
        """   
        if results is not None and all(isinstance(bic, Bicluster) for bic in results):
            self._results = results
        else:
            self._results = set()
        
        self._executionTime = 0
    
    @property
    def executionTime(self):
        """
        Getter and setter methods of the executionTime property.
        """
        return self._executionTime
    
    @executionTime.setter
    def executionTime(self, executionTime):
        self._executionTime = executionTime
    
    @property
    def results(self):
        """
        Getter and setter methods of the results property.
        """
        return self._results
    
    @results.setter
    def results(self, results):
        self._results = results
    
    def __str__(self):
        return '\n'.join(str(bic) for bic in self.results)

class CorrelationModel:
    
    """
    It is a conceptual class that represents the results generated by a correlation method.
    
    :param results: NumPy vector that stores the results of the row pairs of a dataset.
        
    :param executionTime: Time taken to execute the Correlation method.
    :type executionTime: float, optional
    
    """
    
    def __init__(self, name, results, rows, executionTime = None):   
        
        self._name = name
        self._results = results        
        self._executionTime = executionTime
        
        maxPairs = 0        
        for i in range(rows):
            for j in range(i + 1, rows):
                maxPairs += 1
        
        self._geneInteractionsIndex = np.zeros((maxPairs,2))
                
        pattern = 0
        while pattern < maxPairs:
            
            r1 = 0
            r2 = -1
            auxPat = pattern - rows + 1
            
            if auxPat < 0:
                r2 = auxPat + rows

            j = rows - 2
            while r2 == -1:
                auxPat -= j
                r1 += 1
                if auxPat < 0:
                    r2 = (j + auxPat) + (r1 + 1)
                j -= 1
            
            if(r1 < rows and r2 < rows):
                self._geneInteractionsIndex[pattern][0] = r1
                self._geneInteractionsIndex[pattern][1] = r2
            
            pattern += 1
    
    @property
    def name(self):
        """
        Getter and setter methods of the name property.
        """
        return self._name
    
    @name.setter
    def executionTime(self, name):
        self._name = name        
        
    @property
    def executionTime(self):
        """
        Getter and setter methods of the executionTime property.
        """
        return self._executionTime
    
    @executionTime.setter
    def executionTime(self, executionTime):
        self._executionTime = executionTime
    
    @property
    def results(self):
        """
        Getter and setter methods of the results property.
        """
        return self._results
    
    @results.setter
    def results(self, results):
        self._results = results
        
    @property
    def geneInteractionsIndex(self):
        """
        Getter and setter methods of the results property.
        """
        return self._geneInteractionsIndex
    
    @geneInteractionsIndex.setter
    def geneInteractionsIndex(self, geneInteractionsIndex):
        self._geneInteractionsIndex = geneInteractionsIndex
    

class Network:
    """
    This is a conceptual class representing a genetic network after applying a Genetic Network technique.
    
    :param node: List of the nodes of the network.
    :type node: np.array

    :param edge: List of the edges of the network.
    :type edge: np.array

    :param validations: A set of instances from :class:`bioscience.base.models.Validation`.
    :type validations: np.array, optional
    
    :param directed: Represents if a network is a directed graph or not.
    :type directed: boolean, optional

    """
    def __init__(self, nodeOrigin, nodeDestiny, edge, validations=None, directed=None):
        """
        Constructor method
        """
        self._nodesOrigin = nodeOrigin
        self._nodesDestiny = nodeDestiny
        self._edges = edge
        self._sizeNodesOrigin = len(self._nodesOrigin)
        self._sizeNodesDestiny = len(self._nodesDestiny)
        self._sizeEdges = len(self._edges)
        self._validations = validations
        self._directed = directed
    
    @property
    def nodesOrigin(self):
        """
        Getter and setter methods of the nodesOrigin property
        """
        return self._nodesOrigin

    @nodesOrigin.setter
    def nodesOrigin(self, nodesOrigin):
        self._nodesOrigin = nodesOrigin
    
    @property
    def nodesDestiny(self):
        """
        Getter and setter methods of the nodesDestiny property
        """
        return self._nodesDestiny

    @nodesDestiny.setter
    def nodesDestiny(self, nodesDestiny):
        self._nodesDestiny = nodesDestiny

    @property
    def edges(self):
        """
        Getter and setter methods of the edges property
        """
        return self._edges

    @edges.setter
    def edges(self, edges):
        self._edges = edges

    @property
    def sizeNodesOrigin(self):
        return self._sizeNodesOrigin

    @property
    def sizeNodesDestiny(self):
        return self._sizeNodesDestiny

    @property
    def sizeEdges(self):
        return self._sizeEdges

    @property
    def validations(self):
        """
        Getter and setter methods of the validation property
        """
        return self._validations
    
    @validations.setter
    def validations(self, validations):
        self._validations = validations

    @property
    def directed(self):
        return self._directed

    @directed.setter
    def directed(self, directed):
        self._directed = directed
    
    def shared_edges_count(self, gene):
        """
        Function to count the number of edges shared by a given gene (node).
        
        :param gene: The gene (node) to check for shared edges.
        :type gene: str or int
        :return: Number of edges shared by the given gene.
        :rtype: int
        """
        count = 0
        for edge in self.edges:
            if edge.nodeA == gene or edge.nodeB == gene:
                count += 1
        return count

    def __str__(self):
        """
        String representation of the Network object.
        Provides a concise summary of nodes, edges, and graph properties.
        """
        summary = []
        summary.append("Network Summary:")
        summary.append(f"  Directed: {self.directed}")
        summary.append(f"  Nodes (Origin): {self.sizeNodesOrigin}")
        summary.append(f"  Nodes (Destiny): {self.sizeNodesDestiny}")
        summary.append(f"  Edges: {self.sizeEdges}")
        
        # Ejemplo de nodos origen/destino
        if self._sizeNodesOrigin > 0:
            sample_origin = ", ".join(map(str, self.nodesOrigin[:5]))
            summary.append(f"  Origin Nodes: {sample_origin} {'...' if self._sizeNodesOrigin > 5 else ''}")
        if self._sizeNodesDestiny > 0:
            sample_destiny = ", ".join(map(str, self.nodesDestiny[:5]))
            summary.append(f"  Destiny Nodes: {sample_destiny} {'...' if self._sizeNodesDestiny > 5 else ''}")
        
        # Ejemplo de aristas
        if self._sizeEdges > 0:
            try:
                sample_edges = ", ".join(str(e) for e in self.edges[:3])
            except Exception:
                sample_edges = ", ".join(f"({e.nodeA}, {e.nodeB})" for e in self.edges[:3])
            summary.append(f"  Edges: {sample_edges} {'...' if self._sizeEdges > 3 else ''}")
        
        return "\n".join(summary)

    

class NetworkModel:

    """
    This is a conceptual class representing a set of networks generated after applying a Network technique.
    
    :param results: Data structure (set) that stores all networks after running a Network algorithm.
    :type results: set(:class:`bioscience.base.models.Network`), optional
        
    :param executionTime: Time taken to execute the Network method.
    :type executionTime: float, optional
    
    """

    def __init__(self, results = None):     
        """
        Constructor method
        """   
        if results is not None and all(isinstance(net, Network) for net in results):
            self._results = results
        else:
            self._results = set()
        
        self._executionTime = 0
        self._decrease = 0
    
    @property
    def executionTime(self):
        """
        Getter and setter methods of the executionTime property.
        """
        return self._executionTime
    
    @executionTime.setter
    def executionTime(self, executionTime):
        self._executionTime = executionTime
    
    @property
    def results(self):
        """
        Getter and setter methods of the results property.
        """
        return self._results

    @results.setter
    def results(self, results):
        self._results = results

    def sort(self):
        self._results = sorted(self._results, key=lambda net: (len(net.sizeNodes), len(net.sizeEdges)))

    def __str__(self):
        n = len(self.results) if self.results is not None else 0
        if n == 0:
            return f"NetworkModel(empty, executionTime={self.executionTime}s)"

        NetworkModel_list = list(self.results)

        lines = [f"NetworkModel with {n} networks (executionTime={self.executionTime}s):"]
        for i, net in enumerate(NetworkModel_list, 1):
            lines.append(f"[{i}] {str(net)}")
        return "\n".join(lines)

class Node:

    """
    This is a conceptual class represented a node of a genetic network.

    :param info: Aditional information of the node.
    :type info: np.array

    :param id: Id of the gene.
    :type id: integer, optional.

    :param name: Name of the gene.
    :type name: string, optional.

    """
    
    def __init__(self, id = None, name = None):
        """
        Constructor Method
        """
        self._id = id
        self._name = name
    
    @property
    def id(self):
        """
        Getter and setter methods of the id property.
        """
        return self._id

    @id.setter
    def id(self, id):
        self._id = id
    
    @property
    def name(self):
        """
        Getter and setter methods of the name property.
        """
        return self._name
    
    @name.setter
    def name(self, name):
        self._name = name
        
    def __str__(self):
        return f"Node(id={self.id}, name={self.name})"

class Edge:
    """
    This is a conceptual class represented a edge of a genetic network.

    :param nodeA: Beginning node of the edge
    :type nodeA: :class:`bioscience.base.models.Node`
    
    :param nodeB: End node of the edge
    :type nodeB: :class:`bioscience.base.models.Node`

    :param weight: Weight of the edge.
    :type weight: float
    
    :param weightRelatedValues: Values related to the weight such as Spearman or Kendall
    :param weightRelatedValues: np.array

    :param info: Aditional information of the edge.
    :type info: np.array

    """
    def __init__(self, nodeA, nodeB, weight, weightRelatedValues, info):
        """
        Constructor method
        """
        self._nodeA = nodeA
        self._nodeB = nodeB
        self._weight = weight
        self._weightRelatedValues = weightRelatedValues
        self._info = info
    
    @property
    def nodeA(self):
        """
        Getter and setter methods of the nodeA property.
        """
        return self._nodeA

    @nodeA.setter
    def nodeA(self, nodeA):
        self._nodeA = nodeA

    @property
    def nodeB(self):
        """
        Getter and setter methods of the nodeB property.
        """
        return self._nodeB

    @nodeB.setter
    def nodeB(self, nodeB):
        self._nodeB = nodeB

    @property
    def weight(self):
        """
        Getter and setter methods of the weight property.
        """
        return self._weight

    @weight.setter
    def weight(self, weight):
        self._weight = weight
    
    @property
    def weightRelatedValues(self):
        """
        Getter and setter methods of the weightRelatedValues
        """
        return self._weightRelatedValues

    @weightRelatedValues.setter
    def weightRelatedValues(self, weightRelatedValues):
        self._weightRelatedValues = weightRelatedValues
    
    @property
    def info(self):
        """
        Getter and setter methods of the extraCharacteristics property.
        """
        return self._info
    
    @info.setter
    def info(self, info):
        self._info = info
    
    def __str__(self):
        return (
            f"Edge({self._nodeA} → {self._nodeB}, "
            f"weight={self._weight}, "
            f"relatedValues={self._weightRelatedValues}, "
            f"info={self._info})"
        )