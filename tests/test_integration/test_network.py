import bioscience as bs

###################
# 1) Load dataset 
###################
dataset = bs.loadNetwork(db="/Users/aurelio/git-repositories/bioscience/datasets/network1.csv", index_nodeA = 1, index_nodeB = 3, index_weight = 9, separator = ",", skipr = 0, head = 0)

###################
# 2) Checking that everything is saved correctly
###################
#print("Dataset data column 0:")
#print(dataset.data[0])

# print("Dataset data column 1:")
# print(dataset.data[1])

# print("Dataset data column 2:")
# print(dataset.data[2])

# print("Dataset Extra Info:")
# print(dataset.extraInfo[0])

# print("Gene Names")
# print(dataset.geneNames)

# print("Columns Names")
# print(dataset.columnsNames)

# print("Important Columns Name")
# print(dataset.importantColumnsName)

# print("Gene Names Node A")
# print(dataset.geneNamesNodeA[0])

# print("Gene Names Node B")
# print(dataset.geneNamesNodeB[0])


print(dataset.extraInfo)