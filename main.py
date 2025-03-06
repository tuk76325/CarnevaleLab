import timeit   #times the speed of the algorithm
import Bio.AlignIO  #reading fasta file
import numpy
import numpy as np  #data processing
import pandas as pd #data processing
import matplotlib.pyplot as plt #visualization of data
from sklearn.preprocessing import OneHotEncoder #one hot encoding of amino acids
from sklearn.cluster import DBSCAN  #clustering method
from sklearn import metrics #counting number of clusters and outliers (noise points)
import re #to parse the sedIDs of the fasta file to remove the numbers
from sklearn.manifold import TSNE
#NEW COMMENT

def main():
#This block of code opens the text file so that python can read and parse the protein sequences
#     file = open(r'C:\Users\Toshan\PycharmProjects\BioinformaticsDoodnauth\PF00520_rp15.txt') #original fasta file
#     f = open("SeqPF00520_rp15.txt", 'w')
#     fileWrite = open("SeqPF00520_rp15.txt", 'r') #file with all sequences line by line
#     alignment = Bio.AlignIO.read(file, "fasta")
#     listSeq = []
#
# #This creates a list that has a sequence at each index, by this I mean that sequence 1 is the first position in the list, sequence 2 is the next position, and so on.
#     for record in alignment:
#         #print(record.seq)
#         f.write(str(record.seq) + "\n")
#         listSeq.append(str(record.seq).upper())
#     #print(listSeq[0])
#
#     aminoAcidList = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W', 'Y', 'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q', 'X', '-']
# #This block of code creates the binary arrays (vectors) for each character in the aminoAcidList
#     df = pd.DataFrame(aminoAcidList)
#     ohe = OneHotEncoder()
#     arraySequences = ohe.fit_transform(df).toarray()
#     #print(arraySequences)
#
# #This dictionary attaches the amino acid character from aminoAcidList to the corresponding binary vector
#     newDict = {}
#     for i in range(len(aminoAcidList)):
#         newDict.update({aminoAcidList[i] : arraySequences[i]})
#     #print(newDict)
#
# #This block of code currently appends each protein sequences' binary vector to a list called binaryList.
# #Where index 0 is every char in seq1 and index 1 is every char in se2, etc.
# #benchmark this code block = 3.26s
#     binaryList = []
#     tempList = []
#     for item in listSeq:
#         for character in item:
#             tempList.append((newDict[character]))
#         binaryList.append(tempList)
#         tempList = []
#     #print(binaryList[0])
#     #print(len(binaryList))
#
# #This block of code will compute the hamming distance via the matrices gotten from above by changing the 3D matrix to a 2D matrix via flattening
#     matrixCharBySeq = np.asarray(binaryList)
#     x = matrixCharBySeq
#     print(x.shape)
#     a= np.reshape(x, (20081, 1642*22))
#     print(a.shape)
#     print(a[1])
#
#     matrixSeqByChar = np.transpose(a)
#     print(matrixSeqByChar.shape)
#     matrixDot = np.dot(a, matrixSeqByChar)
#     matrixHammings = len(listSeq[0]) - matrixDot
#     print(matrixHammings.shape)
#
#     np.savetxt(f'matrixHammings.txt', matrixHammings)

    matrixHammings = numpy.loadtxt('matrixHammings.txt')

    # print(matrixHammings.flatten())
    # print(len(matrixHammings))

    # counts, bins = np.histogram(matrixHammings, bins= 100)
    #
    #
    # # print(matrixHammings)
    #
    # plt.stairs(counts, bins)
    # plt.show()

    clustering = DBSCAN(eps=205, min_samples=5, metric='precomputed').fit(matrixHammings)
    clusterLabels = clustering.labels_
    print(clustering)
    print(clusterLabels)
    print(type(clusterLabels))

    # Number of clusters in clusterLabels, ignoring noise if present.
    n_clusters_ = len(set(clusterLabels)) - (1 if -1 in clusterLabels else 0)
    n_noise_ = list(clusterLabels).count(-1)

    print("Estimated number of clusters: %d" % n_clusters_)
    print("Estimated number of noise points: %d" % n_noise_)

    # #visualization
    plt.scatter(matrixHammings[:, 0], matrixHammings[:,1], c = clusterLabels, vmin=0, cmap= "plasma") # plotting the clusters
    plt.xlabel("x") # X-axis label
    plt.ylabel("y") # Y-axis label
    plt.show() # showing the plot

#Uniprot to GO (tells you voltage gated postassium channel); Find FTP server with uniprot to GO via google (GO = Gene Ontology)
    #One cluster for every Uniprot code (the name of each sequence in the FASTA file ex. A0A2R2MPD5 for seq 1, 2, 3... etc.)

    # count = 0
    #
    # for i in clustering.labels_:
    #     if (i == -1):
    #         count+=1
    #
    # print(len(np.unique(clustering.labels_)))
    # print(count)

    # newF = open("testFileWrite.txt", 'w')
    # newF.write(str(matrixHammings))

startTime = timeit.default_timer()
main()
print("This algorithm takes: " + str(timeit.default_timer() - startTime) + " seconds")


def getIDs():
    newF = open(r'C:\Users\Toshan\PycharmProjects\BioinformaticsDoodnauth\PF00520_rp15.txt')
    fileIDSet = set()
    fileIDs = open("sequenceIDs.txt", 'w')
    for line in newF:
        if '>' in line:
            print((line))
            seqRE = re.findall('[^>/]+', line)
            fileIDSet.add((seqRE[0]) + "\n")
    print(fileIDSet)
    for line in fileIDSet:
        fileIDs.write(line)
    #upload text file to https://www.uniprot.org/id-mapping
    #https://regexr.com/

#startTime = timeit.default_timer()
#getIDs()
#print("This algorithm takes: " + str(timeit.default_timer() - startTime) + " seconds")

def csvData():
    df = pd.read_csv('Toshan_UniprotSeqIDs.csv')
    df = df.drop(df.index[:2600])
    df = df.drop(df.index[5201:])
    print(df)

def common_starting_string(str1, str2):
    common_string = ''
    for c1, c2 in zip(str1, str2):
        if c1 == c2:
            common_string += c1
        else:
            break
    return common_string.strip()

def group_descriptions(descriptions, min_common_chars=3):
    grouped_descriptions = []
    descriptions = sorted(descriptions)
    i = 0
    while i < len(descriptions) - 1:
        common_str = common_starting_string(descriptions[i], descriptions[i+1])
        if len(common_str) >= min_common_chars:
            grouped_descriptions.append(common_str)
            i += 2
        else:
            grouped_descriptions.append(descriptions[i])
            i += 1
    if i == len(descriptions) - 1:
        grouped_descriptions.append(descriptions[i])
    return grouped_descriptions

# # Read CSV file and parse descriptions
# filename = "UniprotSeqIDs.csv"
# df = pd.read_csv(filename)
#
# # Group similar descriptions
# descriptions = df['Protein names'].tolist()
# grouped_descriptions = group_descriptions(descriptions)
#
# # Count the occurrences of each grouped description
# description_counts = pd.Series(grouped_descriptions).value_counts()
#
# # Print the counted dictionary
# print(description_counts)
#
# description_counts.to_csv('description_counts.txt', header=False, sep='\t')
