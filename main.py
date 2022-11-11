import timeit
import time
import Bio.AlignIO
import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder

def main():
#This block of code opens the text file so that python can read and parse the protein sequences
    file = open(r'C:\Users\Toshan\PycharmProjects\BioinformaticsDoodnauth\PF00520_rp15.txt') #original fasta file
    f = open("SeqPF00520_rp15.txt", 'w')
    fileWrite = open("SeqPF00520_rp15.txt", 'r') #file with all sequences line by line
    alignment = Bio.AlignIO.read(file, "fasta")
    listSeq = []

#This creates a list that has a sequence at each index, by this I mean that sequence 1 is the first position in the list, sequence 2 is the next position, and so on.
    for record in alignment:
        #print(record.seq)
        f.write(str(record.seq) + "\n")
        listSeq.append(str(record.seq).upper())
    #print(listSeq[0])

    aminoAcidList = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W', 'Y', 'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q', 'X', '-']
#This block of code creates the binary arrays (vectors) for each character in the aminoAcidList
    df = pd.DataFrame(aminoAcidList)
    ohe = OneHotEncoder()
    arraySequences = ohe.fit_transform(df).toarray()
    #print(arraySequences)

#This dictionary attaches the amino acid character from aminoAcidList to the corresponding binary vector
    newDict = {}
    for i in range(len(aminoAcidList)):
        newDict.update({aminoAcidList[i] : arraySequences[i]})
    #print(newDict)

#This block of code currently appends each protein sequences' binary vector to a list called binaryList.
#Where index 0 is every char in seq1 and index 1 is every char in se2, etc.
#benchmark this code block = 3.26s
    binaryList = []
    tempList = []
    for item in listSeq:
        for character in item:
            tempList.append((newDict[character]))
        binaryList.append(tempList)
        tempList = []
    #print(binaryList[0])
    #print(len(binaryList))
#This block of code will compute the hamming distance via the matrices gotten from above
    matrixCharBySeq = np.asarray(binaryList)
    x = matrixCharBySeq
    #print(x.shape)
    a= np.reshape(x, (20081, 1642*22))
    #print(type(x[0, :]))
    # print(x[0,:])
    print(np.shape(a))
    print(len(a[0]))
#test that a[0] is == to binarylist[0]
    # # print(matrixCharBySeq[0,0,:])
    # matrixSeqByChar = np.transpose(matrixCharBySeq, axes=(1642, 20081, 22))
    # matrixDot = np.dot(matrixCharBySeq, matrixSeqByChar)
    # matrixHammings = len(listSeq[0]) - matrixDot
    # print(matrixHammings)
    #
    # newF = open("testFileWrite.txt", 'w')
    # newF.write(str(matrixHammings))

startTime = timeit.default_timer()
main()
print("This algorithm takes: " + str(timeit.default_timer() - startTime) + " seconds")


def test():
    newF = open("testFileWrite.txt", 'w')
    a = np.array([[2, 3], [4, 5]])
    b = np.transpose(a)
    c = np.dot(a, b)
    newF.write(str(c))

#think about options to save the data in a txt file or something else
#next step is analysis of data and clustering, after optimization of the algorithm