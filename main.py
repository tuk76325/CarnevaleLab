import Bio.AlignIO
import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder

def main():
#This block of code opens the text file so that python can read and parse the protein sequences
    file = open(r'C:\Users\Toshan\PycharmProjects\BioinformaticsDoodnauth\PF00520_rp15.txt') # original fasta file
    f = open("SeqPF00520_rp15.txt", 'w')
    fileWrite = open("SeqPF00520_rp15.txt", 'r') #file with all sequences line by line
    hammingDis = open("hammingDistanceMatrix.csv", 'w')
    alignment = Bio.AlignIO.read(file, "fasta")
    listSeq = []

#This creates a list that has a sequence at each index, by this I mean that sequence 1 is the first position in the list, sequence 2 is the next position, and so on.
    for record in alignment:
        #print(record.seq)
        f.write(str(record.seq) + "\n")
        listSeq.append(str(record.seq).upper())
    #print(listSeq[0])

    aminoAcidList = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W', 'Y', 'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q', '-']
#This block of code creates the binary arrays (vectors) for each character in the aminoAcidList
    df = pd.DataFrame(aminoAcidList)
    ohe = OneHotEncoder()
    arraySequences = ohe.fit_transform(df).toarray()

#This dictionary attaches the amino acid character from aminoAcidList to the corresponding binary vector
    newDict = {}
    for i in range(len(aminoAcidList)):
        newDict.update({aminoAcidList[i] : arraySequences[i]})
    #print(newDict)

#This block of code currently appends each protein sequence to a list called newList and newListB so that it can be one-hot encoded later
    newList = []
    for character in listSeq[0]:
       newList.append(newDict[character])
    #print(type(newList[0]))

    newListB = []
    for character in listSeq[1]:
        newListB.append(newDict[character])
    #print(type(newListB[0]))

#This block of code computes the Hamming's Distance via he dot product of the two binary arrays (aka vectors) and then subtracting the length
    m = 0
    for i in range(len(newList)):
        m += np.dot(newList[i], newListB[i])
    print(len(newList) - m) #This is the Hamming's Distance


main()


#USE THE INPUT FXN TO INSERT AN INTEGER THAT U GET FROM THE MATRIX TO GET SEQ ID FROM THE COLUMN
#THEN USE THE ROWS TO GET THE THE NUMBER AT WHICH THE BASE IS AT