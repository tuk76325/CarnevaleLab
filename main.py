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

#This block of code currently appends each protein sequences' binary vector to a list called newList and newListB
    binaryList = []
    x=0
    while x < len(listSeq):
        for character in listSeq[x]:
            binaryList.append((newDict[character]))
        #print((binaryList))
        x+=1

    binaryListB = binaryList.copy()


#This block of code computes the Hamming's Distance via he dot product of the two binary arrays (aka vectors) and then subtracting the length
    m = 0
    for i in range(len(newList)):
        m += np.dot(binaryList[i], binaryListB[i])
    print(len(newList) - m) #This is the Hamming's Distance

#The problem now is implementing this for loop to work for every sequence in the fasta file and having it do increment automatically @ O(n)
#After completing this, the next step is to introduce each sequence into a matrix such that the hammings distance is shown for all seqs

main()


#USE THE INPUT FXN TO INSERT AN INTEGER THAT U GET FROM THE MATRIX TO GET SEQ ID FROM THE COLUMN
#THEN USE THE ROWS TO GET THE THE NUMBER AT WHICH THE BASE IS AT
