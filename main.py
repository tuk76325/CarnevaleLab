import Bio.AlignIO
import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder

def hamming_distance(string1, string2):
    distance = 0
    L = len(string1)
    for i in range(L):
        if string1[i] != string2[i]:
            distance += 1
    return distance

def main():
    file = open(r'C:\Users\Toshan\PycharmProjects\BioinformaticsDoodnauth\PF00520_rp15.txt') # original fasta file
    f = open("SeqPF00520_rp15.txt", 'w')
    fileWrite = open("SeqPF00520_rp15.txt", 'r') #file with all sequences line by line
    hammingDis = open("hammingDistanceMatrix.csv", 'w')
    alignment = Bio.AlignIO.read(file, "fasta")
    listSeq = []

    c = 0
    for record in alignment:
        #print(record.seq)
        f.write(str(record.seq) + "\n")
        listSeq.append(str(record.seq).upper())
    #print(listSeq[0])

    aminoAcidList = ['A', 'G', 'I' , 'L', 'P', 'V', 'F', 'W', 'Y', 'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q', '-']
    df = pd.DataFrame(aminoAcidList)
    ohe = OneHotEncoder()
    arraySequences = ohe.fit_transform(df).toarray()
    #print(arraySequences)

    newDict = {}
    for i in range(len(aminoAcidList)):
        newDict.update({aminoAcidList[i] : arraySequences[i]})
    #print(newDict)

    newList = []
    for character in listSeq[0]:
        newList.append(newDict[character])
    print(type(newList[0]))

    newListB = []
    for character in listSeq[1]:
        newListB.append(newDict[character])
    print(type(newListB[0]))

    m = 0
    for i in range((len(aminoAcidList) * len(listSeq[0])) - 1):
        m += newList[i]*newListB[i]
    print(m)

'''for line in listSeq:
        print(line)
        for line2 in listSeq:
            hammingDis.write(str(hamming_distance(line, line2)) + ',')
        hammingDis.write('\n #')'''

main()


#USE THE INPUT FXN TO INSERT AN INTEGER THAT U GET FROM THE MATRIX TO GET SEQ ID FROM THE COLUMN
#THEN USE THE ROWS TO GET THE THE NUMBER AT WHICH THE BASE IS AT