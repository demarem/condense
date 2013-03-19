'''
Reads a nexus formatted tree file and combines identical sequences.

Created: February 22, 2013
Author: Matthew Demarest
Version: 1.0
'''

import re
import argparse
import string
from itertools import combinations
from collections import namedtuple

def parseArgs():
    parser = argparse.ArgumentParser(\
        description='A Python script for condensing Nexus formatted tree files.')
    parser.add_argument("inFile", metavar="inputFile", type=argparse.FileType('rU'),
                help="input file", nargs=1)
    parser.add_argument("outFile", metavar="outputFile", type=argparse.FileType('w'),
                help="output file", nargs=1)

    args = parser.parse_args()
    # debugging argument parsing
    # print args

    return args

def condense(inputFile, outputFile):
    header = saveHeader(inputFile)
    nameToSequence = findIdenticals(inputFile)
    nameToSequenceStr = buildMapping(nameToSequence)
    generateNewHeader(outputFile, header, nameToSequenceStr, len(nameToSequence))
    generateBody(outputFile, nameToSequence)
    generateFooter(inputFile, outputFile)

def generateFooter(inputFile, outputFile):
    for line in inputFile:
        outputFile.write(line)


def generateBody(outputFile, nameToSequence):
    for k, v in iter(sorted(nameToSequence.iteritems())):
        line = "%s\t%s\n" % (''.join(k), v.sequence)
        outputFile.write(line)
    outputFile.write(';\n')

def generateNewHeader(outputFile, header, nameToSequenceMap, numSequences):
    lines = header.split('\n')

    for line in lines:
        if "#nexus" in line.strip().lower():
            line += "\n\n" + nameToSequenceMap + "\n"
        elif "dimensions" and "ntax" in line.strip().lower():
            line = re.sub(r'(?i)NTAX=(\d*)', 'NTAX=' + str(numSequences), line)

        outputFile.write(line + '\n')



def buildMapping(nameToSequence):
    mapping = "[ MAPPING:\n"
    for k, v in iter(sorted(nameToSequence.items())):
        speciesList = ""
        for sp in v.species:
            speciesList += sp + " "
        mapping += "%s -> %s\n" % (''.join(k), speciesList)

    mapping += "]\n"

    return mapping


def findIdenticals(inputFile):
    idGenerator = nextID()
    Sequence = namedtuple('Sequence', 'species sequence')
    nameToSequence = {}
    for line in inputFile:
        if line.strip().lower() in { ';', 'end;' }:
            break

        parts = re.split(r'\s+', line)
        spec = parts[0]
        seq = parts[1]

        if spec and seq:
            wasInserted = False
            for k, v in nameToSequence.items():
                if v.sequence == seq:
                    nameToSequence[k].species.append(spec)
                    wasInserted = True

            if wasInserted == False:
                name = idGenerator.next()
                nameToSequence[name] = Sequence([spec], seq)

    return nameToSequence

def saveHeader(inputFile):
    header = ""
    for line in inputFile:
        header += line
        if line.strip().lower() == "matrix":
            break

    return header


def nextID():
    width = 1
    while True:
        for com in combinations(string.ascii_uppercase, width):
            yield com
        width += 1

def main():
    args = parseArgs()
    condense(args.inFile[0], args.outFile[0])



if __name__ == '__main__':
    main()
