#!/usr/bin/env python
# coding: utf-8

# **CLASS COMMANDLINE**
#
# This class specifies a commandline to use for the program.
#
# Functions:
#
#     __init__(self, inOpts=None) : class constructor




########################################################################
# Class: CommandLine
#
# Purpose:class to represent a command line
#
#
# Author: Chris Nguyen
# History:      cdn 10/18/21 Created
#
########################################################################
import argparse


class CommandLine():

    def __init__(self, inOpts=None):
        '''
        Class constructor that creates the command line.

        '''

        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        # these arguments can be accessed or seen using "args.(argument name)" so args.list would hold the list of values passed to the option
        self.parser.add_argument('-i', '--minMotif', type=int, choices=range(3, 9), action='store', required=True,
                                 help='minimum kmer length')
        self.parser.add_argument('-a', '--maxMotif', type=int, choices=range(3, 9), action='store', required=True,
                                 help='maximum kmer length')
        self.parser.add_argument('-z', '--cutoff', type=int, action='store', required=True, help='zcutoff')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


# **CLASS RANGEERROR**
#
# Class to throw error
#
# Functions:
#
#     __init__(self, msg) : class constructor




class rangeError(Exception):
    def __init__(self, msg):
        '''
        Class constructor for rangeError
        '''
        self.msg = msg


# **CLASS MISSINGFINDER**
#
# This class finds the under-represented kmers in a data set.
#
# Functions:
#
#     __init__(self,infile, kmer, lowEnd, cutoff): class constructor
#     tupleFlip(self,tup): flips a tuple
#     nullMaker(self): makes the null model as a dictionary
#     complementer(self, aString): makes complement of a string
#     expCounter(self, kmer): finds expected counts for a kmer
#     expCountAndZScore(self): calculates z score and updates
#     fullPrint(self): prints out the desired output
#     danSearch(self): searches for each kmer by sliding a window




class missingFinder():

    def __init__(self, infile, kmer, lowEnd, cutoff):
        '''
        Constructor for the missingFinder class.
        Takes in parameters and stores them in member variables.

        Parameters:
            infile (file) - file to read from
            kmer (int) - max length of kmer to look for
            lowEnd (int) - lower end of kmer
            cutoff (int) - z score cutoff
        '''
        self.k = kmer  # high end of kmer range
        self.file = infile  # file to use
        self.low = lowEnd  # low end of kmer
        self.n = 0  # genome size
        self.cutoff = cutoff  # zcutoff
        # makes 3 dictionaries to store values
        self.countDict = self.nullMaker()
        self.exCountDict = self.nullMaker()
        self.zDict = self.nullMaker()

    def tupleFlip(self, tup):
        '''
        Function that flips a tuple around.

        Parameters:
            tup (tuple) - tuple to flip
        '''
        repTup = (tup[1], tup[0])
        return repTup

    def nullMaker(self):
        '''
        Function that makes dictionary of all possible kmers of a specified size

        '''
        import itertools
        import copy

        d1 = {}
        d2 = {}

        # for all values up to a specified size self.k, loop creation
        for j in range(1, self.k + 1):

            # itertools.product makes combinations of the letters
            topCount = itertools.product(['A', 'T', 'C', 'G'], repeat=j)
            kmerGen = (''.join(w) for w in topCount)

            # loops through generator to get values and complements
            for x in kmerGen:
                # complementer returns the complement
                y = self.complementer(x)

                # gets reversed and both parts get tuple'd and put into another tuple
                z = ''.join(reversed(y))
                tup = (x, z)
                key = (tuple(sorted(tup)))

                # puts new keys into dictionary, accounting for repeats
                if key in d2:

                    continue
                else:
                    d2[key] = 0

            # makes a copy and puts it into the larger dictionary, which is returned
            deepd2 = copy.deepcopy(d2)
            d1[f'k = {j}'] = deepd2

            d2.clear()
        return d1

    def complementer(self, aString):
        '''
        Function that makes complement of a string passed in

        Parameters:
            aString (str) - string to complement
        '''
        bString = ''
        for i in (aString):
            if i == 'A':
                bString += 'T'
            elif i == 'T':
                bString += 'A'
            elif i == 'C':
                bString += 'G'
            elif i == 'G':
                bString += 'C'

        return bString

    def expCounter(self, kmer):
        '''
        Function that calculates the expected count for some kmer

        Parameters:
            kmer (str) - kmer to calculate
        '''

        exCount = 0
        # takes the denominator first
        part3F = kmer[1: len(kmer) - 1]  # middle slice
        part3R = self.complementer(part3F)[::-1]  # reverse of slice
        part3 = (part3F, part3R)  # makes tuple to look for
        part3C = 0  # variable to hold the counts of the kmer found

        # looks for it
        for x in self.countDict:
            val = self.countDict[x].get(part3, self.tupleFlip(part3))  # if the key isn't found, the flip is returned
            if type(val) == tuple:  # checks if the flip was returned
                val2 = self.countDict[x].get(val)  # now, if the key really can't be found, 0 is returned
                if val2 == None:  # if val is None, look in a different dictionary
                    continue
                else:
                    part3C = val2  # but if it's found, then it's found. stop looking
                    break
            else:  # if val isn't a tuple, then the key was found successfully
                part3C = val  # so, just set it equal to the value returned and stop looking
                break
        if part3C == 0:
            return part3C  # if any part is 0, then the whole thing is 0, so don't bother with the rest

        # now, do part 1

        part1F = kmer[0: len(kmer) - 1]  # starting slice
        part1R = self.complementer(part1F)[::-1]  # reverse of slice
        part1 = (part1F, part1R)  # makes tuple to look for
        part1C = 0  # variable to hold the counts of the kmer found

        # looks for it
        for x in self.countDict:
            val = self.countDict[x].get(part1, self.tupleFlip(part1))  # if the key isn't found, the flip is returned
            if type(val) == tuple:  # checks if the flip was returned
                val2 = self.countDict[x].get(val)  # now, if the key really can't be found, 0 is returned
                if val2 == None:  # if val is None, look in a different dictionary
                    continue
                else:
                    part1C = val2  # but if it's found, then it's found. stop looking
                    break
            else:  # if val isn't a tuple, then the key was found successfully
                part1C = val  # so, just set it equal to the value returned and stop looking
                break
        if part1C == 0:
            return part1C

            # now, part 2
        part2F = kmer[1: len(kmer)]  # starting slice
        part2R = self.complementer(part2F)[::-1]  # reverse of slice
        part2 = (part2F, part2R)  # makes tuple to look for
        part2C = 0  # variable to hold the counts of the kmer found

        # looks for it
        for x in self.countDict:
            val = self.countDict[x].get(part2, self.tupleFlip(part2))  # if the key isn't found, the flip is returned
            if type(val) == tuple:  # checks if the flip was returned
                val2 = self.countDict[x].get(val)  # now, if the key really can't be found, 0 is returned
                if val2 == None:  # if val is None, look in a different dictionary
                    continue
                else:
                    part2C = val2  # but if it's found, then it's found. stop looking
                    break
            else:  # if val isn't a tuple, then the key was found successfully
                part2C = val  # so, just set it equal to the value returned and stop looking
                break
        if part2C == 0:
            return part2C

            # if the code is here, that means it's passed all the 0 checks. now, calculate E(K) and return that.

        exCount = part1C * part2C / part3C

        return exCount

    def expCountAndZScore(self):
        '''
        Function that updates expected count and z score dictionaries.
        '''
        import math
        # gSize = self.countDict['k = 1'][('A', 'T')] + self.countDict['k = 1'][('C', 'G')] #gSize is size of genome
        # self.n = gSize #sets self.n to genome size
        for outer in self.countDict:  # outer refers to the outer dict keys
            for inner in self.countDict[outer]:  # now we're considering the tuplekey parts
                kmer = inner[0]  # takes just one of the strings

                exCount = self.expCounter(kmer)  # returns expected count
                if exCount == 0:  # if it returns 0, some count was off, and z should be 0.
                    self.exCountDict[outer][inner] = 0
                    continue  # try another tuplekey
                else:  # if it didn't return 0, update the count and do the rest of the calculation
                    self.exCountDict[outer][inner] = exCount  # updates expected count
                    observed = self.countDict[outer][inner]  # takes observed count
                    stanDev = math.sqrt(((1 - (exCount / self.n)) * exCount))  # figures standard dev
                    z = (observed - exCount) / stanDev  # calc z score
                    if z >= self.cutoff:  # if value is greater than the cutoff, destroy it
                        self.zDict[outer].pop(inner)
                    else:
                        self.zDict[outer][inner] = z  # sets z score

    def fullPrint(self):
        '''
        Function that prints output, formatted
        '''
        print(self.file)
        print('N = ' + str(self.n))
        print('{0:16}\t{1:<8s}\t{2:^8s}\t{3:>8s}'.format('sequence:reverse', 'count', 'expected', 'z-score'))
        for x in range(self.k, self.low - 1, -1):
            key = f'k = {x}'
            sortZ = dict(sorted(self.zDict[key].items(), key=lambda item: item[1]))
            for zKey in sortZ:
                # print(str(zKey) + '\t' + str(self.countDict[key][zKey]) + '\t' + str(self.exCountDict[key][zKey]) + '\t' + str(sortZ[zKey]))
                print('{0:8}:{1:8}\t{2:<8d}\t{3:^8.1f}\t{4:>0.2f}'.format(zKey[0], zKey[1], self.countDict[key][zKey],
                                                                          self.exCountDict[key][zKey], sortZ[zKey]))

    def danSearch(self):
        '''
        Function that does the brunt of the searching
        '''

        with open(self.file) as file:
            check = file.read(1)  # reads one char in to see if it's the header
            if check == '>':  # if it's the header, readline() to skip it
                file.readline()
            keySet = ()  # this is the empty set of keys to compare later

            chunk = file.read(self.k - 1)  # reads the first k-length kmer
            chunk = check + chunk  # puts the first base back on
            self.n += len(chunk)  # adds chunk length to self.n
            # do one pass of danSearch
            for x in range(1, (self.k + 1)):  # use 1 - k+1 since x is the rightside bound of the slice
                kmer = chunk[0:x]  # kmer holds slices of the first chunk
                rever = self.complementer(kmer)[::-1]  # flips the complement around
                keyPair = (kmer, rever)  # package it up in a tuple
                keySet += (keyPair,)  # adds the tuple to keySet

            # search the dictionaries for the key
            for pair in keySet:
                for x in self.countDict:
                    if pair in self.countDict[x]:  # since dict1[x] is STILL a dictionary, checks if the pair is a key
                        self.countDict[x][pair] += 1  # if so, increment

            # dansearch for the first chunk is done, so now move to looping until end

            while True:
                keySet = ()  # makes keySet point at empty tuple
                newBase = file.read(1)  # reads in next base
                if newBase == '\n':  # if the character is a newline, keep going
                    continue
                chunk = chunk[1:]  # if everything checks out, chops the head off
                chunk += newBase  # concatenates onto the chunk
                self.n += 1

                # and then, do the whole search process again
                for x in range(1, len(chunk) + 1):
                    kmer = chunk[0:x]
                    rever = self.complementer(kmer)[::-1]
                    keyPair = (kmer, rever)
                    keySet += (keyPair,)

                for pair in keySet:
                    sPair = tuple(sorted(pair))
                    for x in self.countDict:
                        if sPair in self.countDict[x]:
                            self.countDict[x][sPair] += 1
                if chunk == '>':  # if there are other headers, do the whole thing again
                    self.danSearch()
                elif chunk == '':  # otherwise, stop
                    break


# **MAIN**
#
# This contains main, which makes object instances and runs the search.




########################################################################
# File: Main
#
# Purpose:main function for randomizedMotifSearch program
#
#
# Author: Chris Nguyen
# History:      cdn 10/18/21 Created
#
########################################################################
import random
import math


def main(inFile, myCommandLine=None):
    import sys

    try:
        cLine = CommandLine(myCommandLine)
        mini = cLine.args.minMotif
        maxi = cLine.args.maxMotif
        if mini >= maxi:
            raise rangeError('Minimum cannot be greater than maximum, please re-try')
        zcutoff = cLine.args.cutoff
        missing1 = missingFinder(inFile, maxi, mini, zcutoff)
        missing1.danSearch()
        missing1.expCountAndZScore()
        missing1.fullPrint()
    except rangeError as err:
        print(err.msg)


if __name__ == "__main__":
    main('Arthrospira-platensis-NIES-39.fna', ['-i 3', '-a 8', '-z 0'])


