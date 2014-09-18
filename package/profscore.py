#-------------------------------------------------------------------------------
# Name:        profscore
#
# Purpose:     Calculates Valdar and profile-based conservation score, with
#              p-values estimated by random sampling from a larger set of
#              sequences.
#
# Author:      Tet Woo Lee <tw.lee@auckland.ac.nz>
#
# Created:     29/03/2011
#
# Copyright:   (c) Tet Woo Lee 2011-2014
#
# Licence:     GNU General Public License version 3
#
#              This program is free software: you can redistribute it and/or modify
#              it under the terms of the GNU General Public License as published by
#              the Free Software Foundation, either version 3 of the License, or
#              (at your option) any later version.
#
#              This program is distributed in the hope that it will be useful,
#              but WITHOUT ANY WARRANTY; without even the implied warranty of
#              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#              GNU General Public License for more details.
#
#              You should have received a copy of the GNU General Public License
#              along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------

# Tested with Python 2.7.6 (default, Nov 10 2013, 19:24:24) [MSC v.1500 64 bit (AMD64)]
import argparse
import binascii # used for calculating checksum only
from collections import OrderedDict
import csv
from functools import total_ordering
import itertools
import operator
import os
import random
import re
from sets import Set
import sys
import time     # used for getting seed for PRNG
import timeit
import traceback

# Tested with BioPython 1.63 amd64 py2.7
import Bio.SeqIO
import Bio.SubsMat
import Bio.SubsMat.MatrixInfo
import Bio.Alphabet

# Tested with numpy 1.9.0 amd64 py2.7
from numpy import array
from numpy import float64
import numpy



##########################################################################################
# GENERIC HELPER FUNCTIONS
##########################################################################################

def createRecordChecksum(records):
  """
  Generated a crc32 checksum of a list of sequence records.
  Used to uniquely identify a list of sequences.

  Arguments:

    records : list of Bio.SeqRecord objects
        list of records to process

  Return:

    out : number
        32-bit checksum of sequences

  """
  checksum = 0
  for record in records:
      crc32 = binascii.crc32(str(record.seq))
      checksum+=crc32
  return checksum & 0xffffffff # & forces unsigned

def checkTotalWeight(records,weights):
  """
  Return ``True`` if sum of weights adds up to ``1.0``.

  Arguments:

    records : list of Bio.SeqRecord objects
      list of aligned sequences to test weights for
    weights : NamedDict of sequence weights
      sequences weights to check
  """
  totalWeight = 0.0
  for record in records:
      weight = weights[record]
      if verbose>=4: print record.id, 'has weight',weight
      totalWeight+=weight

  if verbose>=3: print 'Check sequence weights: Total weight is',totalWeight
  return abs(totalWeight-1.0)<1e-9

class NamedDict(dict):
  """
  Simple subclass of ``dict`` that contains a ``name`` attributed.
  """
  __slots__ = ('name') # use slots to conserve memory
  def __init__(self,name):
    dict.__init__(self)
    self.name = name #: *string*, name of this dict

def printRecordNames(records):
  """
  Simple function to print record ids.

  Arguments:

      records : list of Bio.SeqRecord objects
        list of records to print record ids
  """
  print " ".join(record.id for record in records)


def isProfileScore(scoreType):
  """
  Return ``True`` if ``scoretype`` is a profile-based score.
  """
  return scoreType == ScoreType.PROFILE or scoreType == ScoreType.PROFILE_UNCORR


def subsetByID(records, nameREstr):
    """
    Return a subset of sequences by regular expression search.

    Arguments:

      records : list of Bio.SeqRecord objects
        list of aligned sequences to process
      nameREstr : string
        regular expression to use to subset sequences, sequences must contain this regular expression


    Return:

      out : list of Bio.SeqRecord objects
        list of sequences containing the regex
    """
    print "Creating a subset of records by searching for regex %s in id..."%nameREstr
    nameRE = re.compile(nameREstr)
    sample = []
    for record in records:
        result = nameRE.search(record.id)
        if verbose>=2: print "Sequence %s %s."%(record.id,"contains regex" if result is not None else "does not contain regex")
        if result is not None:
            sample.append(record)
    return sample

def enum(*sequential, **named):
  """
  Simple method to create enum-like objects.

  Usage: ``enum('APPLE','ORANGE','PEAR')``

  Source: Alec Thomas, stack overflow http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python

  """
  enums = dict(zip(sequential, range(len(sequential))), **named)
  reverse = dict((value, key) for key, value in enums.iteritems())
  enums['reverse_mapping'] = reverse
  return type('Enum', (), enums)


class Dummy:
  """
  Dummy object class for storing data.
  """
  pass

##########################################################################################
# ENUMS
##########################################################################################

ScoreType = enum('VALDAR','PROFILE','PROFILE_UNCORR')
"""
Simple enum-like class for specifying score type.

Possible values:

  ``VALDAR``
    Valdar score
  ``PROFILE``
    Profile-based score, corrected for self-comparisons
  ``PROFILE_UNCORR``
    Profile-based score, not corrected for self-comparisons

"""

FilePosition = enum('START','HEADER','DATA','END')
"""
Simple enum-like class for specifying file position when reading randomization file.

Possible values:

  ``START``
    Initial state, before reading magic start word
  ``HEADER``
    Reading header with parameterrs
  ``DATA``
    Reading data site by site
  ``END``
    After reading magic end word

"""

OutputPositions = enum('ALL_SIMPLENUM','ALL','PRESENT','REFERENCE')
"""
Simple enum-like class for specifying output positions when reading reports.

Possible values:

  ``ALL``
    output all positions in alignment
  ``PRESENT``
    output only present positions in subset
  ``REFERENCE``
    output only positions present in reference sequence

"""


##########################################################################################
# SUBSTITUTION MATRIX FUNCTIONS
##########################################################################################

def readSubsMat(filename):
    """
    Read a substitution matrix from file in aaindex2 format.

    Arguments:

      filename : string
          path to the file to read, file should be in aaindex2 format

    Return:

      out : Bio.SubsMat.SeqMat object
          substitution matrix read from file

    """
    infile = open(filename,'r')
    print 'Reading substitution matrix:',filename,'...'
    row = -1 # -1 = reading header, >=0 reading matrix
    desc  = '' # description from 'D' line
    for line in infile:
        if row<0:
            if line[0] == 'D': # D = description
                desc = line[1:].strip()
                print 'Description:', desc
            if line[0] == 'M': # M = matrix info (i.e. letters)
                line = line.replace('=',' ').replace(',',' ')
                Mtokens = line.split()
                rowstoken = Mtokens.index('rows')+1
                colstoken = Mtokens.index('cols')+1
                rows = Mtokens[rowstoken]
                cols = Mtokens[colstoken]
                if verbose: print 'Rows:', rows
                if verbose: print 'Cols:', rows
                if rows!=cols:
                    raise Exception('rows and cols must be identical')
                seqmat = Bio.SubsMat.SeqMat(mat_name = desc)
                seqmat.alphabet.letters = rows
                if verbose: print 'Reading matrix...'
                row = 0
        else:
            entries = line.split()
            if entries[0] == '//':
                if verbose: print 'Reached end of matrix.'
                seqmat._correct_matrix()
                if verbose:
                  seqmat.print_mat( format="%5.1f", bottomformat="%5s", alphabet=seqmat.alphabet.letters)
                break
            row_letter = rows[row]
            col = 0
            for entry in entries:
                col_letter = cols[col]
                seqmat[col_letter,row_letter] = float64(entry)
                col+=1
            row+=1

    infile.close()
    return seqmat

def linearTransform(subsMat, diagonal = None, gap = 0.0, gapgap = None):
    """
    Linear transform a substitution matrix, so that min is 0.0 and max is 1.0.
    Note that the matrix is transformed in place.

    Arguments:

      subsMat : Bio.SubsMat.SeqMat object
        the substitution matrix to transform
      diagonal : number, optional
        value given to diagonal elements in matrix, if ``None`` (default) diagonal to set to "rounded average of diagonal elements in the unmodified matrix" (Valdar 2002),
        applied before transformation
      gap : number, optional
        score given to aa -> gap substitutions, defaults to ``0.0``, applied after transformation
      gapgap : number, optional
        score for gap -> gap substitutions, if ``None``, same as gap score, applied after transformation

    Return:

      out : Bio.SubsMat.SeqMat object
        Transformed substitution matrix.
    """
    print 'Linear transforming matrix:',subsMat.mat_name

    if gapgap is None: gapgap = gap

    # transform diagonal to "rounded average of diagonal elements in the unmodified matrix" (Valdar 2002)
    aaletters = 'ARNDCQEGHILKMFPSTWYV' # use only aa letters
    if diagonal is None:
        diagonalSum = 0.0
        n = 0
        for a in aaletters:
            m = subsMat[a,a]
            diagonalSum+=m
            n+=1
        diagonalScore = round(diagonalSum/n,0)
        subsMat.mat_name+=', average diagonal, linear transform, gap %.3f, gap-gap %.3f'%(gap,gapgap)
    else:
        diagonalScore = diagonal
        subsMat.mat_name+=', diagonal %.3f, linear transform, gap %.3f, gap-gap %.3f'%(diagonal,gap,gapgap)

    if verbose>=2: print 'Diagonal score set to:',diagonalScore
    for a in subsMat.alphabet.letters:
        subsMat[a,a] = diagonalScore
    min = subsMat[aaletters[0],aaletters[0]]
    max = min
    for a in aaletters:
        for b in aaletters:
            if a<b:
                m = subsMat[a,b]
                if m>max: max = m
                if m<min: min = m
    if verbose>=2: print 'Matrix min:',min
    if verbose>=2: print 'Matrix max:',max
    subsMat.alphabet.letters+='-'
    for a in subsMat.alphabet.letters:
        for b in subsMat.alphabet.letters:
            if a<=b:
                if a=='-' and b=='-':
                    Mut = gapgap
                elif a=='-':
                    Mut = gap
                else:
                    m = subsMat[a,b]
                    Mut = (m - min) / (max - min)
                subsMat[a,b]  = Mut
                subsMat[b,a]  = Mut
    print "New matrix:",subsMat.mat_name
    if verbose: subsMat.print_mat( format="%6.3f", bottomformat="%6s", alphabet=subsMat.alphabet.letters)


class VectorSubsMat():
    """
    Generate a substitution matrix in vector form from a (standard) substitution matrix.
    If ``vectorSubsMat`` is an instance of this class, then each lookup of ``vectorSubsMat[aa1]``
    (where ``aa1`` is amino acid in upper case one letter code)
    provides a *k*-dimensional vector giving scores for all possible substitutions.

    If ``vector = vectorSubsMat[aa1]``, a particular substitution score get be found by
    ``vector[index]`` where ``index`` is  the index of a (second) aa.
    The index of an aa given by ``vectorSubsMat.alphabetLUT[ord(aa2)]``.

    Arguments:

      subsMat : Bio.SubsMat.SeqMat object
        the substitution matrix to convert

    """
    def __init__(self, subsMat):
        print 'Creating vector substitution matrix from:',subsMat.mat_name
        self.mat_name = 'vectors of '+subsMat.mat_name #: *string*, name of this substitution matrix
        self.alphabet = Bio.Alphabet.Alphabet() #: *Bio.Alphabet.Alphabet object*, alphabet of this substitution matrix
        self.alphabet.letters = subsMat.alphabet.letters #: *string*, letters of this alphabet
        self.subsMat = subsMat #: *Bio.SubsMat.SeqMat object*, source substitution matrix

        self.scoreLUT = {} #: *dict*, look-up table mapping amino acid letters to vectors, access by indexing ``VectorSubsMat`` object, i.e. ``vectorSubsMat[aa] = scoreVector``

        for a in subsMat.alphabet.letters:
            scores = numpy.zeros(len(subsMat.alphabet.letters),dtype=float64)
            for bi in range(len(subsMat.alphabet.letters)):
                b = subsMat.alphabet.letters[bi]
                if a<=b: # use lower triangle only
                    m = subsMat[a,b]
                else:
                    m = subsMat[b,a]
                scores[bi] =  m
            #VectorSubsMat[a] = scores
            self.scoreLUT[a] = scores


        # generate lookup table to get alphabet_index from ordinal
        maxOrdinal = 0
        for letter in self.alphabet.letters:
            ordinal = ord(letter)
            if ordinal>maxOrdinal: maxOrdinal = ordinal

        self.alphabetLUT = list(-1 for i in range(maxOrdinal+1))
        """
        *list*, look-up table of amino acid letters to substitution vector indexes.

        Usage:
          ``alphabetLUT[ord(aa)]`` where ``aa`` is the amino acid
          in upper-case one letter code e.g. ``A``, ``Y``, ``M`` etc.
        """
        for letterIndex in range(len(self.alphabet.letters)):
            letter = self.alphabet.letters[letterIndex]
            ordinal = ord(letter)
            self.alphabetLUT[ordinal] = letterIndex

        self.dimensions = len(self.alphabet.letters) #: *number*, number of dimensions in vector (i.e. aa characters)


    def __getitem__(self, key):
        return self.scoreLUT[key]



##########################################################################################
# CALCULATION OF SEQUENCE DISTANCES AND WEIGHTS
##########################################################################################

def evolDistance(recordI, recordJ, subsMat):
    r"""
    Calculate the evolutionary distance between two aligned sequence records given a substitution matrix.
    As per Valdar, 2001 PMID: 11093265

    .. math::
            D(s_i,s_j) = 1 - \frac{\displaystyle\sum_{x \in aligned_{ij}}M\Big (s_{i}(x),s_{j}(x)\Big )}{\displaystyle n(aligned_{ij})}

    Arguments:

      recordI : Bio.SeqRecord object
        first sequence
      recordJ : Bio.SeqRecord object
        second sequence
      subsMat : Bio.SubsMat.SeqMat object
        substitution matrix, transformed

    Returns:

      out : number
        evolutionary distance
    """
    seqI = recordI.seq
    seqJ = recordJ.seq
    if len(seqI)!=len(seqJ):
        raise ValueError("Sequences must be same length %s %i vs %s %i"%(recordI.id,len(seqI),recordJ.id,len(seqJ)))
    sum = float64(0.0)
    n = 0
    for i in range(len(seqI)):
        sj = seqI[i]
        sk = seqJ[i]
        if sj=='-' and sk=='-': continue
        Mut = subsMat[sj,sk]
        #print sj,sk,Mut
        sum+=Mut
        n+=1
    result = 1.0-sum/n
    if verbose>=3: print "Calculated evolutionary distance between %s and %s: %f"%(recordI.id,recordJ.id,result)
    return result

def createUniformDistanceMatrix(records, uniformDistance = 1.0):
    """
    Create a simple distance matrix with a uniform value for pairs of aligned sequences.

    Arguments:

      records : list of Bio.SeqRecord objects
        list of aligned sequences (only size of list is used by this function)
      uniformDistance : number, optional
        the value to use as the uniform distance between each pair of sequences, defaults to ``1.0``

    Returns:

      out : NamedDict
        A dict of giving uniform distance for each pair of records *ij* i.e ``out[recordI, recordJ] = uniformDistance``.
        Value of ``out[recordI, recordI]`` is undefined.
    """
    if verbose>=1: print "Creating uniform evolutionary distance matrix..."
    uniformDistance = float64(uniformDistance)
    distanceMatrix = NamedDict("uniform distance of %f"%uniformDistance)
    n = 0
    for recordI in records:
        for recordJ in records:
            if recordJ>recordI:
                distanceMatrix[recordI, recordJ] = uniformDistance
                distanceMatrix[recordJ, recordI] = uniformDistance
                n+=1
    print 'Created uniform distance matrix for %i pairwise comparisons'%n
    return distanceMatrix

def calculateDistanceMatrix(records, subsMat):
    r"""
    Calculate a matrix of evolutionary distances between pairs of aligned sequences given a substitution matrix.
    Distance calculations as per Valdar, 2001 PMID: 11093265

    .. math::
            D(s_i,s_j) = 1 - \frac{\displaystyle\sum_{x \in aligned_{ij}}M\Big (s_{i}(x),s_{j}(x)\Big )}{\displaystyle n(aligned_{ij})}

    Arguments:

      records : list of Bio.SeqRecord objects
        list of aligned sequences to process
      subsMat : Bio.SubsMat.SeqMat object
        substitution matrix, transformed

    Returns:

      out : NamedDict
        A dict giving evolutionary distance for each pair of records ij i.e ``out[recordI, recordJ] = distance(recordI, recordJ)``.
        Matrix is symmetrical and ``out[recordI, recordI]`` is undefined.

    """
    if verbose>=1: print "Calculated evolutionary distance matrix..."
    startTime = timeit.default_timer()
    distanceMatrix = NamedDict("evolutionary distance")
    n = 0
    for recordI in records:
        for recordJ in records:
            if recordJ>recordI:
                dist = evolDistance(recordI,recordJ,subsMat)
                distanceMatrix[recordI, recordJ] = dist
                distanceMatrix[recordJ, recordI] = dist
                n+=1
    print 'Calculated distance matrix for %i pairwise comparisons'%n
    endTime = timeit.default_timer()
    if verbose>=1: print 'Calculation time was %.3f secs.'%(endTime-startTime)

    return distanceMatrix

def calculateWeights(records, distanceMatrix, silent = False):

    r"""
    Calculate weight of each sequence *|s_i|* as mean evolutionary distance between *|s_i|* and all other sequences in the alignment.
    As per Valdar, 2001 PMID: 11093265

    .. math::
          w_i = \frac{\displaystyle\sum_{j\neq i}^{N}D(s_i,s_j)}{\displaystyle N-1}
    .. |s_i| replace:: s\ :sub:`i`

    Arguments:

      records : list of Bio.SeqRecord objects
        list of aligned sequences to process
      distanceMatrix : NamedDict of evolutionary distances
        pre-calculated distance matrix giving pairwise distances between all records to process
      silent : boolean, optional
        if ``True`` do not produce any output (used when called by randomization functions), defaults to ``False``

    Returns:

      out : NamedDict
        A dict of giving weight for each sequence record, i.e. ``out[recordI] = weightI``
    """
    if not silent and verbose: print "Calculating sequence weights from evolutionary distance matrix..."
    startTime = timeit.default_timer()
    weights = NamedDict("weighted by %s"%distanceMatrix.name)
    N = float64(len(records)-1)
    totalWeight = float64(0.0)
    for recordI in records:
        sumDistance = float64(0.0)
        for recordJ in records:
            if recordI==recordJ: continue
            sumDistance+=distanceMatrix[recordI, recordJ]
        weight = sumDistance/N
        weights[recordI] = weight
        totalWeight += weight
    # normalise to total 1
    for recordI in records:
        weight = weights[recordI]
        weight = weight/totalWeight
        weights[recordI] = weight
        if not silent and verbose>=2: print "Weight of sequence %s: %f"%(recordI.id,weight)
    endTime = timeit.default_timer()
    if verbose>=1: print 'Calculation time was %.3f secs.'%(endTime-startTime)
    return weights



##########################################################################################
# CONSERVATION SCORING
##########################################################################################

class ConservationScores():
    """
    A data storage class to store calculated sequence conservations scores.

    Arguments:

      records : list of Bio.SeqRecord objects
        list of aligned sequences to process
      scoreName : string
        name of this scoring system, added to ``self.metadata``
      scoreAgainst : string
        details of scoring (e.g. scoring against what matrix/profile), added to ``self.metadata``
      weightBy: string
        details of weighting, added to ``self.metadata``

    """

    def __init__(self, records, scoreName, weightBy, scoreAgainst):
        if len(records)==0: raise ValueError("Empty sequence records list provided!")
        self.records = records #: *list of Bio.SeqRecord objects*, list of records provided to this score object
        self.M = len(records[0].seq) #: *number*, alignment length, i.e. total number of sites in alignment
        self.metadata = OrderedDict() #: *OrderedDict*, metadata for this profile

        self.scores = None #: *list of numbers*, list of scores indexed by 0-based site number in alignment, ``None`` until scores calculated

        self.metadata = OrderedDict() #: *OrderedDict*, metadata for this profile
        self.metadata['score name'] = scoreName
        self.metadata['score sequences'] = "%i sequences (%08X checksum)"%(len(self.records),createRecordChecksum(self.records))
        self.metadata['score sequence ids'] = " ".join(record.id for record in self.records)
        self.metadata['score number of sites'] = "%i sites"%self.M
        self.metadata['score weighting'] = weightBy
        self.metadata['score against'] = scoreAgainst

    def allocateScores(self):
        """
        Allocate an array of length ``self.M`` for storing scores, with all values initially set to ``nan``.
        """
        self.scores = numpy.zeros(self.M, dtype=float64)
        self.scores[:] = numpy.nan




##########################################################################################
# VALDAR CONSERVATION SCORING
##########################################################################################

class ValdarScore(ConservationScores):
  r"""
  Calculate the Valdar conservation score for all sites in an alignment, given a list of sequence
  records, a substitution matrix and a distance matrix or list of sequence weights.
  As per Valdar, 2001 PMID: 11093265

  Scores are automatically calculated when the class in constructed. Once constructed,
  this class will contain list of Valdar scores for each site in the alignment in ``self.scores``.

  Inherits from ``ConservationScores``.

  Valdar score:

  .. math::
     V_x = \frac{\displaystyle \sum_{i}^{N}\sum_{j>i}^{N}w_i w_j M\Big (s_{ix},s_{jx}\Big ) } {\displaystyle \sum_{i}^{N}\sum_{j>i}^{N}w_i w_j }

  Arguments:

    records : list of Bio.SeqRecord objects
      list of aligned sequences to process
    subsMat : Bio.SubsMat.SeqMat object
      substitution matrix, transformed
    distanceMatrix : NamedDict of evolutionary distances, optional, may be ``None``
      pre-calculated distance matrix giving pairwise distances between all records (if not provided, calculated automatically; not used if weights provided)
    weights : NamedDict of sequence weights, optional, may be ``None``
      pre-calculated weights of all sequences (if not provided, calculated automatically from evolutionary distances)

  """
  def __init__(self, records, subsMat, distanceMatrix = None, weights = None):
    if verbose: print "Calculating Valdar conservation scores..."

    if not distanceMatrix and not weights: raise ValueError("Either distanceMatrix or weights must be provided")
    self.subsMat = subsMat #: *Bio.SubsMat.SeqMat object*, substitution matrix used to calculated scores

    self.weights = weights #: *NamedDict of sequence weights*, weights of sequences used to calculated scores

    # if necessary calculate distance matrix
    if not distanceMatrix and not self.weights:
        if verbose: print "Distance matrix not provided, calculating..."
        distanceMatrix = calculateDistanceMatrix(self.records,self.subsMat)
    if not self.weights:
        if verbose: print "Sequence weights not provided, calculating..."
        self.weights = calculateWeights(self.records, distanceMatrix)

    ConservationScores.__init__(self,records,"Valdar score",self.weights.name,subsMat.mat_name)

    self.allocateScores()
    startTime = timeit.default_timer()
    for site in xrange(self.M):
        self.scores[site] = self.valdarBySite(site, self.records, self.subsMat, weights = self.weights)
    endTime = timeit.default_timer()
    if verbose>=1: print 'Finished calculating scores for all sites %.3f secs.'%(endTime-startTime)

  @staticmethod
  def valdarBySite(site, records, subsMat, weights, silent = False):
      """
      Calculate the Valdar conservation score for a given site in an alignment, given a list of sequence
      records, a substitution matrix, a site number and a distance matrix or list of sequence weights.
      Static function allows it to be called without creating a ``ValdarScore`` object.
      As per Valdar, 2001 PMID: 11093265

      Arguments:

        site : number
          site in alignment to calculate score for (note 0-based index)
        records : list of Bio.SeqRecord objects
          list of aligned sequences to process
        subsMat : Bio.SubsMat.SeqMat object
          substitution matrix, transformed
        weights : NamedDict of sequence weights
          pre-calculated weights of all sequences
        silent : boolean, optional
          if ``True`` do not produce any output (used when called by randomization functions), defaults to ``False``

      Returns:

        out : number
          Valdar score for the site.
      """
      if not weights: raise ValueError("Sequence weights must be provided")

      if len(records)==0: raise ValueError("Empty sequence records list provided!")


      normaliser = float64(0.0) # normaliser (denominator) = sum of weights
      score = float64(0.0)
      # do weighted sum of scores
      for recordI in records:
          for recordJ in records:
              if recordJ>recordI:
                  seqI = recordI.seq
                  seqJ = recordJ.seq
                  weightI = weights[recordI]
                  weightJ = weights[recordJ]
                  weightIweightJ = weightI * weightJ
                  normaliser += weightIweightJ

                  sj = seqI[site]
                  sk = seqJ[site]
                  Mut = subsMat[sj,sk]
                  score+=weightIweightJ*Mut

      # renormalise by sum of weights
      normscore = score/normaliser
      if not silent and verbose>=4: print "site %i: score %f = %f/%f"%(site+1,normscore,score,normaliser)
      return normscore


##########################################################################################
# SEQUENCE PROFILE GENERATION AND SCORING
##########################################################################################

class AAScore:
    """
    Simple data storage class to store and order AA scores for finding consensus.
    Sorting ``AAScore`` objects results in ordering by score.

    """
    def __init__(self, letter, score, pMax):
        self.letter = letter #: *character*, letter for this aa
        self.score = score #: *number*, absolute score for this aa
        self.pMax = pMax #: *number*, score for this aa as a proportion of the max score at the site
    def __cmp__(self, other):
        return cmp(self.score,other.score)
    def __str__(self):
        return '%s:%f'%(self.letter,self.score)

class SequenceProfile():
    r"""
    Generate a SequenceProfile from set of alignment records.
    The sequence profile is the weighted average vector of amino acid scores for that position.

    The profile-based scoring system score breaks Valdar score into two parts, sequence profile generation (``SequenceProfile``) and profile scoring (``ProfileScore``).

    Profile automatically calculated when the class in constructed. Once constructed,
    this class will contain list (for each site in the alignment) of profile vectors as ``self.profile``.

    Sequence profile:

    .. math::
      \bar{P}_x = \frac{\displaystyle \sum_{i}^{N}w_iP_{ix}}{\displaystyle \sum_{i}^{N}w_i}

    Arguments:

      records : list of Bio.SeqRecord objects
        list of aligned sequences to process
      vectorSubsMat : VectorSubsMat object
        substitution matrix, transformed and converted to vector form
      distanceMatrix : NamedDict of evolutionary distances, optional, may be ``None``
        pre-calculated distance matrix giving pairwise distances between all records (if not provided, calculated automatically; not used if profileWeights provided)
      profileWeights : NamedDict of sequence weights, optional, may be ``None``
        pre-calculated weights of all sequences used for weighting profile (if not provided, calculated automatically from evolutionary distances),
        weights for sequences used should sum to 1.0

    """
    def __init__(self, records, vectorSubsMat, distanceMatrix = None, profileWeights = None):

        if verbose: print 'Generating sequence profile...'

        if not distanceMatrix and not profileWeights: raise ValueError("Either distanceMatrix or profileWeights must be provided")
        if len(records)==0: raise ValueError("Empty sequence records list provided!")
        self.records = records #: *list of Bio.SeqRecord objects*, list of records used to generate profile
        self.vectorSubsMat = vectorSubsMat #: *VectorSubsMat object*, vector substitution matrix used to generate profile

        if profileWeights:
            self.profileWeights = profileWeights #: *NamedDict of sequence weights*, weights of sequences
        else:
            #calculate profileWeights, normalised to total of 1.0
            self.profileWeights = calculateWeights(self.records, distanceMatrix)

        if not checkTotalWeight(self.records,self.profileWeights):
            raise Exception('Total weight should be 1.0')


        self.M = len(self.records[0].seq)  #: *number*, alignment length, i.e. total number of sites in alignment

        self.consensus = None
        """
          *list (for each site) of list (for multiple consensus residues) of AAScore objects*, ``None`` prior to call of ``findConsensus()``;
          stores consensus residues, as ``AAScore`` objects for each site, i.e. ``self.consensus[site]`` gives list of ``AAScore`` objects as consensus residues for the site

        """
        self.consensusDetails = None #: *string*, ``None`` prior to call of ``findConsensus()``;  details of how consensus determined, i.e. string giving cutoff values

        # generate full profile
        self.profile = [] #: *list of numpy.array*, list of vectors for the profile, indexed by 0-based site number in alignment

        sequences = len(self.records)
        self.metadata = OrderedDict() #: *OrderedDict*, metadata for this profile
        self.metadata['profile name'] = "Sequence profile"
        self.metadata['profile sequences'] = "%i sequences (%08X checksum)"%(len(self.records),createRecordChecksum(self.records))
        self.metadata['profile sequence ids'] = " ".join(record.id for record in self.records)
        self.metadata['profile number of sites'] = "%i sites"%self.M
        self.metadata['profile weighting'] = self.profileWeights.name
        self.metadata['profile substitution matrix'] = self.vectorSubsMat.mat_name
        
        startTime = timeit.default_timer()
        for i in range(self.M):
            profile_x = self.generateProfileAtSite(i, self.records, self.vectorSubsMat, self.profileWeights)
            self.profile.append(profile_x)
        endTime = timeit.default_timer()
        if verbose>=1: print 'Finished calculating profile for all sites %.3f secs.'%(endTime-startTime)

    @staticmethod
    def generateProfileAtSite(site, records, vectorSubsMat, profileWeights, silent = False):
        """
        Generate the sequence profile at a site in the alignment.
        The sequence profile is the weighted average vector of amino acid scores for that position.
        Static function allows it to be called without creating a ``SequenceProfile`` object.

        Arguments:

          site : number
            position in the alignment to calculate profile (note: 0-based)
          records : list of Bio.SeqRecord objects
            list of aligned sequences to process
          vectorSubsMat : VectorSubsMat object
            substitution matrix, transformed and converted to vector form
          profileWeights : NamedDict of sequence weights
            pre-calculated weights of all sequences used for weighting profile (weights for sequences used assumed to sum to 1.0)
          silent : boolean, optional
            if ``True`` do not produce any output (used when called by randomization functions), defaults to ``False``

        Returns:

          out : numpy.array
            a vector giving the sequence profile at requested site
        """
        profile_x = numpy.zeros(vectorSubsMat.dimensions, dtype=float64)

        for recordI in records:
            si = recordI.seq[site]
            weight = profileWeights[recordI]
            profile_x = profile_x + weight*vectorSubsMat[si] #add vector of scores to profile vector
        if not silent and verbose>=4: print "site %i: profile %s"%(site+1, ', '.join(['%s %f'%(vectorSubsMat.alphabet.letters[k],value) for k,value in enumerate(profile_x)]))
        return profile_x

    def findConsensus(self, cutoffAbs, cutoffProportion = 1.0, excludeGap = True):
        """
        Find consensus residues for each position in profile.
        For a residue to be a consensus residue at a site, it must be represented at that site in the alignment
        and must have a score at least equal to the cutoffs.
        Result is stored in attribute ``self.consensus``, with criteria given in ``self.consensusDetails``.

        Arguments:

          cutoffAbs : number
            absolute cutoff, residues must have absolutes scores >= this to be considered a consensus residue
          cutoffProportion : number, optional
            proportional cutoff, residues must have scores >= this proportion of the max score to be considered a consensus residue
          excludeGap : boolean, optional
            if ``True`` (default) prevents gaps from being considered consensus residues

        """
        self.consensus = []
        self.consensusDetails = 'Consensus: residues with score >= %.3f or >= %.2f of max score'%(cutoffAbs, cutoffProportion)
        if verbose: print self.consensusDetails
        letters = self.vectorSubsMat.alphabet.letters
        for site in range(self.M):
            scoresAtSite = self.profile[site]
            consensusAtSite = []
            maxScore = max(scoresAtSite)
            if maxScore>=cutoffAbs:
                for letter in letters:
                    if excludeGap and letter=='-': continue
                    aaPresent = False
                    # filter residues by presence in alignment
                    for recordI in self.records:
                        si = recordI.seq[site]
                        if si==letter:
                            aaPresent = True
                            break
                    if not aaPresent: continue
                    letterIndex = self.vectorSubsMat.alphabetLUT[ord(letter)]
                    score = scoresAtSite[letterIndex]
                    pMax = score/maxScore
                    if score>=cutoffAbs or pMax>=cutoffProportion:
                        aaScore = AAScore(letter, score, pMax)
                        consensusAtSite.append(aaScore)
            consensusAtSite.sort(reverse=True)
            self.consensus.append(consensusAtSite)
            if verbose>=3: print "site %i: consensus %s"%(site+1,", ".join("%s %.3f (%.2f)"%(aaScore.letter,aaScore.score,aaScore.pMax) for aaScore in consensusAtSite))

    def __getitem__(self, key):
        return self.profile[key]


class ProfileScore(ConservationScores):
  r"""
  Calculate sequence profile-based conservation scores for a set of alignment records given a sequence profile.

  The profile-based scoring system score breaks Valdar score into two parts, sequence profile generation (``SequenceProfile``) and profile scoring (``ProfileScore``).
  The score will be identical to the Valdar score if corrected for self-comparisons. This code does not assume that the sequence profile and profile scoring are
  generated from the same set of sequences, although it has not been tested or validated when different sets of sequences are used for each step.

  Scores are automatically calculated when the class in constructed. Once constructed,
  this class will contain list of conservation scores for each site in the alignment in ``self.scores``.

  Inherits from ``ConservationScores``.

  Uncorrected score:

  .. math ::
    \widetilde{C_x} = \frac{\displaystyle \sum_{i}^{N}w_i\bar{P}_{x}(s_{ix})}{\displaystyle \sum_{i}^{N}w_i}

  Corrected score:

  .. math ::
      C_x = \frac{\displaystyle \sum_{i}^{N}w_i\bar{P}_{x}(s_{ix}) -   \displaystyle \sum_{i}^{N}w_i w_i M\Big (s_{ix},s_{ix}\Big ) } {\displaystyle \sum_{i}^{N}w_i - \displaystyle \sum_{i}^{N}w_i w_i }


  Arguments:

    records : list of Bio.SeqRecord objects
      list of aligned sequences to process
    sequenceProfile : SequenceProfile object
      calculated sequence profile
    distanceMatrix : NamedDict of evolutionary distances, optional, may be ``None``
      pre-calculated distance matrix giving pairwise distances between all records (if not provided, calculated automatically; not used if weights provided)
    scoringWeights : NamedDict of sequence weights, optional, may be ``None``
      pre-calculated weights of sequences to be used for weighting the score (if not provided, calculated automatically from evolutionary distances),
      weights for sequences used should sum to 1.0
    correctForSelfComparisons : boolean, optional
      defaults to ``True``, if ``True``, score is corrected for self-comparisons

  """
  def __init__(self, records, sequenceProfile, distanceMatrix = None, scoringWeights = None, correctForSelfComparisons = True):

    correctedStr = "corrected" if correctForSelfComparisons else "uncorrected"

    if verbose: print 'Generating %s profile-based conservation scores...'%(correctedStr)

    if not distanceMatrix and not scoringWeights: raise ValueError("Either distanceMatrix or weights must be provided")
    #if len(records)==0: raise ValueError("Empty sequence records list provided!")

    self.sequenceProfile = sequenceProfile #: *SequenceProfile object*, sequenceProfile used to calculate scores

    self.scoringWeights = None #: *NamedDict of sequence weights*, weights of sequences used to calculated scores
    if not scoringWeights:
        self.scoringWeights = calculateWeights(records, distanceMatrix)  #: *NamedDict of sequence weights*, weights of sequences
    else: self.scoringWeights = scoringWeights

    if not checkTotalWeight(records,self.scoringWeights):
        raise Exception('Total weight should be 1.0')

    ConservationScores.__init__(self,records,"Profile-based score (%s for self-comparisons)"%(correctedStr),self.scoringWeights.name,sequenceProfile.metadata['profile name'])

    self.allocateScores()
    vectorSubsMat = self.sequenceProfile.vectorSubsMat
    profileWeights = self.sequenceProfile.profileWeights
    
    startTime = timeit.default_timer()
    for site in xrange(self.M):
        profileAtSite = self.sequenceProfile[site]
        self.scores[site] = ProfileScore.scoreAgainstProfileBySite(site, self.records, profileAtSite, vectorSubsMat, self.scoringWeights, profileWeights, correctForSelfComparisons = correctForSelfComparisons)
    endTime = timeit.default_timer()
    if verbose>=1: print 'Finished calculating profile scores for all sites %.3f secs.'%(endTime-startTime)


  @staticmethod
  def scoreAgainstProfileBySite(site, records, profileAtSite, vectorSubsMat, scoringWeights, profileWeights = None, correctForSelfComparisons = True, silent = False):
    """
    Calculate conservation score at a site in the alignment.
    The score is based on comparison to ``sequenceProfile``.
    Static function allows it to be called without creating a ``ProfileScore`` object.

    Arguments:

        site : number
            position in the alignment to calculate profile (note: 0-based)
        records : list of Bio.SeqRecord objects
            list of aligned sequences to process
        profileAtSite : k-dimensional numpy array
            sequence profile at the site i.e weighted average vector of amino acid scores at site
        vectorSubsMat : VectorSubsMat object
            substitution matrix, transformed and converted to vector form; used for self-comparison correction
        scoringWeights : NamedDict of sequence weights
            pre-calculated weights of sequences to be used for weighting the score, weights for sequences used assumed to sum to 1.0
        profileWeights : NamedDict of sequence weights, optional
            pre-calculated weights of sequences to be used for weighting the profile; if absent assumed to be the same as sequence weights
        correctForSelfComparisons : boolean, optional
            defaults to ``True``, if ``True``, score is corrected for self-comparisons
        silent : boolean, optional
            if ``True`` do not produce any output (used when called by randomization functions), defaults to ``False``

    Returns:

      out : number
        conservation score for this site
    """

    mutAll = float64(0.0) # summed mut score for all sequences (including self-comparisons)
    weightAll = float64(1.0) # summed weights for all sequences (including self-comparisons), assumed to be 1.0
    mutDiagonal = float64(0.0) # summed mut score for self-comparisons
    weightDiagonal = float64(0.0) # summed weights for self-comparisons

    if profileWeights is None : profileWeights = scoringWeights

    if not correctForSelfComparisons: profileWeights = [] # cheat to skip self-comparison calculation block in loop

    # Two kinds of sequence weight used:
    # scoringWeights
    #  = Weight of each sequence used when generating scores, i.e. contribution of a sequence to conservation score given profile
    # profileWeight
    #  = Weight of each sequence used when generating profile, i.e. contribution of a sequence to the sequence profile
    # both stored so that the sum of weights for all sequences is 1.0
    #
    # As scoring scheme is designed for the same set of sequence set used for generating profile and subset, these two will be generally by identical.
    # However, coded differently to allow for possibly use of different profile vs scoring sequences

    for recordI in records:
      sij = recordI.seq[site] # aa at site i in seq j
      letterIndex = vectorSubsMat.alphabetLUT[ord(sij)] # index of this aa in vector

      weightI = scoringWeights[recordI] # scoring weight of seq j
      m_ij = profileAtSite[letterIndex] # score of seq j from profile vector

      if recordI in profileWeights:
        # seq j was in profile, calculate self-specific mut needed for self-comparison comparison
        weightP_I = profileWeights[recordI] # profile weight of seq j
        mutDiagonal += weightI * weightP_I * vectorSubsMat[sij][letterIndex]
        weightDiagonal += weightI * weightP_I  # sum all diagonal weights

      mutAll+=weightI * m_ij

    if correctForSelfComparisons:
      # corrected score = mutNonDiagonal/weightNonDiagonal = (mutAll - mutDiagonal) / (weightAll - weightNonDiagonal) = (mutAll - mutDiagonal) / (1.0 - weightNonDiagonal)
      score = (mutAll - mutDiagonal) / (weightAll - weightDiagonal)
      if not silent and verbose>=4: print "site %i: score %f, all %f/%f, diagonal %f/%f, non-diagonal %f/%f"%(site+1,score,mutAll,weightAll,mutDiagonal,weightDiagonal,mutAll - mutDiagonal,weightAll-weightDiagonal)
    else:
      # uncorrected score = mutAll/weightAll; weightAll is normalised to 1.0
      score = mutAll
      if not silent and verbose>=4: print "site %i: score %f, all %f/%f"%(site+1, score, mutAll,1.0)
    return score


##########################################################################################
# SEQUENCE SAMPLING
##########################################################################################

def randomSampleReplacement(records, number):
    """
    Return a sample of sequences by sampling with replacement.

    Arguments:

      records : list of Bio.SeqRecord objects
        list of aligned sequences to process
      number : number
        number of sequences to sample

    Return:

      out : list of Bio.SeqRecord objects
        random sample of sequences
    """
    sampleIndexes = []
    lastIndex = len(records)-1
    while len(sampleIndexes)<number:
        index = random.randint(0, lastIndex)
        sampleIndexes.append(index)
    return [records[sampleIndex] for sampleIndex in sampleIndexes]

def randomSampleNoReplacement(records, number):
    """
    Return a sample of sequences by sampling without replacement.

    Arguments:

      records : list of Bio.SeqRecord objects
        list of aligned sequences to process
      number : number
        number of sequences to sample

    Return:

      out : list of Bio.SeqRecord objects
        random sample of sequences
    """
    if number>len(records): raise ValueError("Cannot sample more sequences than number provided!")
    sampleIndexes = Set()
    lastIndex = len(records)-1
    while len(sampleIndexes)<number:
        index = random.randint(0, lastIndex)
        sampleIndexes.add(index)
    return [records[sampleIndex] for sampleIndex in sampleIndexes]


##########################################################################################
# SCORE PROBABILITY DISTRIBUTION
##########################################################################################

class AlignmentPDF:
    r"""
    A conservation score probability distribution for the alignment. Generated by repeatedly randomly
    sampling a certain number of sequences from a larger set of sequences (without replacement).

    When constructed, this object will store and validate arguments provided in constructor.
    The object can then be used to:

      * Open an saved PDF and check for compatibility and consistency with this object
        using ``readFromFile(checkFile = True)``.
      * Calculate a PDF by performing randomization iterations using ``randomize()``, possibly
        appending existing iterations from an input file.
      * Open an saved PDF and look up densities site-by-site using ``prepareLookup()`` and
        ``lookupDensities()``.

    Arguments:

        records : list of Bio.SeqRecord objects
            list of aligned sequences to process
        sampleSize : number
            sample size for random sampling operation
        scoreType : ScoreType.VALDAR, ScoreType.PROFILE or ScoreType.PROFILE_UNCORR
            type of score to calculate calculate
        subsMat : Bio.SubsMat.SeqMat or VectorSubsMat object
            substitution matrix, transformed and converted to vector form if profile-based score being calculated
        distanceMatrix : NamedDict of evolutionary distances, optional, may be ``None``
            pre-calculated distance matrix giving pairwise distances between all records (if not provided, calculated automatically)
    """

    startMagicWord = 'ProfileConsScore_PDFRand01_START' #: *string*, magic word at start of a saved ``ScoreProbDist``
    endMagicWord = 'ProfileConsScore_PDFRand01_END' #: *string*, magic word at end of a saved ``ScoreProbDist``
    dataMagicWord = 'note header ends' #: *string*, magic word for beginning of data section in file
    sequenceIdHeader = 'sequence ids' #: *string*, header used for sequence ids
    randomizationsHeader = 'total randomizations' #: *string*, header used for number of randomizations

    def __init__(self, records, sampleSize, scoreType, subsMat, distanceMatrix = None):

        self.records = records #: *list of Bio.SeqRecord objects*, list of records provided to this score object
        self.sampleSize = sampleSize #: *number*, sample size of random sample
        self.scoreType = scoreType #: *ScoreType*, type of score

        self.distanceMatrix = distanceMatrix #: *NamedDict of evolutionary distances*, matrix of distances between all sequences

        if isinstance(subsMat,VectorSubsMat):
            self.vectorSubsMat = subsMat #: *VectorSubsMat object*, vector substitution matrix to use for scoring (if profile-based scoring used)
            self.subsMat = self.vectorSubsMat.subsMat
        else:
            self.vectorSubsMat = None
            self.subsMat = subsMat #: *Bio.SubsMat.SeqMat object*, subsitution matrix to use for scoring (if Valdar scoring used)

        if self.distanceMatrix is None:
            self.distanceMatrix = calculateDistanceMatrix(self.records,self.subsMat)

        self.M = len(records[0].seq) #: *number*, alignment length, i.e. total number of sites in alignment
        self.totalRandomizations = 0 #: *number*, total number of randomizations in this PDF

        self.curSite = 0 #: *number*, current site being read or processed

        self.csvReader = None #: *csv.CSVReader* object, stored for use when reading randomization file
        self.filePosition = None #: one of *FilePosition* values, used when reading randomization file

    def prepareMetadata(self):
        """
        Generate and return metadata for this object.

        Return:

          out : OrderedDict
            OrderedDict containing metadata for this object.

        """
        metadata = OrderedDict()


        metadata['sequence number'] = "%i"%len(self.records)
        metadata['sequence checksum'] = "%08X"%createRecordChecksum(self.records)
        metadata[self.sequenceIdHeader] = " ".join(record.id for record in self.records)
        metadata['sequence number of sites'] =  "%i"%self.M

        metadata['sample size'] =  "%i"%self.sampleSize
        metadata['score type'] =  ScoreType.reverse_mapping[self.scoreType]
        if self.vectorSubsMat:
            metadata['substitution matrix'] = self.vectorSubsMat.mat_name
        else:
            metadata['substitution matrix'] = self.subsMat.mat_name
        metadata['distance matrix type'] = self.distanceMatrix.name
        metadata[self.randomizationsHeader] = self.totalRandomizations

        metadata[self.dataMagicWord] = "data follows immediately below"

        return metadata

    def randomize(self, outputHandle, numRandomizations, inputHandle = None):
        """
        Perform a certain number of randomizations to generate a PDF.

        In each randomization iteration, a random sample set of sequences is drawn
        for the list of records. Conservation scores are then calculated for this sample
        as per the method and options provided in the constructor.

        To conserve memory, randomizations are performed site-by-site. Scores for all
        randomizations for a given site are sorted in ascending order, binned (where
        identical scores are repeatedly obtained) and saved to disk.

        If an input file handle is given, randomizations will be
        read from that file and then appended to the new randomizations. The input file
        is checked for compatibility prior to starting.

        Once calculated, the data must be reread from the output file to lookup values
        from the PDF (see ``prepareLookup``).

        Arguments:

          outputHandle : output file handle
            handle for randomization output file used for writing randomizations;
            for each site, data written immediately after randomizations for that site performed.
            if ``None``, data is not written (for benchmarking)
          numRandomizations : number
            number of randomizations to perform
          inputHandle : input file handle, optional
            handle for randomization input file used for reading existing randomizations;
            these are combined with the (additional) randomizations performed here prior to
            writing; if ``None`` (default), no existing randomizations will be read

        """

        if outputHandle is not None:
          csvWriter = csv.writer(outputHandle)
          csvWriter.writerow([self.startMagicWord])

        if inputHandle:
          self.readFromFile(inputHandle) # read any existing randomizations
        self.totalRandomizations += numRandomizations
        reportMetadata = self.prepareMetadata()

        if verbose>=1: print 'Writing header...'
        for item in reportMetadata.items():
          if verbose>=2: print ': '.join([str(value) for value in item])
          if outputHandle is not None: csvWriter.writerow(item)

        if verbose>=1 : print 'Performing randomizations using %s score...'%ScoreType.reverse_mapping[self.scoreType]
        if verbose>=1 : print ' Total %i sequences, sampling %i sequences.'%(len(self.records),self.sampleSize)
        if verbose>=1 : print ' Total %i sites.'%self.M


        scores = numpy.zeros(numRandomizations, dtype=float64)
        allSamples = [None for x in xrange(numRandomizations)]
        allWeights = [None for x in xrange(numRandomizations)]

        seed = long(time.time() * 256)
        random.seed(seed)

        for randomization in xrange(numRandomizations):
          # move these out of inner loops below to speed up running (at cost of more memory)
          sample = randomSampleNoReplacement(self.records,self.sampleSize)
          weights = calculateWeights(sample, self.distanceMatrix, silent = True)
          allSamples[randomization] = sample
          allWeights[randomization] = weights

        rndStartTime = timeit.default_timer()

        for site in xrange(self.M):
            site1 = site + 1 # 1-based index, used for output, internally 0-based used
            scores[:] = numpy.nan
            iterStartTime = timeit.default_timer()

            if self.scoreType == ScoreType.VALDAR:
                # randomize with Valdar scoring
                for randomization in xrange(numRandomizations):
                    sample = allSamples[randomization]
                    weights = allWeights[randomization]
                    scores[randomization] = ValdarScore.valdarBySite(site, sample, self.subsMat, weights, silent = True)

            elif self.scoreType == ScoreType.PROFILE or self.scoreType == ScoreType.PROFILE_UNCORR:
                correctForSelfComparisons = (self.scoreType == ScoreType.PROFILE)
                for randomization in xrange(numRandomizations):
                    sample = allSamples[randomization]
                    weights = allWeights[randomization]
                    profileAtSite = SequenceProfile.generateProfileAtSite(site, sample, self.vectorSubsMat, weights, silent = True)
                    scores[randomization] = ProfileScore.scoreAgainstProfileBySite(site, sample, profileAtSite, self.vectorSubsMat, weights, correctForSelfComparisons = correctForSelfComparisons, silent = True)
            else: raise AssertionError("Invalid score type.")
            if verbose>=4:
              print "  Random scores:"
              print scores

            sitePDF = SitePDF.createFromScores(scores)
            iterEndTime = timeit.default_timer()
            if verbose>=1: print '  Completed randomizations for site %i, %.3f secs.'%(site1,iterEndTime-iterStartTime)
            if inputHandle:
              prevSitePDF = self.readNextSite()
              sitePDF = SitePDF.createByCombining((prevSitePDF,sitePDF))
              assert(site1==self.curSite)
            pointDensity = sitePDF.totalDensity()
            if self.totalRandomizations!=pointDensity:
              raise ValueError("Unexpected number of randomizations, have %i expect %i"%(pointDensity,self.totalRandomizations))

            if outputHandle is not None: csvWriter.writerow([str(item) for item in itertools.chain([site1],sitePDF.array)])
        if outputHandle is not None:
          csvWriter.writerow([self.endMagicWord])
          outputHandle.close()
        if inputHandle: inputHandle.close()
        rndEndTime = timeit.default_timer()
        if verbose>=1: print ' Completed all randomizations, %.3f secs.'%(rndEndTime-rndStartTime)

    def readFromFile(self,inputHandle, checkFile = False):
        """
        Read a file containing a saved PDF, checking for compatibility with ``self``. Data in file can
        then be read to check for consistency (if ``checkFile = True``),
        or read site-by-site using ``readNextSite`` to obtain the PDF for each site.

        PDF given by ``inputHandle`` is tested for compatibility with arguments provided to constructor.

        Arguments:

          inputHandle : input file handle
            handle for randomization input file
          checkFile : boolean, optional
            if ``True``, entire file is read and checked for compatibility and consistency, if function returns
            without throwing exceptions then the file passes these checks;
            if ``False`` (default), header is read and checked for compatibility, and then reading stopped
            prior to reading data, site-by-site reading can then be carried out using ``readNextSite`` to
            return the PDF for each site

        """
        inputHandle.seek(0) # reset file position
        self.csvReader = csv.reader(inputHandle)
        thisMetadata = self.prepareMetadata()

        self.filePosition = FilePosition.START

        for index,row in enumerate(self.csvReader):
          if len(row)==0: continue # skip blank lines
          if self.filePosition == FilePosition.START:
            if row[0] != self.startMagicWord:
              raise ValueError("Randomization input file (line %i) does not start with magic word %s - correct format?"%(index,self.startMagicWord))
            self.filePosition = FilePosition.HEADER
            fileMetadata = {}
          elif self.filePosition == FilePosition.HEADER:
            if len(row)<2:
              raise ValueError("Randomization input file header (line %i): Expected pairs of values - correct format?"%index)
            key = row[0]
            value = row[1]
            fileMetadata[key] = value
            if key == self.dataMagicWord:
              self.filePosition = FilePosition.DATA
              for key,expectedValue in thisMetadata.iteritems():
                if key not in fileMetadata:
                  raise ValueError("Randomization input file header (line %i): Value for key %s does not exist"%(index,key))
                value = fileMetadata[key]
                if value!=expectedValue and key != self.randomizationsHeader and key!=self.sequenceIdHeader:
                  raise ValueError("Randomization input file header (line %i): Value for key %s does not match expected (%s vs %s)"%(index,key,value,expectedValue))
              self.fileRandomizations = int(fileMetadata[self.randomizationsHeader]) # store number of read randomizations
              if not checkFile:
                # not checking file, i.e. reading to add to this
                self.totalRandomizations = self.fileRandomizations
                if verbose: print "Randomization input file has %i randomizations."%self.totalRandomizations
              self.curSite = 0 # 1-based, but incremented after successfully reading a site
              break
        if checkFile:
          if self.filePosition == FilePosition.DATA:
            while self.readNextSite(True):
              pass
            if self.filePosition != FilePosition.END:
              raise ValueError("Randomization input file (line %i) does not end with magic word %s- correct format?"%(index,self.endMagicWord))

    def readNextSite(self, checkFile = False):
        """
        Read a data from input file site-by-site. Should be called after ``readFromFile`` to allow header
        to be read and ensure file pointer is in the correct position.

        Arguments:

          checkFile : boolean, optional
            if ``True``, returns boolean depending of whether site can be read correctly, used internally
            by ``readFromFile(checkFile = True)``;
            if ``False`` (default), read the next site and return ``SitePDF`` object for the site

        Return:

          out : SitePDF object or boolean
            if ``checkFile==False``, ``SitePDF`` for next site, position of site given by ``self.curSite``;
            if ``checkFile==True``, ``True`` if site can be read correctly

        """
        if self.filePosition == FilePosition.DATA:
          while True:
            row = self.csvReader.next()
            if len(row)==0: continue # skip blank lines
            else: break
          if row[0] == self.endMagicWord:
            if self.curSite!=self.M:
              raise ValueError("Randomization input file data: Reached end of file, last site %i expected %i"%(self.curSite,self.M))
            self.filePosition = FilePosition.END
            return None
          try:
            fileSite = int(row[0])
            self.curSite+=1 # increment after successfully reading a site
            siteData = itertools.islice(row,1,None) # skip first element = site number
            if self.curSite!=fileSite:
              raise ValueError("Got site %i expected %i"%(fileSite,self.curSite))
            if checkFile:
              SitePDF.convertFromStrList(siteData,False) # check data is correct format
              return True
            else:
              sitePDF = SitePDF.convertFromStrList(siteData,True)
              pointDensity = sitePDF.totalDensity()
              if pointDensity != self.fileRandomizations:
                raise ValueError("Read unexpected number of randomizations, read %i expect %i"%(pointDensity,self.fileRandomizations))
              return sitePDF
          except Exception as e:
              raise ValueError("Randomization input file data: Error parsing data %s"%(str(e)))

        else:
          raise ValueError("Randomization input file data: Unexpected file position")

    def combineFromMultiple(self, outputHandle, inputHandles):
      csvWriter = csv.writer(outputHandle)
      csvWriter.writerow([self.startMagicWord])

      if verbose>=1 : print 'Combining randomizations from %i files...'%len(inputHandles)
      alignmentPDFs = []
      for index,inputHandle in enumerate(inputHandles):
        if verbose>=1 : print 'Checking file %i for compatibility...'%index
        alignmentPDF = AlignmentPDF(self.records, self.sampleSize, self.scoreType,
                                    self.vectorSubsMat if self.vectorSubsMat else self.subsMat,
                                    self.distanceMatrix)
        alignmentPDF.readFromFile(inputHandle, checkFile = True) # check whole file
        alignmentPDF.readFromFile(inputHandle) # reseek to data start position
        self.totalRandomizations += alignmentPDF.totalRandomizations
        if verbose>=1 : print 'Read %i randomizations from file %i...'%\
                              (alignmentPDF.totalRandomizations,index)
        alignmentPDFs.append(alignmentPDF)
      reportMetadata = self.prepareMetadata()

      if verbose>=1 : print 'Resulting file will have %i randomizations...'%self.totalRandomizations

      if verbose>=1: print 'Writing header...'
      for item in reportMetadata.items():
        if verbose>=2: print ': '.join([str(value) for value in item])
        if outputHandle is not None: csvWriter.writerow(item)

      for site in xrange(self.M):
          site1 = site + 1 # 1-based index, used for output, internally 0-based used
          if verbose>=1 : print 'Combining site %i...'%site1
          sitePDFs = []
          for index,alignmentPDF in enumerate(alignmentPDFs):
            if verbose>=3 : print 'Reading site %i from file %i...'%(site1,index)
            sitePDF = alignmentPDF.readNextSite()
            assert(site1==alignmentPDF.curSite)
            sitePDFs.append(sitePDF)
          sitePDF = SitePDF.createByCombining(sitePDFs)
          pointDensity = sitePDF.totalDensity()
          if self.totalRandomizations!=pointDensity:
            raise ValueError("Unexpected number of randomizations, have %i expect %i"%(pointDensity,self.totalRandomizations))
          csvWriter.writerow([str(item) for item in itertools.chain([site1],sitePDF.array)])
      csvWriter.writerow([self.endMagicWord])
      outputHandle.close()
      for index,inputHandle in enumerate(inputHandles):
        inputHandle.close()
      if verbose>=1: print ' Completed file combining.'

    def prepareLookup(self, inputHandle):
        """
        Prepare a file for looking-up values from the PDF site-by-site.
        PDF density for a site can then be found by calling ``lookupDensities``.

        Arguments:

          inputHandle : input file handle
            handle for randomization input file

        """
        self.readFromFile(inputHandle)

    def lookupDensities(self, site, score):
        """
        Lookup densities in a prepared PDF. Must be called after ``prepareLookup()``.

        Returns density in the PDF where function values ≥ ``score``,
        as well as the total function density. As data is obtained by reading from current file position,
        must be called in increasing order of site numbers.

        Arguments:

          site : number
            index (0-based) of site to lookup density for, must be greater than ``self.curSite``
          score : number
            score to return density for

        Return:

          out : tuple of (number,number)
            (density where function ≥ ``score``, total density of function)

        """
        site1 = site+1 # 1-based index, used in randomization file
        while True:
          # read correct site, may need to skip lines as not all sites looked up
          sitePDF = self.readNextSite()
          if verbose>=3: print "Read site %i from file."%self.curSite
          assert(site1>=self.curSite)
          if site1==self.curSite: break
        density = sitePDF.density(score,lambda a,b: operator.ge(a,b) or abs(a-b) < 1e-12) # greater or APPROX equal to allow for rounding error (important for comparisons to 1.0)
        totalDensity = sitePDF.totalDensity()
        if verbose>=2: print "Looked up density for site %i is %i/%i when random score>=%.4f."%(site1,density,totalDensity,score)
        assert(self.totalRandomizations==totalDensity)
        return (density, totalDensity)

class SitePDF:
    """
    Probability density function for a single site (in an alignment).
    Stored as an array of discrete points (``PDFPoint`` objects), each containing a value and a density.

    Initialized as an empty PDF (no points), use class methods to create a ``PDFPoint``:
      * From a an array of scores, using ``createFromScores()``.
      * By converting from a list of strings (e.g. read from a csv file) using ``convertFromStrList()``.
      * By combining two existing ``PDFPoint`` objects, usng ``createByCombining()``

    """
    __slots__ = ('array') # use slots to save memory
    def __init__(self):
        self.array = [] #: *list of PDFPoint objects*, array of points making up of this PDF

    @classmethod
    def createFromScores(cls, scores):
        """
        Create a new ``SitePDF`` object from an list/array of scores.

        Arguments:

            scores : list of numbers
                list of scores to convert

        Return:

            out : SitePDF object
                SitePDF containing points created from the supplied scores
        """
        this = cls()

        # sort scores and the summarize by collapsing equal scores into single points
        scores.sort()
        curPDFPoint = None

        for i,score in enumerate(scores):
            if curPDFPoint and curPDFPoint.approxEqual(score):
                curPDFPoint.increment() # same value, increment current point's density
            else:
                # different value, need new point
                curPDFPoint = PDFPoint(score,1)
                this.array.append(curPDFPoint)

        return this

    @classmethod
    def convertFromStrList(cls, strData, store = True):
        """
        Create a new ``SitePDF`` object by converting from a list of strings to a series
        of points.

        Format is [``'value1:density', 'value2:density2', ...]``,
        if ``:density`` is absent, density is assumed to be 1.

        This format for a point (``'value(:density)'``) is produced by ``str(pdfPoint)``.
        To produce a list of strings for all points in this format, use ``[str(item) for item in self.array]``.

        Arguments:

            strData : list of string
                string data to convert
            store : boolean, optional
                if ``True`` (default) read and store the data, otherwise just read without storing (for checking data format)

        Return:

            out : SitePDF object
                SitePDF containing points created from supplied data  or ``None`` if ``store`` is ``False``
        """
        if store: this = cls()
        else: this = None

        for element in strData:
            if ':' in element:
              valueStr, densityStr = element.split(':')
              density = int(densityStr)
            else:
              valueStr = element
              density = 1
            value = float64(valueStr)

            if store:
              curPDFPoint = PDFPoint(value,density)
              this.array.append(curPDFPoint)

        return this

    @classmethod
    def createByCombining(cls, sitePDFs):
        """
        Create a new ``SitePDF`` object by combining multiple ``SitePDF`` objects.

        Arguments:

            sitePDF : iterable of SitePDF objects
                objects to combine

        Return:

            out : SitePDF object
                SitePDF containing points created by combining the supplied objects
        """
        this = cls()
        combinedArray = reduce(lambda x, y: x+y, (sitePDF.array for sitePDF in sitePDFs))
        combinedArray.sort() # combine then sort points

        lastPoint = combinedArray [0]
        for i in range(1,len(combinedArray)): # iterate sorted points
            point = combinedArray[i]
            point.combineIfEqual(lastPoint) # check if equal value to previous point, if so, combine points leaving one point object with zero density
            lastPoint = point
        this.array = [point for point in combinedArray if point.density>0] # create new array excluding any points with zero density
        return this

    def density(self, score, oper):
      """
      Return summed density of this object for values that return ``True`` for ``oper(value,score)``.
      e.g. ``density(0.7,operator.ge)`` will return density where value is ≥ 0.7.
      IMPORTANT: Keep in mind that equality comparisons for floating point numbers may be unsafe. To
      ensure that comparisons to 1.0 (absolute conservation) succeed, use a greater than or 
      approximately equal operator (e.g. ``lambda a,b: operator.ge(a,b) or abs(a-b) < 1e-12``)


        Arguments:

            score : number
                comparison number
            oper : operator
                operator to use for comparison, e.g. ``operator.ge``

        Return:

          out : number
              summed density
      """
      return sum(point.density for point in self.array if oper(point.value,score))

    def totalDensity(self):
      """
      Sum total density of this object.

        Return:

          out : number
              total density in this object
      """
      return sum(point.density for point in self.array)

@total_ordering
class PDFPoint:
    """
    Create a single point of a probability density function containing a value and a density.
    Can be compared (``==, <`` etc.) to numbers and other ``PDFPoint`` objects to allow ordering
    and testing whether this ``PDFPoint`` represents same value.

    Note that values are considered equal if their difference is less than a very small value to allow
    for discrete nature of IEEE floating point representation (see ``approxEqual()``).

    Allows conversion to string in the format ``value:density`` or ``value`` (if ``density=1``) using
    ``str(pdfPoint)``.

    Arguments:

        value : number
            value of the point
        density : number, optional
            density of the point, stored as an integer indicating
            number of times this value obtained by randomization

    """
    __slots__ = ('value','density')
    def __init__(self, value, density = 0):
        self.value = value #: *number*, value of this point
        self.density = density #: *number*, density of this point (integer)

    def increment(self):
        """
        Increment the density of this point by one.
        """
        self.density += 1

    def combineIfEqual(self, other):
        """
        Combine this ``PDFPoint`` with another if the two have an (approximately) equal value.

        Arguments:

            other : PDFPoint
                other ``PDFPoint`` to possibly combine with

        """
        if self.approxEqual(other.value):
            self.density += other.density
            other.density = 0
    def checkComparable(self,other):
        """
        Check whether this ``PDFPoint`` can be compared to ``other``. Returns the value of
        ``other`` to compare against if ``other`` is another ``PDFPoint`` or same type as
        ``self.value``, returns ``None`` otherwise.

        Arguments:

            other : object
                other object to test

        """
        if type(other) is type(self):
            return other.value
        elif type(other) is type(self.value):
            return other
        else:
            return None
    def approxEqual(self, otherValue):
        """
        Return ``True`` if difference between ``self.value`` and ``otherValue`` is less than
        ``1e-12``, used to test equality given discrete nature of IEEE floating point
        representation and possible round-off error when calculating conservation scores.

        Arguments:

            otherValue : object
                other object to test

        """
        return abs(self.value-otherValue) < 1e-12
    def __str__(self):
        if self.density>1:
          return "%s:%i"%(str(self.value),self.density)
        else: return str(self.value)
    def __eq__(self, other):
        otherValue = self.checkComparable(other)
        if not otherValue: return False
        return self.approxEqual(otherValue) # return true if other is approx equal
    def __lt__(self, other):
        otherValue = self.checkComparable(other)
        if otherValue is None: raise TypeError('Comparing PDFPoint to unsupported type')
        return self.value<otherValue and not self.approxEqual(other.value)


##########################################################################################
# REPORT OUTPUT
##########################################################################################

def getAllSites(records):
    """
    Return a list of all sites in the alignment.

    Arguments:

      records : list of Bio.SeqRecord objects
        list of aligned sequences to process
    """
    M = len(records[0].seq) # alignment length
    allSites = list(i for i in range(M))
    return allSites

def getPresentSites(records):
    """
    Return a list of all present sites in the alignment.
    A present site must be represented by a non-gap character
    in at least one sequence.

    Arguments:

      records : list of Bio.SeqRecord objects
        list of aligned sequences to process
    """
    presentSites = []
    M = len(records[0].seq) # alignment length
    for site in range(M):
        sitePresent = False
        for record in records:
            if record.seq[site] != '-':
                sitePresent = True
                break
        if sitePresent:
            presentSites.append(site)
    return presentSites

class SequenceNumerator:
  """
  Iterator for sequence numbering.
  Skips from -1 to +1 (no 0).
  Allows subnumbering for gap positions, i.e. 1,1a,1b,1c.

  Arguments:

    start : number
      starting number
  """
  subnums = 'abcdefghijklmnopqrstuvwxyz' #: *string*, suffixes to use for subnumbering
  def __init__(self, start):
      self.reset(start)

  def __iter__(self):
      return self

  def reset(self,start):
    """
    Reset numbering to new starting position.

    Arguments:

      start : number
        starting number
    """
    self.num = None
    self.nextnum = start
    self.subnum = 0

  def next(self):
    """
    Return next (whole) number.
    """
    self.num  = self.nextnum
    self.nextnum = self.nextnum+1
    self.subnum = 0
    if self.nextnum==0: self.nextnum = 1 # no zeroes, -1 to +1
    return self.num

  def subnext(self):
    """
    Return next sub number, e.g. 1a, 1b etc.
    Must be called after a call to ``next()``.
    """
    if self.num is None: raise ValueError("SequenceNumerator cannot request a subnumber before a number requested")
    ret = "%i%s"%(self.num,self.__class__.subnums[self.subnum])
    self.subnum +=1
    return ret


def generateGenericReport(outputHandle, outputObject,
  records = None, alignmentPDF = None, referenceId = None,
  outputPositions = OutputPositions.PRESENT,
  simpleNumbering = False, numberingStart = 1,
  blockSequenceOutput = False, showCompleteProfile = False, pAsFraction = False):
    """
    Generate a conservation score report

    Arguments:

      path : string
        output path for the report
      outputObject : ConservationScores or SequenceProfile object
        scores or profile to output in the report
      records : list of Bio.SeqRecord objects, optional
        list of aligned sequences to output, use to override records stored in ``outputObject``
      alignmentPDF : ScoreProbDist object, optional
        probability distribution of the scores at each site, used to assign p-values, if absent not p-values written;
        must be prepared for lookup by precalling ``prepareLookup``
      referenceId : string, optional
        sequence record id to use as the reference sequence, used for numbering site positions, if absent first second used for numbering
      outputPositions : OutputPositions.ALL, OutputPositions.PRESENT or OutputPositions.REFERENCE, optional
        if ``OutputPositions.ALL`` scores for all positions written,
        if ``OutputPositions.PRESENT`` (default) only present positions written (i.e. at least one non-gap at site),
        if ``OutputPositions.REFERENCE`` only positions in reference sequence written (i.e. at least one non-gap at site),
      simpleNumbering : boolean, optional
        if ``True`` simple numbering (1,2,3...) used, otherwise (default) numbering is for reference sequence
      numberingStart : number, optional
        number used to start numbering, if absent numbering starts from 1
      blockSequenceOutput : boolean, optional
        if ``True`` sequences output as a block (column) of text, with one block for all scored sequences and
        one block for sequences not scored (but added to output), otherwise (default) each sequence is written to individual columns
      showCompleteProfile : boolean, optional
        if ``True`` complete sequence profile shown (i.e. score for each amino acid), otherwise (default) profile not shown; only valid for ``ProfileScore`` or ``SequenceProfile``
      pAsFraction : boolean, optional
        if ``True`` p-values produced as fraction e.g. 130/1000, otherwise (default) produced as a proportion
    """
    # get record names for header
    #recordNames = numberingBy

    if verbose>=0: print "Writing report..."

    reportMetadata = OrderedDict() # details about report
    reportMetadata.update(outputObject.metadata)

    sequenceProfile = None
    if isinstance(outputObject,SequenceProfile):
      sequenceProfile = outputObject
    if isinstance(outputObject,ProfileScore):
      sequenceProfile = outputObject.sequenceProfile
      reportMetadata.update(outputObject.sequenceProfile.metadata)

    if records is None: records = outputObject.records
    if referenceId is None: referenceId = records[0].id

    # determine index of record to number by
    referenceIndex = [index for index,record in enumerate(records) if record.id == referenceId ]
    if len(referenceIndex) == 0:
      raise Exception("Couldn't find reference sequence %s"%referenceId)
    else: referenceIndex = referenceIndex[0]
    if verbose>=1: print "Reference sequence is %s, index %i."%(referenceId,referenceIndex)

    referenceRecord = records[referenceIndex]

    # find output sites
    if outputPositions == OutputPositions.ALL:
        sitesOfInterest = getAllSites(records)
        reportMetadata['report site selections'] = 'all sites'
    elif outputPositions == OutputPositions.REFERENCE:
      sitesOfInterest = getPresentSites([records[referenceIndex]])
      reportMetadata['report site selections'] = 'sites in %s'%referenceId
    else:
      sitesOfInterest = getPresentSites(records)
      reportMetadata['report site selections'] = 'present sites'
    if verbose>=1: print "Writing report with %i sites."%(len(sitesOfInterest))

    reportMetadata['report sites'] = '%i sites'%len(sitesOfInterest)
    reportMetadata['report sequences shown'] = " ".join(record.id for record in records)


    # get numbering
    if simpleNumbering:
      sequenceNumerator = SequenceNumerator(numberingStart)
      sequenceNums = [sequenceNumerator.next() for site in sitesOfInterest]
      reportMetadata['report numbering by'] = 'simple numbering'
    else:
      # number by reference sequence

      for index,site in enumerate(sitesOfInterest):
        if referenceRecord.seq[site]!='-':
          break # find first non-gap residue in reference

      firstNonGap = site # number of first site
      firstNonGapInSOI = index # index of first site within sitesOfInterest
      firstNonGapNum = numberingStart

      if verbose>=3: print "First non-gap site is %i, numbered %i."%(site+1,numberingStart)

      # determine sequence numbering
      # number: -3, -2, -1, 1, 2, 3, 3a, 3b, 3c, 4
      # a,b,c used if there are gaps in reference sequence

      lastGapNum = numberingStart - 1 if numberingStart < 0 else -1 # sites up to first non-gap always given negative numbers
      sequenceNumerator = SequenceNumerator(lastGapNum+1-firstNonGapInSOI) # init numerator to first site number, so that sites up to first non gap -x... -2,-1
      sequenceNums = []
      for site in sitesOfInterest:
        if site==firstNonGap: sequenceNumerator.reset(firstNonGapNum) # jump to first numbering position (only has effect if this is >1)
        if site<firstNonGap or referenceRecord.seq[site]!='-':
          nextNum = sequenceNumerator.next()
        else:
          nextNum = sequenceNumerator.subnext()
        sequenceNums.append(nextNum)
      if verbose>=4:
        print "Site numbering (index, refseq, site num):"
        for index,siteNum in enumerate(sequenceNums): print index+1,referenceRecord.seq[sitesOfInterest[index]],siteNum
      reportMetadata['report numbering by'] = 'reference sequence %s starting from %i'%(referenceId,numberingStart)
    reportMetadata['report numbering extent'] = '%s to %s'%(str(sequenceNums[0]),str(sequenceNums[-1]))

    reportRows = [OrderedDict() for site in sitesOfInterest]

    reportMetadata['note sequences'] = "sequences marked * were not scored"

    refSeqHeader = '%sRef: %s'%('' if referenceRecord in outputObject.records else '*', referenceId)

    # add data for numbering
    for siteIndex, site in enumerate(sitesOfInterest):
      reportRow = reportRows[siteIndex]
      reportRow['alignment num'] = site+1 #1-based
      reportRow['report num'] = sequenceNums[siteIndex]

    # add data for reference seq
    for siteIndex, site in enumerate(sitesOfInterest):
      reportRow = reportRows[siteIndex]
      reportRow[refSeqHeader] = referenceRecord.seq[site]

    # add data for other sequences
    if blockSequenceOutput:
      # write sequences as a block of text
      recordsNotInScoreHeader = '*Not scored:'
      recordsInScoreHeader = 'Scored:'

      indexS = 1
      indexNS = 1
      # generate headers
      for record in records:
        if record in outputObject.records:
          # write scored sequence names to header
          recordsInScoreHeader += ' %i=%s'%(indexS,record.id)
          indexS += 1
        else:
          # write not-scored sequence names to header
          recordsNotInScoreHeader += ' %i=%s'%(indexNS,record.id)
          indexNS += 1

      # write sequences
      for siteIndex, site in enumerate(sitesOfInterest):
        reportRow = reportRows[siteIndex]

        for record in records:
          if record in outputObject.records:
            data = reportRow.get(recordsInScoreHeader,'')
            data += record.seq[site]
            reportRow[recordsInScoreHeader] = data

          else:
            data = reportRow.get(recordsNotInScoreHeader,'')
            data += record.seq[site]
            reportRow[recordsNotInScoreHeader] = data

    if not blockSequenceOutput: # write sequences as individual columns
      # write sequences
      for siteIndex, site in enumerate(sitesOfInterest):
        reportRow = reportRows[siteIndex]

        for record in records:
          # write not scored first
          if record not in outputObject.records: keyname = '*'+record.id # add * to indicate not scored
          else: keyname = record.id
          while keyname in reportRow: keyname+=" " # deal with duplicates by adding spaces to name

          reportRow[keyname] = record.seq[site]

    # output sequence profile
    if sequenceProfile:
      if showCompleteProfile:
        reportMetadata['note profile'] = "showing complete profile vector, score for each amino acid listed according to one letter code"

      if sequenceProfile.consensus:
          reportMetadata['note consensus'] = "showing consensus determined from profile"
          reportMetadata['consensus details'] = sequenceProfile.consensusDetails

      letters = sequenceProfile.vectorSubsMat.alphabet.letters
      for siteIndex, site in enumerate(sitesOfInterest):
        reportRow = reportRows[siteIndex]
        if showCompleteProfile:
          scoresAtSite = sequenceProfile[site]
          for letter in letters:
              letterIndex = sequenceProfile.vectorSubsMat.alphabetLUT[ord(letter)]
              reportRow[letter] = repr(scoresAtSite[letterIndex])
        if sequenceProfile.consensus:
            consensusAtSite = sequenceProfile.consensus[site]
            for i,aa in enumerate(consensusAtSite):
              consensusHeader = 'consensus %i'%(i+1)
              reportRow[consensusHeader] = aa.letter
              reportRow[consensusHeader+" score"] = repr(aa.score)

    # output scores
    if isinstance(outputObject,ConservationScores):
      conservationScores = outputObject
      reportMetadata['note scores'] = "showing score for each position under score"
      if alignmentPDF:
        #reportMetadata['note p-value'] = "showing p-value under pvalue, determined from %i randomizations"%alignmentPDF.totalRandomizations
        reportMetadata['note p-value'] = "showing p-value from randomizations under pvalue"
        for key,value in alignmentPDF.prepareMetadata().iteritems(): # add randomization metadata
          if key!=AlignmentPDF.dataMagicWord:
            reportMetadata['randomization '+key] = value
      for siteIndex, site in enumerate(sitesOfInterest):
        reportRow = reportRows[siteIndex]
        score = conservationScores.scores [site]
        reportRow["score"] = repr(score)
        if alignmentPDF:
          density,totalDensity = alignmentPDF.lookupDensities(site,score)
          if pAsFraction: reportRow["pvalue"] = "%i/%i"%(density,totalDensity)
          else: reportRow["pvalue"] = str(density*1.0/totalDensity)

    # write metadata
    csvWriter = csv.writer(outputHandle)
    reportMetadata['note header ends'] = "data follows immediately below"
    if verbose>=1: print 'Writing metadata...'
    for item in reportMetadata.items():
      if verbose>=2: print ': '.join([str(value) for value in item])
      csvWriter.writerow(item)

    # write rows
    if verbose>=1: print 'Writing data...'

    referenceRow = None
    for reportRow in reportRows: # use largest row as reference row
      if not referenceRow or len(reportRow)>len(referenceRow):
        referenceRow = reportRow


    csvDictWrite = csv.DictWriter(outputHandle, referenceRow)
    if verbose>=4: print 'Reference row (headings) has %i columns.'%len(referenceRow)
    csvDictWrite.writeheader()
    for index,reportRow in enumerate(reportRows):
      if verbose>=4: print 'Row %i has %i columns.'%(index+1,len(reportRow))
      csvDictWrite.writerow(reportRow)

    outputHandle.close()
    return



#@profile
def main():
  """
  Execute main script code function.

  Verbosity levels:

    0/None : default, main messages

    1 : a few extra messages showing progress

    2 : outputs top level calculations

    3 : outputs mid-level calculations

    4 : outputs detailed calculations
  """

  global verbose
  global version

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                   description='Calculates sequence conservation scores',
                                   epilog="")
  parser.add_argument("infile",
    help="input alignment to process, fasta format")

  parser.add_argument("-v", "--verbose", action="count",help="produce verbose output, increase number of times specified to increase verbosity")
  parser.add_argument('--version', action='version', version='%(prog)s '+"%s"%version)
  parser.add_argument('--overwrite', action='store_true', help="force overwriting of output files")

  scoreGroup = parser.add_argument_group('score parameters', 'parameters for score calculation')

  scoreGroup.add_argument("-s", "--subset", default=".*",
    help="subset of sequences to analyse, selected by regular expression search for SUBSET")
  scoreGroup.add_argument("-t", "--type", default="PROFILE",
    help="type of conservation score calculated, VALDAR for Valdar score, PROFILE for (corrected) profile-based score, PROFILE_UNCORR for uncorrected profile-based score",
    choices=("VALDAR","PROFILE","PROFILE_UNCORR"))
  scoreGroup.add_argument("-m","--subsmat", default="Valdar2001",
    help="substitution matrix to use, aaindex2 format")
  scoreGroup.add_argument("--diagonal",type=float,
    help="score to use for diagonal elements in substitution matrix (identities), set prior to matrix transformation, if none the rounded average of diagonal elements in the unmodified matrix used")
  scoreGroup.add_argument("--gap",type=float,default=0.0,
    help="score to use for gaps, set after matrix transformation")
  scoreGroup.add_argument("--gapgap",type=float,
    help="score to use for gap-gap substitutions, set after matrix transformation, if absent GAP score used")
  scoreGroup.add_argument("--uniform",action='store_true',
    help="if present, apply weighting to all sequences, otherwise sequences weighted by evolutionary distance")

  consGroup = parser.add_argument_group('consensus parameters', 'parameters for determining consensus residues (profile-based scores only)')

  consGroup.add_argument("-c","--cutoff", type=float,
    help="score cutoff for consensus residues, if absent consensus not determined")
  consGroup.add_argument("--propcutoff", type=float,default=1.0,
    help="proportional score cutoff for consensus residues consensus may also be this proportion of the highest scoring residue")

  rndGroup = parser.add_argument_group('randomization parameters', 'parameters for randomizations')

  rndGroup.add_argument("-n","--numrnd", type=int,default=0,
    help="number of (additional) randomizations, if none no randomizations performed, must be none if combining multiple input files")
  rndGroup.add_argument("--rndsubset", default=".*",
    help="subset of sequences to randomize from, selected by regular expression search for RNDSUBSET, default is all sequences")
  rndGroup.add_argument("-i","--rndinfile", nargs='*',
    help="randomization input file(s) for reading previously performed randomizations, if multiple files specified these will be combined into output file")
  rndGroup.add_argument("-p","--rndoutfile",
    help="randomization output file, required when performing additional randomizations or combining multiple input files, must differ from RNDINFILE")

  outputGroup = parser.add_argument_group('output parameters', 'parameters for generating score report')


  outputGroup.add_argument("-o","--outfile",
    help="output file to produce, csv format, if not present no report generated")
  outputGroup.add_argument("-r","--refid",
    help="reference sequence for numbering, if absent first sequence in subset used")
  outputGroup.add_argument("--outsubset",
    help="subset of sequences to output in report, selected by regular expression search for OUTSUBSET, if absent score SUBSET used")
  outputGroup.add_argument("--outpos", default="PRESENT",
    help="positions to output in report, PRESENT for non-gap postions, ALL for all positions or REFERENCE for positions in reference sequence only",
    choices=("PRESENT","ALL","REFERENCE"))
  outputGroup.add_argument("--simplenum", action='store_true',
    help="if present use simple sequential numbering of sequence positions (1,2,3...), otherwise numbering is by reference sequence")
  outputGroup.add_argument("--numstart", type=int, default=1,
    help="start sequence numbering at this value")
  outputGroup.add_argument("--blockoutput", action='store_true',
    help="if present output sequences in as a single block of text, otherwise default to one column per sequence")
  outputGroup.add_argument("--completeprofile", action='store_true',
    help="if present output complete sequence profile, i.e. profile scores for all positions")
  outputGroup.add_argument("--pfraction", action='store_true',
    help="if present output p-values as fractions instead of proportions")


  args = parser.parse_args()

  print ""
  print "Parameters:"

  verbose = args.verbose
  overwrite = args.overwrite
  print " Input alignment file is %s."%(args.infile)
  if verbose: print " Verbosity level %i and will %soverwrite files."%(verbose,"" if overwrite else "NOT ")

  # get score arguments
  scoreOptions = Dummy()
  scoreOptions.subsetRE = args.subset
  scoreOptions.scoreType = getattr(ScoreType,args.type)
  scoreOptions.subsMat = args.subsmat
  scoreOptions.gap = float(args.gap)
  scoreOptions.gapgap = float(args.gapgap) if args.gapgap is not None else scoreOptions.gap
  scoreOptions.diag = args.diagonal
  if scoreOptions.diag: scoreOptions.diag = float(scoreOptions.diag) # convert to float if not none
  scoreOptions.uniform = args.uniform

  print " Calculating scores from subset \"%s\" using %s score."%(scoreOptions.subsetRE,ScoreType.reverse_mapping[scoreOptions.scoreType])
  print " Using subsitution matrix %s."%(scoreOptions.subsMat)
  print " Score for gaps is %.2f, gap-gap is %.2f and diagonal is %s."%(scoreOptions.gap,scoreOptions.gapgap,"%.2f"%scoreOptions.diag if scoreOptions.diag else "average")
  print " %s."%("Uniform weighting of all sequences" if scoreOptions.uniform else "Weighting by evolutionary distance")

  # get consensus arguments
  consOptions = None
  if isProfileScore(scoreOptions.scoreType):
    if args.cutoff is not None:
      consOptions = Dummy()
      consOptions.cutoffAbs = float(args.cutoff)
      consOptions.cutoffProp = float(args.propcutoff)

  if consOptions:
    print " Consensus residues must have score of >=%.2f or >=%.2f of highest score."%(consOptions.cutoffAbs,consOptions.cutoffProp)
  else: print " Not determining consensus residues."

  # get randomization arguments
  rndOptions = Dummy()
  rndOptions.number = args.numrnd
  rndOptions.infile = args.rndinfile
  rndOptions.outfile = args.rndoutfile
  rndOptions.subsetRE = args.rndsubset

  print " %i randomizations requested from subset \"%s\"."%(rndOptions.number,rndOptions.subsetRE)
  if rndOptions.infile:
    print " Using randomization input files %s."%(" ".join(rndOptions.infile))
  if rndOptions.outfile:
    print " Using randomization output file %s."%(rndOptions.outfile)

  # get report arguments
  reportOptions = Dummy()
  reportOptions.outfile = args.outfile
  reportOptions.refSeqID = args.refid
  reportOptions.subsetRE = args.outsubset if args.outsubset else scoreOptions.subsetRE
  reportOptions.outputPositions = getattr(OutputPositions,args.outpos)
  reportOptions.simpleNumbering = args.simplenum
  reportOptions.numStart = int(args.numstart)
  reportOptions.blockOutput = args.blockoutput
  reportOptions.showCompleteProfile = args.completeprofile
  reportOptions.pAsFraction = args.pfraction
  reportOptions.rndlookupfile = args.rndoutfile
  if reportOptions.rndlookupfile is None and args.rndinfile is not None:
   reportOptions.rndlookupfile = args.rndinfile[0]
  if rndOptions.outfile:
    print " Using randomization lookup file %s."%(rndOptions.outfile)
  if not isProfileScore(scoreOptions.scoreType):
    reportOptions.showCompleteProfile = None # only profile scores can use this option
  if reportOptions.outfile:
    print " Writing report to %s."%(reportOptions.outfile)
    print " Will write scores from subset \"%s\" as %s with %s as reference sequence."%\
      (reportOptions.subsetRE, "a block" if reportOptions.blockOutput else "columns", reportOptions.refSeqID if reportOptions.refSeqID else "first sequence")
    print " Numbering from %i using %s numbering."%\
      (reportOptions.numStart,"simple" if reportOptions.simpleNumbering else "reference sequence")
    print " Will write %s positions%s%s."%\
      (OutputPositions.reverse_mapping[reportOptions.outputPositions], ", consensus residues" if consOptions else "",
      ", complete profile" if reportOptions.showCompleteProfile else "")
    if reportOptions.rndlookupfile: print " p-values will be looked up from randomization file %s and shown as %s."\
      %(reportOptions.rndlookupfile,"fractions" if reportOptions.pAsFraction else "proportions")

  #sys.exit()

  print ""
  print "Reading/opening files:"

  # read input alignment
  ALLrecords = [] # stores alignment records

  try:
    print "Reading input alignment %s..."%args.infile
    if not os.path.exists(args.infile):
        msg = "Input alignment file does not exist."
        print msg
        raise ValueError(msg)
    infile = open(args.infile, "rU")
    for index,record in enumerate(Bio.SeqIO.parse(infile, "fasta")):
        ALLrecords.append(record)
        if verbose>=2: print "Read sequence %i: %s"%(index+1,record.id)
    infile.close()
    print "Read input alignment containing %i sequences."%len(ALLrecords)
    if verbose:
      print "List of input sequences:"
      printRecordNames(ALLrecords)

  except Exception as e:
    print "Error reading input file!"
    if verbose:
      print traceback.format_exc()
    sys.exit(-1)

  # ensure output file is writable
  outfile = None
  try:
    if reportOptions.outfile:
      print "Opening output file %s..."%reportOptions.outfile
      if not overwrite and os.path.exists(reportOptions.outfile) and os.stat(reportOptions.outfile).st_size>0:
          msg = "Output file exists and contains data. Will not replace."
          print msg
          raise ValueError(msg)
      outfile = open(reportOptions.outfile, "wb")
      print "Output file opened."
  except Exception as e:
    print "Error opening output file!"
    if verbose:
      print traceback.format_exc()
    sys.exit(-1)

  # read and transform subs matrix
  try:
    #print "Reading substitution matrx %s..."%scoreOptions.subsMat
    subsMat = readSubsMat(scoreOptions.subsMat)
  except Exception as e:
    print "Error reading substitution matrix!"
    if verbose:
      print traceback.format_exc()
    sys.exit(-1)


  # read input randomization file
  rndinfile = None
  rndinfiles = None
  try:
    if rndOptions.infile:
      rndinfiles = []
      for path in rndOptions.infile:
        print "Opening randomization input file %s..."%path
        if not os.path.exists(path):
            msg = "Randomization input file does not exist."
            print msg
            raise ValueError(msg)
        rndinfile = open(path, "rb")
        rndinfiles.append(rndinfile)
      if len(rndinfiles)==1:
        rndinfiles = None # use rndinfile
      elif len(rndinfiles)>1:
        rndinfile = None # use rndinfiles
        if rndOptions.number>0:
          msg = "Randomizations requested and multiple input files provided: unsupported."
          print msg
          raise ValueError(msg)
  except Exception as e:
    print "Error opening randomization input file(s)!"
    if verbose:
      print traceback.format_exc()
    sys.exit(-1)

  # ensure output randomization file is writable
  rndoutfile = None
  try:
    if rndOptions.number>0 or rndinfiles:
      if not rndOptions.outfile:
        msg = "Randomizations or combining input files requested and no randomization output file provided."
        print msg
        raise ValueError(msg)
      print "Opening randomization output file %s..."%rndOptions.outfile
      if not overwrite and os.path.exists(rndOptions.outfile) and os.stat(rndOptions.outfile).st_size>0:
          msg = "Randomization output file exists and contains data. Will not replace."
          print msg
          raise ValueError(msg)
      rndoutfile = open(rndOptions.outfile, "wb")
  except Exception as e:
    print "Error opening randomization output file!"
    if verbose:
      print traceback.format_exc()
    sys.exit(-1)


  # Calculate scores for sequence subset
  print ""
  print "Processing:"

  print "Subsetting sequences by score regex %s..."%scoreOptions.subsetRE
  scoreSubsetRecords = subsetByID(ALLrecords, scoreOptions.subsetRE)

  print "Got a subset of %i sequences."%len(scoreSubsetRecords)
  if verbose:
    print "List of sequences in subset:"
    printRecordNames(scoreSubsetRecords)

  print "Transforming substitution matrix..."
  linearTransform(subsMat,gapgap=scoreOptions.gapgap,gap=scoreOptions.gap,diagonal=scoreOptions.diag)
  if isProfileScore(scoreOptions.scoreType):
    print "Converting substitution matrix to vector form..."
    vectorSubsMat = VectorSubsMat(subsMat)

  if scoreOptions.uniform:
    # use simple uniform distances
    distanceMatrix = createUniformDistanceMatrix(ALLrecords,1.0)
  else:
    # find full distance matrix, i.e. pairwise distance of all sequences
    print "Calculating pairwise evolutionary distances..."
    distanceMatrix = calculateDistanceMatrix(ALLrecords,subsMat)

  print "Calculating scores for subset..."

  # find weights for each sequence from distance matrix
  print "Determining sequence weights for subset..."
  weights = calculateWeights(scoreSubsetRecords, distanceMatrix)

  if scoreOptions.scoreType == ScoreType.VALDAR:
    print "Determining Valdar conservation scores for subset..."
    theScores = ValdarScore(scoreSubsetRecords, subsMat, weights = weights)
  elif scoreOptions.scoreType == ScoreType.PROFILE or scoreOptions.scoreType == ScoreType.PROFILE_UNCORR:
      print "Determining sequence profile for subset..."
      theSequenceProfile = SequenceProfile(scoreSubsetRecords, vectorSubsMat, profileWeights = weights)
      if scoreOptions.scoreType == ScoreType.PROFILE:
        print "Determining corrected profile-based conservation scores for subset..."
        correctForSelfComparisons = True
      else:
        print "Determining uncorrected profile-based conservation scores for subset..."
        correctForSelfComparisons = False
      theScores = ProfileScore(scoreSubsetRecords, theSequenceProfile , scoringWeights = weights, correctForSelfComparisons=correctForSelfComparisons)
  else: raise AssertionError("Invalid score type.")

  if verbose>=3:
    for site,score in enumerate(theScores.scores):
      print "Site %i: score %.4f"%(site+1,score)

  if consOptions:
    print "Determining consensus residues..."
    theSequenceProfile.findConsensus(consOptions.cutoffAbs, consOptions.cutoffProp)

  # Do randomizations
  alignmentPDF = None
  if rndOptions.number>0 or rndinfile or rndinfiles:
    print ""
    print "Determining score PDFs:"
    sampleSize = len(scoreSubsetRecords)
    print "Subsetting sequences by randomization regex %s..."%rndOptions.subsetRE
    rndSubsetRecords = subsetByID(ALLrecords, rndOptions.subsetRE)
    if verbose:
      print "List of sequences in randomization subset:"
      printRecordNames(rndSubsetRecords)
    print "Alignment PDFs from randomly sampling %i sequences from %i sequences..."%(sampleSize,len(rndSubsetRecords))
    rndSubsMat = vectorSubsMat if isProfileScore(scoreOptions.scoreType) else subsMat
    alignmentPDF = AlignmentPDF(records = rndSubsetRecords, sampleSize = sampleSize, scoreType = scoreOptions.scoreType, subsMat = rndSubsMat, distanceMatrix = distanceMatrix)
    if rndinfile:
      print "Checking randomization input file..."
      alignmentPDF.readFromFile(rndinfile, True)
    if rndOptions.number>0:
      print 'Starting randomizations...'
      alignmentPDF.randomize(rndoutfile,rndOptions.number,rndinfile)
    elif rndinfiles:
      # combine multiple files
      alignmentPDF.combineFromMultiple(rndoutfile,rndinfiles)


  # Write report
  if outfile:
    print ""
    print "Preparing report:"

    if reportOptions.subsetRE==scoreOptions.subsetRE:
      reportSubsetRecords = scoreSubsetRecords
    else:
      print "Subsetting sequences by report regex %s..."%reportOptions.subsetRE
      reportSubsetRecords = subsetByID(ALLrecords, reportOptions.subsetRE)
      if verbose:
        print "List of sequences in output subset:"
        printRecordNames(reportSubsetRecords)

    if alignmentPDF:
      print "Preparing randomization file %s for lookup:"%reportOptions.rndlookupfile
      lookupFile = open(reportOptions.rndlookupfile, "rb")
      alignmentPDF.prepareLookup(lookupFile)

    generateGenericReport(outfile, theScores, records = reportSubsetRecords, alignmentPDF = alignmentPDF, referenceId = reportOptions.refSeqID,
      outputPositions = reportOptions.outputPositions, simpleNumbering = reportOptions.simpleNumbering, numberingStart = reportOptions.numStart,
      blockSequenceOutput = reportOptions.blockOutput, showCompleteProfile = reportOptions.showCompleteProfile, pAsFraction = reportOptions.pAsFraction)
  else:
    print "No report requested. Did you forget the -o option?"

  print "Complete!"

verbose = 0
version = "v1.0.9"

if __name__ == "__main__":
    main()