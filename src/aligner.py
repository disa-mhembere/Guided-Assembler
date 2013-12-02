#!/usr/bin/python

# aligner.py
# Created on 2013-11-26.
# Email: disa@jhu.edu, slee320@jhu.edu
# Copyright (c) 2013. All rights reserved.

import dynamic_bwt
import reference
import numpy as np
import random

class Aligner():

  def __init__(self, ref_str):
    """
    Object used for alignment
    """
    self.dbwt = dynamic_bwt.BWT(ref_str) # contains SA, LF, get_LF
    self.ref = reference.Reference(ref_str)


  def align(self, read):
    """
    Find an approximate match of the read in the reference

    @param read: Is the read 
    """ 
    T = self.ref.R

    K = 3
    tol = 5
    parts,poffs = partition(read,K)

    minEditDist = 10**6
    minIndels = 10**6
    cutoff = []
    transcripts = []

    hits = [3,8,11]
    poff = 0
    partlen = len(read)/K
    for part in parts:    # Go through partitions to find matches
      for hit in hits:    # Query each partition exactly
        left = max(0,hit-poff-tol)
        right = min(len(T),hit-poff+len(read)+tol)

        #Find approximate match
        M,cut,transcript,mindels = kEdit(read,T[left:right],tol)

        if M == -1:   # This happens if no match found below tolerance
          continue

        if M < minEditDist:   # Search for one with min edits
          minEditDist = M
          minIndels = mindels
          cutoff = [cut+left]
          transcripts = [transcript]
        elif M == minEditDist:      
          if mindels < minIndels:   # Trying to reduce # of indels
            minIndels = mindels
            cutoff = [cut + left]
            transcripts = [transcript]
          elif mindels == minIndels:
            if cut+left not in cutoffs:   # Avoid already seen ones
              cutoff.append(cut+left)
              transcripts.append(transcript)
      poff += partlen

    if len(cutoff)>1:    # If multiple choices, just get random one
      idx = random.choice(range(0,len(cutoff),1))
      cutoff = cutoff[idx]
      transcripts = transcripts[idx]
    else:
      cutoff = cutoff[0]
      transcripts = transcripts[0]

    print "READ:",read
    print "REF:",T
    print cutoff
    print transcripts


    ### NOW UPDATE REFERENCE...
    j = 0
    for nt in transcripts:
      self.ref.match(cutoff,read[j])
      cutoff+=1
      j+=1

    print self.ref.match_count




  def alter_bwt(self,):
    raise NotImplementedError("Alter BWT as necessary")



def partition(read, k):
  """
  Partition the read into k parts

  @param read: Is the read to partition
  @param k: Number of parts for the partition
  """

  assert len(read)%k is 0
  partlen = len(read)/k
  parts = []
  poffs = []
  for i in xrange(k):
    parts.append(read[i*partlen:(i+1)*partlen])
    poffs = poffs + [i*partlen]
  return parts, poffs

def kEdit(p,t,k):
  """
  Find an approximate match of p in t with up to k edits

  @param p: query string
  @param t: reference string
  @param k: max number of edits
  @param M: minimum edit distance
  @param spots: offsets in t where edits hit
  """

  m = len(p)+1      # Number of rows in DP matrix
  n = len(t)+1      # Number of cols in DP matrix

  #Initialize DP matrix
  D = [[[0,'']]*n for i in xrange(m)]
  for i in range(m)[1:]:
    D[i][0]=[i,'I']
  for j in range(n)[1:]:
    D[0][j] = [0,'D']

  # Dynamic programming through the matrix
  for i in range(m)[1:]:
    for j in range(n)[1:]:

      # s indicates a character match or not (0 is match)
      s = 1
      if p[i-1] == t[j-1]:
        s = 0

      d = min([D[i-1][j][0]+1,D[i][j-1][0]+1,D[i-1][j-1][0]+s])

      if d == D[i-1][j][0]+1:
        ch = 'I'
      if d == D[i][j-1][0]+1:
        ch = 'D'
      if d == D[i-1][j-1][0]+s:     # Write down transcript chars
        if s == 0:
          ch = 'M'
        else:
          ch = 'R'
      D[i][j] = [d,ch]

  # BACKTRACE

  finrow = []
  for blob in D[-1]:
    finrow = finrow + [blob[0]]
  M = min(finrow)

  if M > k:
    return -1,[],[],[]

  minIndex = []
  for r in xrange(n):
    if finrow[r]==M:
      minIndex.append(r)

  minIndels = 10**6
  optP = ''

  wackyChanges = {"A":"B","C":"D","G":"H","T":"U"}

  for spot in minIndex:
    row = m-1
    col = spot
    ep = ''

    indels = 0

    while row > 0 or col > 0:     # Trace back through matrix

      c = D[row][col][1]

      if c == 'I':
        ep = wackyChanges[p[row-1]] + ep
        row = row - 1
        indels += 1
      elif c == 'D':
        ep = '-' + ep
        col = col - 1
        indels += 1
      else:
        ep = p[row-1] + ep
        row = row - 1
        col = col - 1

    if indels < minIndels:
      minIndels = indels
      optP = ep

  del D

  cutoff = 0
  for i in xrange(len(optP)):
    if optP[i]=="-":
      cutoff+=1
    else:
      break

  return M,cutoff,optP[cutoff:],minIndels


def test():
  aligner = Aligner("ACTGTTGGAAAACCTTGTTGTACCCGGGTTAAACCCCC")
  aligner.align("CCTTGTTGT") # Should be exact match of idx 11 if read lengths are 5

if __name__ == "__main__":
  test()