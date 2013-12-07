#!/usr/bin/python

# aligner.py
# Created on 2013-11-26.
# Email: disa@jhu.edu, slee320@jhu.edu
# Copyright (c) 2013. All rights reserved.

import dynamic_bwt
import reference
import numpy as np
import random
import Queue
import assembler
import math
import reference
class Aligner():

  def __init__(self, ref_str, target):
    """
    Object used for alignment of a target object that can generate multiple reads
    related to target string. i.e. reads to be aligned

    @param ref_str: the reference string
    @param target: the target object (Not string)
    """
    self.dbwt = dynamic_bwt.dBWT(ref_str) # contains SA, LF, get_LF
    self.ref = reference.Reference(ref_str)
    self.targ = target

    self.tol = 0      # Set error tolerance here...

    self.updates = Queue.Queue()


  def align(self):
    """
    Find an approximate match of the read in the reference.
    Pulls reads from the target_reads class.
    """
    T = self.ref.R
    read = self.targ.get_read()

    K = 2
    parts,poffs = partition(read,K)

    minEditDist = 10**6
    minIndels = 10**6
    cutoff = []
    transcripts = []

    poff = 0
    partlen = len(read)/K
    for part in parts:    # Go through partitions to find matches
      hits = self.dbwt.match(str(part))
      if not hits:
        continue
      if None in hits: pdb.set_trace()
      for hit in hits:    # Query each partition exactly
        left = max(0,hit-poff-self.tol)
        right = min(len(T),hit-poff+len(read)+self.tol)

        #Find approximate match
        M,cut,transcript,mindels = kEdit(read,T[left:right],self.tol)

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
            if cut+left not in cutoff:   # Avoid already seen ones
              cutoff.append(cut+left)
              transcripts.append(transcript)
      poff += partlen

    if len(cutoff)>1:    # If multiple choices, just get random one
      idx = random.choice(range(0,len(cutoff),1))
      cutoff = cutoff[idx]
      transcripts = transcripts[idx]
    elif not cutoff:      # If NO choices, return flags to SKIP
      return None
    else:
      cutoff = cutoff[0]
      transcripts = transcripts[0]

    if minEditDist == 0: self.targ._update_seen(cutoff)




    ### NOW UPDATE REFERENCE...
    lenny = len(self.ref.R)
    for nt in transcripts:
      if cutoff < lenny:
        if self.ref.match(cutoff,nt): self.updates.put([cutoff,nt])
      if nt not in "BDHU":
        cutoff += 1


  def alter_bwt(self,no_opt):
    """
    Alters the BWT dynamically based on the updates pending in the
    queue created in the align() method.

    @param no_opt: boolean if true do no optimizations i.e use a static bwt & no dbwt
    """
    changes = {"B":"A","D":"C","H":"G","U":"T"}

    while self.updates.qsize() > 0:
      T = self.ref.R
      up = self.updates.get()   #up = [spot,nt]
      spot = up[0]
      ch = up[1]

      #CASE: NO OPTIMIZATION
      if no_opt:
              # MAKE APPROPRIATE BWT UPDATES
        if ch in "ACGT" and T[spot] != ch:
          T = T[0:spot] + ch + T[spot+1:]
          print "SWAP IN",ch,"AT",spot
        elif ch in "BDHU":
          T = T[0:spot] + changes[ch] + T[spot:]
          self.ref.insert_idx(spot)
          print "INSERT",changes[ch],"AT",spot
        else:
          T = T[0:spot] + T[spot+1:]
          self.ref.del_idx(spot)
          print "DELETEE",spot
        self.dbwt = dynamic_bwt.dBWT(T)
        self.ref.R = T

      #CASE: OPTIMIZATION
      else:
      # MAKE APPROPRIATE BWT UPDATES
        if ch in "ACGT" and T[spot] != ch:
          T = T[0:spot] + ch + T[spot+1:]
          self.dbwt.replace_one(ch,spot)
          print "SWAP IN",ch,"AT",spot
        elif ch in "BDHU":
          T = T[0:spot] + changes[ch] + T[spot:]
          self.ref.insert_idx(spot)
          self.dbwt.insert_one(changes[ch],spot)
          print "INSERT",changes[ch],"AT",spot
        else:
          T = T[0:spot] + T[spot+1:]
          self.ref.del_idx(spot)
          self.dbwt.delete_one(spot)
          print "DELETEE",spot
        self.ref.R = T


def partition(read, k):
  """
  Partition the read into k parts

  @param read: Is the read to partition
  @param k: Number of parts for the partition
  """

  assert len(read)%k == 0, "(Read length %% %d) cannot equal zero!" %k
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
  from target_reads import Target
  targ = Target(p=0.03, read_length=6, T="ACTGTTGGGAGCCTTGTTGTCCAGGGTTAAACCCCC", coverage=20)
  #aligner = Aligner("ACTGTTGGAAAACCTTGTTGTACCCGGGTTAAACCCCC", targ)
  ref = reference.Reference("ACTGTTGGAAAACCTTGTTGTACCCGGGTTAAACCCCC")
  aligner = assembler.assemble("ACTGTTGGAAAACCTTGTTGTACCCGGGTTAAACCCCC",targ,math.ceil(0.7*20),0.9,False)
  #aligner.align()  # should contain 3 large contigs: ACTGTTGG, CCTTGTTGT, GGGTTAAACCCCC
  print "Here are your contigs:", aligner.ref.get_contigs(1) # TODO edit this thresh is higher
  #aligner.ref.plot_maxes(True) # Something wrong here even with coverage 20!

if __name__ == "__main__":
  test()