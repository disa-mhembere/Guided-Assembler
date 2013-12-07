#!/usr/bin/python

# target_reads.py
# Created on 2013-11-26.
# Email: disa@jhu.edu, slee320@jhu.edu
# Copyright (c) 2013. All rights reserved.

import numpy as np
import random
from copy import copy

class Target(object):
  def __init__(self, p, read_length, T, seed=None, coverage=5):
    """
    An Object that enables the partitioning &generation of synthetically mutated
    data given a full Target sequence read

    @param p: the probability of having a polymorphism at any single letter within the read
    @param read_length: read lengths that we will produce
    @param T: the target string which we are trying to assemble using R, the ref
    @param seed: seed for random
    @param coverage: max redundancy at any single nt position
    """
    assert (p <= 1 and p >= 0), "Probability must be in range [0:1] inclusive"
    assert isinstance(T, str), "Target (T) must be a string"
    assert not read_length > len(T), "Target (T) must be shorter than read length"

    self.T = T.upper() # All upper case
    self.p = p
    self.read_length = read_length

    self.coverage = coverage

    if seed: random.seed(seed) # for result reproducibility
    self.seen = [0]*len(T) # list of seen start indexes

  def get_read(self, **kwargs):
    """
    Get a single read from the Target sequence for a streaming read model

    @keyword p: the probability of SNP occuring
    @keyword read_length: the read length
    @keyword max_trials: the maximum number of times to rand redraw the idx if we keep finding over-covered idxs
    @return: a string with with a read obtained at a random position
    """
    if kwargs.has_key("p"): p = kwargs["p"]
    else: p = self.p

    if kwargs.has_key("read_length"): p = kwargs["read_length"]
    else: read_length = self.read_length

    assert set(kwargs.keys()).issubset(set(["p", "read_length", "max_trials"]))\
                                            ,"Unknown keyword argument in input"

    # pick idx in target string where to start from
    idx = np.random.random_integers(0, high=len(self.T)-read_length)  # The range is inclusive
    num_trials = 0 # Never let this loop go beyond 5 attempts to find a read
    while(self.seen[idx] > self.coverage and num_trials < 5): # don't over-cover
      idx = np.random.random_integers(0, high=len(self.T)-read_length) # The range is inclusive
      num_trials += 1

    read = self.T[idx:idx+read_length]

    #read = self.mutate(read) # A SNPs at random indexes

    #self._update_seen(idx)
    return read

  def mutate(self, read):
    """
    Use p to determine how to add SNPs to the returned string

    @param read: the read we want mutated
    @return: a string with some SNPs added
    """

    # Get indices to mutate; each has probability p of mutating
    mut_idx = []
    for i in xrange(len(read)):
      U = random.random()
      if U < self.p:
        mut_idx.append(i)

    read = bytearray(read)
    incoming = copy(read) # TODO: Testing

    # Mutate into something else
    comps = {"A":"CGT","C":"AGT","G":"ACT","T":"ACG"}
    for idx in mut_idx:
      read[idx] = random.choice(comps[chr(read[idx])])

    #if not read == incoming: print "Diff! %s != %s" % (incoming, read) # TODO: Testing
    return str(read)

  def _update_seen(self, idx):
    """
    Update the reference to seen start indicies so we dont continuously sample
    the same intervals. Private method.

    @param idx: The start index of the target string from where we just extracted a read
    """
    assert isinstance(idx, int), "Index must be an integer"

    self.seen[idx:idx+self.read_length] = map(pp, self.seen[idx:idx+self.read_length])

  def get_read_list(self, **kwargs):
    """
    Return a list of even-coverage randomly sampled strings with possible SNPs
    given T. For the non-streaming version of the read splitting.

    Optional:
    --------
    @keyword p: the probability of SNP occuring
    @keyword read_length: the read length
    @keyword save: boolean save or don't to disk
    @keyword save_fn: the filename you want to use to write to disk
    """

    if kwargs.has_key("p"): p = kwargs["p"]
    else: p = self.p

    if kwargs.has_key("read_length"): p = kwargs["read_length"]
    else: read_length = self.read_length

    assert set(kwargs.keys()).issubset(set(["p", "read_length", "save", "save_fn"]))\
                                            , "Unknown keyword argument in input"
    mutated_reads = list()

    for cov_idx in xrange(self.coverage+1):
      reads = [ self.T[l:l+read_length] for l in xrange(cov_idx, len(self.T)-read_length+1) ]
      reads = mutated_reads.extend(map(self.mutate, reads))

    if kwargs.get("save", False):
      import datetime as dt
      save_fn = kwargs.get("save_fn", "mutated_reads%s"(str(dt.datetime.now())))

      np.save(save_fn, mutated_reads)
    return mutated_reads

def pp(var):
  """
  Used in mapping operation for lists to ++ an index in the list
  @param var: some integer
  """
  return var + 1

def test():
  seq = "acgttttacccgggttac"
  T = Target(p=0.01, read_length=6, T=seq, seed=1234, coverage=3)

  #read_list = map(str,T.get_read_list())
  #print read_list

  for i in xrange(100):
    print T.get_read()

  print T.seen

def main():
  test()

if __name__ == "__main__":
  main()