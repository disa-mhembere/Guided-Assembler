#!/usr/bin/python

# target_reads.py
# Created on 2013-11-26.
# Email: disa@jhu.edu, slee320@jhu.edu
# Copyright (c) 2013. All rights reserved.

import numpy as np

class Target(object):
  def __init__(self, p, read_length, T):
    """
    An Object that enables the partitioning &generation of synthetically mutated
    data given a full Target sequence read

    @param p: the probability of a polymorphism in any single letter with the read
    @param read_length: read lengths that we will produce
    @param T: the target string which we are trying to assemble using R, the ref
    """
    assert (p <= 1 and p >= 0), "Probability must be in range [0:1] inclusive"
    assert isinstance(T, str), "Target (T) must be a string"
    assert not read_length > len(T), "Target (T) must be shorter than read length"

    self.T = T.upper() # All upper case
    self.p = p
    self.read_length = read_length

    self.seen = list() # list of seen start indexes

  def get_read(self, **kwargs):
    """
    Get a single read from the Target read

    @todo: SL

    @keyword p: the probability of SNP occuring
    @keyword read_length: the read length
    @return: a string with with a read obtained at a random position
    """
    if kwargs.has_key("p"): p = kwargs["p"]
    else: p = self.p

    if kwargs.has_key("read_length"): p = kwargs["read_length"]
    else: read_length = self.read_length

    assert set(kwargs.keys()).issubset(set(["p", "read_length"]))\
                                            , "Unknown keyword argument in input"

    # pick idx in target string where to start from
    idx = np.random.random_integers(0, high=len(self.T)-read_length) # The range is inclusive
    read = self.T[idx:read_length]

    read = self.mutate(read, p) # TODO: SL

    self._update_seen(idx)
    return read

  def mutate(self, read, p):
    """
    Use p to determine how to add SNPs to the returned string

    @todo: SL

    @param: p the probability of SNP occurring # TODO SL verify
    @return: a string with some possible SNPs added
    """
    raise NotImplementedError("***TODO: SL****")


  def _update_seen(self, idx):
    """
    Update the reference to seen start indicies so we dont continuously sample
    the same intervals. Private method.

    """
    assert isinstance(idx, int), "Index must be an integer"
    self.seen.append(idx)

    self.seen

  def get_read_list(self, **kwargs):
    """
    Return a list of randomly sampled strings with possible SNPs given T
    """

    if kwargs.has_key("p"): p = kwargs["p"]
    else: p = self.p

    if kwargs.has_key("read_length"): p = kwargs["read_length"]
    else: read_length = self.read_length

    assert set(kwargs.keys()).issubset(set(["p", "read_length"]))\
                                            , "Unknown keyword argument in input"

    raise NotImplementedError("TODO: SL ")

def test():
  T = Target(0.1, 6, "acgttttacccgggttac")
  T._update_seen(4)
  print "Seen so far = ", T.seen
  try:
    print T.get_read()
    print get_ris
  except Exception, msg:
    print msg

def main():
  test()


if __name__ == "__main__":
  main()