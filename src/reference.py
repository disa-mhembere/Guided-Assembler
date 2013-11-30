#!/usr/bin/python

# reference.py
# Created on 2013-11-26.
# Email: disa@jhu.edu, slee320@jhu.edu
# Copyright (c) 2013. All rights reserved.
import numpy as np

class Reference(object):
  def __init__(self, len_R):
    """
    Object to hold a reference string and accompanying metadata associated with
    keeping track of how many matches/mismatches occur at each position

    @param len_R: Is the length of the reference string. Note the actual
    string is obtainable via the BWT so no need to hold it here.
    """
    self.look_up = {"A":0, "C":1, "G":2, "T":3}

    # keeps track of the letters that have landed at a particular index.
    self.match_count = np.zeros((len_R, 4))# 4 is for 'ACGT' --STRICTLY in that order!

  def match(self, idx, char, cnt=1):
    """
    If there is a match in the reference at some position

    @param idx: the index in R where char matched
    @param char: the char that matched
    @param cnt: the number of times we should record char matched. Default=1
    """
    self.match_count[idx, self.look_up[char]] += cnt

  def build_hist(self, coverage, show=False, save=False, save_fn="max_hist_plot"):
    """
    Build a histogram to determine what the maxes look & visualize match_count
    Might be used to determine a resonable threshold

    @param coverage: the average coverage for an single nt
    @param show: Show visualization with match maxes
    @param save_fn: Save to disk with this file name or else it will be the default

    @return: the histogram array
    """
    import matplotlib.pyplot as plt

    maxes = self.match_count.max(1) # get maxes along 1st dim

    h = plt.hist(maxes, bins=self.match_count.shape[0]) # figure out where the majority
    if show: plt.show()
    if save: plt.savefig(save_fn, dpi=160, frameon=False)

    return h[0]

  def get_contigs(self, thresh):
    """
    Use thresholding to determine what is a contig using the match_count

    @param thresh: anything under this value will not be included in the
    """
    contigs = list()
    curr_contig = ""

    for idx in xrange(self.match_count.shape[0]):
      if not curr_contig: curr_contig_st_idx = idx # start position of contig in R

      # could replace with: `mx_idxs = np.where(self.match_count[idx] == self.match_count[idx].max())` # and figure out a tie breaker # TODO: DM,SL
      mx_idx = self.match_count[idx].argmax() # **NOTE: CAUTION - using argmax means the lexic 1st will win in the case of a tie!

      if self.match_count[idx, mx_idx] >= thresh:
        curr_contig += "ACGT"[mx_idx] # Note assumption of ACGT here again
      else:
        if curr_contig: # If its not empty
          contigs.append((curr_contig, curr_contig_st_idx))
          curr_contig = ""

    # Add last contig if it exists
    if curr_contig: contigs.append(curr_contig, curr_contig_st_idx)

    return contigs

  def plot_maxes(self, show=False):
    """
    Plot maxes to try to visualize where contigs will lie
    """
    import matplotlib.pyplot as plt
    maxes = self.match_count.max(1) # get maxes along 1st dim

    plt.figure()
    plt.plot(maxes) # visualize where contigs may lie

    if show: plt.show()

def test():
  seq = "ACGTTTTACCCGGGTTAC"
  ref = Reference(len(seq))

  # pretend read_length = 4, coverage = 5
  # should return: match CGTT @ idx = 1 and CCGG @ idx = 9

  ref.match(1, "C", 5); ref.match(2, "G", 5); ref.match(3, "C", 5); ref.match(4, "C", 5)
  ref.match(9, "C", 4);  ref.match(10, "C", 3); ref.match(11, "G", 4); ref.match(12, "G", 5);

  # bogus matches
  ref.match(1, "A", 1); ref.match(2, "T", 2); ref.match(3, "G", 1); ref.match(4, "A", 2)
  ref.match(9, "T", 2);  ref.match(10, "A", 1); ref.match(11, "A", 3); ref.match(12, "G", 5);

  print ref.get_contigs(thresh=3)

  print ref.match_count
  ref.build_hist(coverage=5, show=True)

  ref.plot_maxes(show=True)


if __name__ == "__main__":
  test()