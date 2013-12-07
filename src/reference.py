#!/usr/bin/python

# reference.py
# Created on 2013-11-26.
# Email: disa@jhu.edu, slee320@jhu.edu
# Copyright (c) 2013. All rights reserved.
import numpy as np

class Reference(object):
  def __init__(self, R, data_counts=9):
    """
    Object to hold a reference string and accompanying metadata associated with
    keeping track of how many matches/mismatches occur at each position

    @param R: Is the reference string
    @param data_counts: the counts related to which letter landed at a particular index
    """
    self.look_up = {"A":0,"C":1,"G":2,"T":3,"B":4,"D":5,"H":6,"U":7,"-":8}
    self.R = R # reference string
    self.dc = 9
    # keeps track of the letters that have landed at a particular index.
    self.match_count = np.zeros((len(R), self.dc))# 4 is for 'ACGTBDHU-' --STRICTLY in that order!

  def match(self, idx, char, cnt=1, thresh=10):
    """
    If there is a match in the reference at some position

    @param idx: the index in R where char matched
    @param char: the char that matched
    @param cnt: the number of times we should record char matched. Default=1
    @param thresh: threshold to indicate a change in the refernce should be made

    @return: True or False on (....something... TODO: SL)
    """

    if char in "BDHU-": thresh = 30

    assert idx >= 0 and idx < self.match_count.shape[0], "Out of bounds with index %d"%idx

    maxElement = self.match_count[idx,:].max()
    self.match_count[idx, self.look_up[char]] += cnt

    # If you exceed the threshold AND are a change, indicate
    # that a change to the reference should be made...
    if char in "BDHU-" and self.match_count[idx,self.look_up[char]] > thresh:
      self.match_count[idx,self.look_up[char]] = 0
      return True
    elif self.match_count[idx,self.look_up[char]] > max(maxElement,thresh)\
       and self.R[idx] != char:
      return True
    else:
      return False

  def build_hist(self, coverage, show=False, save=False, save_fn="max_hist_plot"):
    """
    Build a histogram to determine what the maxes look & visualize match_count
    Might be used to determine a resonable threshold

    @param coverage: the average coverage for an single nt
    @param show: Show visualization with match maxes
    @param save_fn: Save to disk with this file name or else it will be the default

    @return: the histogram array
    """
    #import matplotlib
    #matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    maxes = self.match_count.max(1) # get maxes along 1st dim

    h = plt.hist(maxes, bins=self.match_count.shape[0]) # figure out where the majority

    plt.ylabel("Frequency")
    plt.xlabel("Count per index")
    plt.title("Frequency count histogram")

    if show: plt.show()
    if save: plt.savefig(save_fn, dpi=160, frameon=False)

    return h[0]

  def get_contigs(self, thresh):
    """
    Use thresholding to determine what is a contig using the match_count

    @param thresh: anything under this value will not be included in the
    @return: the contigs found in the string
    """
    contigs = list()
    curr_contig = ""

    for idx in xrange(self.match_count.shape[0]):
      if not curr_contig: curr_contig_st_idx = idx # start position of contig in R

      # could replace with: `mx_idxs = np.where(self.match_count[idx] == self.match_count[idx].max())` # and figure out a tie breaker # TODO: DM,SL
      mx_idx = self.match_count[idx].argmax() # **NOTE: CAUTION - using argmax means the lexic 1st will win in the case of a tie!


      if self.match_count[idx, mx_idx] >= thresh and mx_idx < 4:
        curr_contig += "ACGT"[mx_idx] # Note assumption of ACGT here again
      else:
        if curr_contig: # If its not empty
          contigs.append((curr_contig, curr_contig_st_idx))
          curr_contig = ""

    # Add last contig if it exists
    if curr_contig: contigs.append((curr_contig, curr_contig_st_idx))

    return contigs

  def plot_maxes(self, show=False, save=False, save_fn="contig_plot"):
    """
    Plot maxes to try to visualize where contigs will lie

    @param show: boolean on if you want to show the image on the screen
    """
    #import matplotlib
    #matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    maxes = self.match_count.max(1) # get maxes along 1st dim

    plt.figure()
    plt.ylabel("Max Frequency")
    plt.title("Contig plot")
    plt.xlabel("Reference index")
    plt.plot(maxes, lw=3) # visualize where contigs may lie

    if show: plt.show()
    if save: plt.savefig(save_fn, dpi=160, frameon=False)

  def get_consensus(self, thresh):
    """
    Get the fraction of positions in the reference that have come to consus on the
    nt at that position.

    @param thresh: the number that defines what nt count is sufficient to be confident
    about the result at a single position
    @return: number that gives the fraction of positions that have come to consensus
    """
    return (np.where(self.match_count.max(1) >= thresh)[0].shape[0])/float(self.match_count.shape[0])

  def del_idx(self, idx):
    """
    Delete an index within the count array
    @param idx: the index we want deleted
    """
    assert idx >= 0 or idx < self.match_count.shape[0], "Out of bounds with index %d"%idx
    self.match_count = np.delete(self.match_count, idx, axis=0)

  def insert_idx(self, idx):
    """
    Insert an index.
    @param idx: the index where we want to insert the charactere
    """
    assert idx>=0 or idx < self.match_count.shape[0],"Out of bounds with index %d"%idx

    data = np.zeros((1, self.dc))

    if not idx == self.match_count.shape[0]: # last position insert gets no data
      self.match_count[idx, 4:8] = data[0,:4] # Give Prior BDHU to new ACGT

    self.match_count = np.insert(self.match_count, obj=idx, values=data ,axis=0)

  def zero_idx(self, idx):
    """
    Zero out an index within the count array
    @param idx: the index we want zeroed
    """
    assert idx >= 0 or idx < self.match_count.shape[0], "Out of bounds with index %d"%idx

    zeroRow = np.zeros((1, self.dc))
    self.match_count[idx,:] = zeroRow

def test(show=True):
  ref = Reference("ACGTTTTACCCGGGTTAC")

  # pretend read_length = 4, coverage = 5
  # should return: match CGTT @ idx = 1 and CCGG @ idx = 9

  ref.match(1, "C", 5); ref.match(2, "G", 5); ref.match(3, "C", 5); ref.match(4, "C", 5)
  ref.match(9, "C", 4);  ref.match(10, "C", 3); ref.match(11, "G", 4); ref.match(12, "G", 5);

  # bogus matches
  ref.match(1, "A", 1); ref.match(2, "T", 2); ref.match(3, "G", 1); ref.match(4, "A", 2)
  ref.match(9, "T", 2);  ref.match(10, "A", 1); ref.match(11, "A", 3); ref.match(12, "G", 5);

  contigs = ref.get_contigs(thresh=3)

  print contigs
  ref.del_idx(1)
  ref.insert_idx(12)

  print  ref.get_contigs(thresh=3)

  "Current concensus = %.3f %%" % (ref.get_consensus(thresh=3)*100.0)

  #ref.build_hist(coverage=5, show=show)

  #ref.plot_maxes(show=show)

  return ref, contigs

if __name__ == "__main__":
  test()
