#!/usr/bin/python

# assembler.py
# Created on 2013-11-26.
# Email: disa@jhu.edu, slee320@jhu.edu
# Copyright (c) 2013. All rights reserved.

import argparse
import sys
import os

import target_reads as tr
import reference as rf
import aligner as al
from math import ceil
import utils
import random

def assemble(reference, target, threshold, min_consensus):
  """
  Run the assembler on data


  @param reference:
  @param target:
  @param threshold:
  @param min_consensus:
  """
  aligner = al.Aligner(reference, target)


  failCount = 0
  failThresh = 200
  while(True):
    if failCount > failThresh:
      aligner.tol += 1
      failCount = 0
    indicators,transcript = aligner.align() # STUB
    if indicators == -1:
      failCount+=1
      continue
    if aligner.ref.get_consensus(threshold) >= min_consensus: break
    else:
      if True in indicators: aligner.alter_bwt(indicators,transcript) # STUB

  print "Sequence assembly complete!"
  print aligner.ref.match_count
  return aligner

def eval_acc(target_seq, contigs):
  if not contigs:
    print "No contigs found!"
    return

  recon_target = ""
  recon_target += " "*contigs[0][1]

  last = 0
  matches = 0
  for contig in contigs:
    contig_idx = contig[1]
    recon_target += " "*(contig_idx-last)
    recon_target += contig[0]
    last = len(recon_target)

  recon_target = recon_target.strip()

    # # Eval accuracy
    # for i, c in enumerate(recon_target[contig_idx:]):
    #   if c == target_seq[i+contig_idx]:
    #     matches += 1

  utils.edta(target_seq,recon_target)

  #print "Reconstrution complete with %.3f%% accuracy ..." % ((matches/float(len(target_seq)))*100)
  print "Target = %s" % target_seq
  print "Rarget = %s" % recon_target

def randStrings(n,corrupt):
  """
  Generate a random reference string and a related target string

  @param n: length of toy strings
  @param corrupt: Percentage of corrupted nt
  """

  R = ''
  for i in xrange(n):
    R = R + random.choice("ACGT")

  T = R

  corr = ceil(n*corrupt)
  ml = n/2 - int(corr*0.35)
  # Corrupted tides appear first, last, and middle
  indices = range(0,int(0.15*corr),1)+range(ml,ml+int(0.7*corr),1)+range(n-int(0.15*corr),n,1)
  T = bytearray(R)

  # Mutate into something else
  for idx in indices:
    T[idx] = random.choice("ACGT")

  return R,str(T)


def main():
  parser = argparse.ArgumentParser(description="Run the assembler and determine where contigs lie using a dynamic bwt index")
  parser.add_argument("ref", action="store", help="The file path to the reference sequence")
  parser.add_argument("targ", action="store", help="The file path to the target sequence")

  ## Optional ##
  parser.add_argument("-r", "--read_length", action="store", type=int, help="How long each read is/should be")
  parser.add_argument("-t", "--test", action="store_true", help="Run the mini-test to see if the world is not broken")
  parser.add_argument("-T", "--threshold", action="store", type=float, default=0.7, help="Fraction used to consider a value as part of \
                          a contig. i.e if it is 0.6 then must have ceil(0.6*coverage) matches to add the index to a contig sequence. Default=0.7")
  parser.add_argument("-s", "--split_target", action="store_true", help="Run with a target that must be split and have SNPs added")
  parser.add_argument("-c", "--coverage", action="store", type=int, default=10, help="The ideal average coverage. Default=10")
  parser.add_argument("-p", "--prob", action="store", type=float, default=0.05, help="prob of an SNP occuring in an target. Default=0.05")
  parser.add_argument("-O", "--output_filename", action="store", help="If we want output written to disk instead of stdout") # TODO DM
  parser.add_argument("-m", "--min_consensus", action="store", type=float, default=0.75, help="Minimum fraction of consensus for all\
                          positions when assembly is performed. Default=0.75")
  parser.add_argument("-n", "--test_length",action="store",type=int,help="How long the test strings should be")
  parser.add_argument("-C", "--corruption", action="store",type=float,help="How much mutation in the test string")

  parser.add_argument("-e", "--eval_acc", action="store_true", help="Given the correct target result evaluate the accuracy of the assembly")
  parser.add_argument("-P", "--plot", action="store_true", help="Display ALL plots visually. *Note: Causes os.system('pause') until figure is closed")  # TODO DM
  parser.add_argument("-S", "--save_figs", action="store_true", help="Save figures to disk? Boolean")  # TODO DM
  parser.add_argument("-F", "--figure_name", action="store", help="Save figures to disk with this file name. If this is set the -S is unnessary")# TODO DM


  result = parser.parse_args()


  if result.test:
    # target, contigs = rf.test(False)
    # eval_acc(target.R, contigs)
    # sys.exit(0) # should terminate after test

    ref_seq,targ_seq = randStrings(result.test_length,result.corruption)
    print "REF:",ref_seq
    print "TAR:",targ_seq
    if result.split_target:
      targ = tr.Target(result.prob, result.read_length, targ_seq, result.coverage)
    else:
      assert False, "Target must be split at this point!"

    aligner = assemble(rf.Reference(ref_seq).R, targ, ceil(result.threshold*result.coverage), result.min_consensus)

    print "REF:",aligner.ref.R
    print "TAR:",targ_seq

    aligner.ref.build_hist(result.coverage, True)
    aligner.ref.plot_maxes(True)

    if result.eval_acc:
      print "Test?"
      eval_acc(targ_seq, aligner.ref.get_contigs(result.threshold))
    exit()




  assert os.path.exists(result.ref), "File %s does not exits! Check the name and try again" % result.ref
  assert os.path.exists(result.targ), "File %s does not exits! Check the name and try again" % result.targ

  ref_seq = open(result.ref).readline().strip().upper()
  targ_seq = open(result.targ).readline().strip().upper()


  print "REF:",ref_seq
  print "TAR:",targ_seq

  if result.split_target:
    targ = tr.Target(result.prob, result.read_length, targ_seq, result.coverage)
  else:
    assert False, "Target must be split at this point!"

  aligner = assemble(rf.Reference(ref_seq).R, targ, ceil(result.threshold*result.coverage), result.min_consensus)

  print "REF:",aligner.ref.R
  print "TAR:",targ_seq

  if result.eval_acc:
    eval_acc(targ_seq, aligner.ref.get_contigs(result.threshold))


if __name__ == "__main__":
  main()