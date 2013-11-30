#!/usr/bin/python

# assembler.py
# Created on 2013-11-26.
# Email: disa@jhu.edu, slee320@jhu.edu
# Copyright (c) 2013. All rights reserved.

import pdb
import argparse
import sys
import os

import target_reads as tr
import reference as rf
import aligner as al
from math import ceil

def assemble(reference, target, threshold, min_consensus):
  """
  Run the assembler on data


  @param reference:
  @param target:
  @param threshold:
  @param min_consensus:
  """
  aligner = al.Aligner(reference, target)

  while(True):
    aligner.align() # STUB
    if aligner.ref.get_consensus() >= min_consensus: break
    else: aligner.alterbwt() # STUB

  print "Sequence assembly complete!"
  return aligner

def eval_acc(target_seq, contigs):
  reconstructed_target = ""

  pdb.set_trace()
  idxs, contig = zip(*contigs)


  # TODO

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
  parser.add_argument("-O", "--output", action="store", help="If we want output written to disk instead of stdout") # TODO DM
  parser.add_argument("-m", "--min_consensus", action="store", type=float, default=0.75, help="Minimum fraction of consensus for all\
                          positions when assembly is performed. Default=0.75")

  parser.add_argument("-e", "--eval_acc", action="store_true", help="Given the correct target result evaluate the accuracy of the assembly")
  parser.add_argument("-P", "--plot", action="store_true", help="Display ALL plots visually. *Note: Causes os.system('pause') until figure is closed")  # TODO DM
  parser.add_argument("-S", "--save_figs", action="store_true", help="Save figures to disk? Boolean")  # TODO DM
  parser.add_argument("-F", "--figure_name", action="store", help="Save figures to disk with this file name. If this is set the -S is unnessary")# TODO DM


  result = parser.parse_args()

  if result.test:
    target, contigs = rf.test(False)
    eval_acc(target, contigs)

    print "Mini-test todo" # TODO: DM IMPLEMENT-ME
    sys.exit(0) # should terminate after test

  assert os.path.exists(result.ref), "File %s does not exits! Check the name and try again" % result.ref
  assert os.path.exists(result.targ), "File %s does not exits! Check the name and try again" % result.targ

  ref_seq = open(result.ref).readline().strip().upper()
  targ_seq = open(result.targ).readline().strip().upper()

  if result.split_target:
    targ = tr.Target(result.p, result.read_length, ref_seq, result.coverage)
  else:
    assert False, "Target must be split at this point!"

  aligner = assemble(rf.Reference(ref_seq), targ, ceil(result.threshold*result.coverage), result.min_consensus)

  if result.eval_acc:
    eval_acc(targ_seq, aligner.get_contigs)


if __name__ == "__main__":
  main()