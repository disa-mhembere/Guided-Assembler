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


def assemble(ref_seq, targ_seq):
  """
  Run the assembler on data
  """

  ref = parse_ref_data(fn)




def main():
  parser = argparse.ArgumentParser(description="Run the assembler and determine where contigs lie")
  parser.add_argument("ref", action="store", help="")
  parser.add_argument("read_length", action="store", type=int, help="How long each read is/should be")
  parser.add_argument("-t", "--test", action1="store_true", help="Run the mini-test to see if the world is not broken")
  parser.add_argument("-s", "--split_target", action="store_true", help="Run with a target that must be split and have SNPs added")

  parser.add_argument("-p", "--prob", action="store_true", default=0.1, help="prob of an SNP occuring in an target")
  parser.add_argument("-O", "--output", action="store", help="If we want to output written to disk instead of stdout")

  result = parser.parse_args()

  assert os.path.exists(result.ref), "File %s does not exits! Check the name and try again" % result.ref
  if result.test:
    print "Mini-test todo" # TODO: DM
    sys.exit(0) # should terminate after test


  ref_seq = open(result.ref).readline()


  if result.syn:
    targ = tr.Target(result.p, result.read_length, None) # TODO: DM FIXME
    pass # take the refer

  assemble(ref_seq, targ_seq)

if __name__ == "__main__":
  main()