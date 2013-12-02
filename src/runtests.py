#!/usr/bin/python

# runtests.py
# Created by Disa Mhembere on 2013-12-02.
# Email: disa@jhu.edu
# Copyright (c) 2013. All rights reserved.

from subprocess import call
import sys

def runtests(num_tests, scriptname, out_fn):

  csv = "Test #, Corruption, %% Error, Total Edit distance\n"
  f = open(out_fn, "ab")
  f.write(csv)
  f.close()

  C = [0.6]
  for c in C:
    for i in xrange(num_tests):
      print "Running "
      call(["python", scriptname, "-s", "-p 0.01", "-r 10", "-e", "-c 20", "-m 0.9", "-t", "-n 60", "-C %d"%c, ">>", out_fn])


  f = open(out_fn, "ab")
  f.write(csv)
  f.close()


if __name__ == "__main__":
  assert len(sys.argv) > 1, "You must specify what file to write to .."
  runtests(10, "assembler.py", sys.argv[1])