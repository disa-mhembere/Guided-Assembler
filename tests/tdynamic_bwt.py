#!/usr/bin/python

# tdynamic_bwt.py
# Created by Disa Mhembere on 2013-12-06.
# Email: disa@jhu.edu
# Copyright (c) 2013. All rights reserved.

from src.dBWt import dBWT # TODO: Add src to path

def test_move_row():
  """ Test the functionality of dbwt moverow """
  print "Testing dBWT.moverow ..."
  f = dBWT("CGTAACGT")

  print "Before ..."
  print "F:", f.F
  print "L:", f.L

  # Test Moverow
  print "Move row 1 to 4..."
  f.moverow(f.F, f.L, 1, 4)

  print "After 1 ..."
  print "F:", f.F
  print "L:", f.L
  assert f.F == ["$", "A", "C", "C", "A", "G", "G", "T", "T"] \
      and f.L == ["T", "A", "A", "$", "T", "C", "C", "G", "G"], "Equiv Failure!"

  # Test move from begin to somewhere
  print "Test move row 0 to 6..."

  f.moverow(f.F, f.L, 0, 6)
  print "F:", f.F
  print "L:", f.L
  assert f.F == ["A", "C", "C", "A", "G", "G", "$", "T", "T"] \
      and f.L == ["A", "A", "$", "T", "C", "C", "T", "G", "G"], "Equiv Failure!"

  # Test move somewhere to end
  print "Move row 2 to 7"
  f.moverow(f.F, f.L, 2, 8)
  print "F:", f.F
  print "L:", f.L

  # Test move j > jp
  print "Move row 5 to 1"
  f.moverow(f.F, f.L, 5, 1)
  print "F:", f.F
  print "L:", f.L

  assert f.F == ["A", "$", "C", "A", "G", "G", "T", "T", "C"] \
      and f.L == ["A", "T", "A", "T", "C", "C", "G", "G", "$"], "Equiv Failure!"

  # Test move somewhere to end
  print "Move row 8 to 0"
  f.moverow(f.F, f.L, 8, 0)
  print "F:", f.F
  print "L:", f.L

  assert f.F == ["C", "A", "$", "C", "A", "G", "G", "T", "T"] \
    and f.L == ["$", "A", "T", "A", "T", "C", "C", "G", "G"], "Equiv Failure!"

  print "Test successfully completed!"