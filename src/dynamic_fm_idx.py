#!/usr/bin/python

# dynamic_fm_idx.py
# Created by Disa Mhembere on 2013-11-18.
# Email: disa@jhu.edu
# Copyright (c) 2013. All rights reserved.

from dynamic_bwt import dBWT

class dFMidx(dBWT):

  def __init__(self, isa=False):
    super(dBWT, self).__init__(isa)

    #TODO