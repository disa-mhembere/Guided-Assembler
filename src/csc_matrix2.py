#!/usr/bin/python

# csc_matrix2.py
# Created by Disa Mhembere on 2013-11-21.
# Email: disa@jhu.edu
# Copyright (c) 2013. All rights reserved.

import scipy
from scipy.sparse.csc import csc_matrix
import numpy as np
from exceptions import IndexError

class csc_matrix2(csc_matrix):
  """
  Sub-class of lil_matrix that allows permits popping rows off
  """
  def pop_row(self, ):
    if self.shape[0] == 0:
      raise IndexError('Cannot pop a matrix with rows = 0')

    self.rows = np.delete(self.rows, self.shape[0]-1, 0) # clean up
    self.data = np.delete(self.data, self.shape[0]-1, 0) # clean up
    self._shape = (self._shape[0]-1, self.shape[1])

  def append_col(self, sp_mat=None, init=True):
    """
    Add a column to the sparse matrix given another sp_mat to append.
    Used when adding a new letter to the alphabet.

    @param sp_mat: the sparse matrix to be appended to the self object
    """
    self._shape = (self._shape[0], self.shape[1]+1)
    if sp_mat is not None:
      self[:,-1] = sp_mat

    elif init:
      self[self.shape[0]-1, self.shape[1]-1] = 1

  def append_row(self, ):
    """
    Append a row to the bottom of a lil_matrix2 object
    """
    self.rows = np.append(self.rows, 0)
    self.rows[-1] = []
    self.data = np.append(self.data, 0)
    self.data[-1] = []
    self._shape = (self._shape[0]+1, self.shape[1])

# Note we never remove columns -- too expensive