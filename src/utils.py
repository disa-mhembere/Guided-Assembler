#!/usr/bin/python

# utils.py
# Created on 2013-11-30.
# Email: disa@jhu.edu, slee320@jhu.edu
# Copyright (c) 2013. All rights reserved.

def Override(interface_class):
  """
  Method decorator for overriding class method name checking.
  Adapted from http://stackoverflow.com/questions/1167617/in-python-how-do-i-indicate-im-overriding-a-method
  """
  def overrider(method):
    assert method.__name__ in dir(interface_class), "No similar named method '%s' in super class '%s' ..." % (method.__name__, interface_class.__name__)
    return method
  return overrider

def edta(X,Y):
  """
  Compute the Edit distance between two sequences

  @param X: any arbritrary string
  @param Y: any arbritrary string
  """
  m = len(X)+1
  n = len(Y)+1

  D = [[[0,'']]*n for i in range(m)]

  #Initialize top row/column
  for i in range(m)[1:]:
    D[i][0]=[i,'I']
  for j in range(n)[1:]:
    D[0][j]=[j,'D']


  for i in range(m)[1:]:
    for j in range(n)[1:]:

      #Match character??
      s = 1
      if X[i-1] == Y[j-1]:
        s = 0

      d = min([D[i-1][j][0]+1,D[i][j-1][0]+1,D[i-1][j-1][0]+s])

      if d == D[i-1][j][0]+1:
        ch = 'I'
      if d == D[i][j-1][0]+1:
        ch = 'D'
      if d == D[i-1][j-1][0]+s:
        if s == 0:
          ch = 'M'
        else:
          ch = 'R'

      D[i][j] = [d,ch]

  print D[m-1][n-1][0]

  row = m-1
  col = n-1

  eX = ''
  eY = ''

  while row > 0 or col > 0:
    c = D[row][col][1]

    if c == 'I':
      eX = X[row-1] + eX
      eY = '-' + eY
      row = row - 1
    elif c == 'D':
      eX = '-' + eX
      eY = Y[col-1] + eY
      col = col - 1
    else:
      eX = X[row-1] + eX
      eY = Y[col-1] + eY
      row = row - 1
      col = col - 1

  print eX
  print eY

def pp(var):
  """
  Used in mapping operation for lists to ++ an index in the list
  @param var: some integer
  @return: plus one to value of var
  """
  return var + 1

if __name__ == "__main__":
  print "No main implementation for", __file__
