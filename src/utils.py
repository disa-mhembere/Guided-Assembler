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

if __name__ == "__main__":
  print "No main implementation for", __file__