#!/usr/bin/env python3
import re
from numpy import array, dot
from numpy.linalg import norm
from itertools import product

def periodic_distance(a_prim, r0, r1, ishift=[-1, 0, 1]):
    return min(
        norm(dot(a_prim, r1 - r0 + array([i0, i1, i2])))
        for [i0, i1, i2] in product(ishift, ishift, ishift)
    )

def is_same_position(a_prim, r, r_list, eps=1e-1):
    for ri in r_list:
        if periodic_distance(a_prim, r, ri) < eps:
            return True
    return False
            
def label2symbol(label):
    res = re.match("([A-Z][a-z]*).*", label)
    return res.group(1)
