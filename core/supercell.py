#!/usr/bin/env python3

from itertools import product
from numpy import dot, array, empty, argmax, rint
from numpy.linalg import norm, inv, det
from . import tool

eps = 1e-3

def search_supercell(a_prim, rion, kion, nmax=2):

    def func_score(it):
        ns, a_temp = it
        p0 = abs(tool.dot_cos(a_temp[0], a_prim[0]))
        p1 = abs(tool.dot_cos(a_temp[1], a_prim[1]))
        p2 = abs(tool.dot_cos(a_temp[2], a_prim[2]))
        return (p2, p1, p0)

    nlist = [
        array([i0, i1, i2])
        for (i0, i1, i2) in product(
            range(-nmax, nmax + 1),
            range(-nmax, nmax + 1),
            range(-nmax, nmax + 1)
        )
        if (i0, i1, i2) != (0, 0, 0)
    ]

    a_temp = empty([3, 3])  # Synthesized lattice vector

    v_prim = abs(det(a_prim)) # volume of primitive cell
    v_min = None  # Minimum volume of lattice cell (a_temp)

    result = []  # Available set of (n0, n1, n2)

    for n0 in nlist:
        a_temp[0] = dot(n0, a_prim)

        for n1 in nlist:
            a_temp[1] = dot(n1, a_prim)
            if abs(tool.dot_cos(a_temp[0], a_temp[1])) > eps:
                continue

            for n2 in nlist:
                a_temp[2] = dot(n2, a_prim)
                if abs(tool.dot_cos(a_temp[0], a_temp[2])) > eps:
                    continue
                if abs(tool.dot_cos(a_temp[1], a_temp[2])) > eps:
                    continue

                v = rint(det(a_temp) / v_prim) # Volume of the cell
                
                if 0 < v:
                    if (v_min is None) or (v < v_min):
                        v_min = v
                        result = []
                    if v == v_min:
                        result += [([n0, n1, n2], a_temp.copy())]

    priority_score = [
        func_score(ns) for ns in result
    ]
    
    ns, a_orth = result[argmax(priority_score)]

    M = inv(ns).T

    tlist = [
        array([t0, t1, t2])
        for (t0, t1, t2) in product(
            range(-nmax, nmax + 1),
            range(-nmax, nmax + 1),
            range(-nmax, nmax + 1),
        )
    ]

    new_rion = []
    new_kion = []

    for r, k in zip(rion, kion):

        for t in tlist:
            R = dot(M, r + t)
            if (all(0 <= R) and all(R < 1)):
                new_rion += [R[:]]
                new_kion += [k]

    return a_orth, new_rion, new_kion
