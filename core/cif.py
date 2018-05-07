#!/usr/bin/env python3
import re
from . import tool
import collections
from numpy import pi, cos, sqrt, empty, array, vstack

au_aa = 0.52917721067 # Bohr radius in Angstrom

def split_cif(text):
    mark = ''
    for i, char in enumerate(text):
        if not mark:
            if char in ['"', "'", ';']:
                start, mark = i + 1, char
            elif not char.isspace():
                start, mark = i, "*"
                
        elif ((mark == "*" and char.isspace()) or
              (mark in ["'", '"', ';'] and char == mark) or
              (mark in ["'", '"'] and char in ['\n', '\r'])):
            mark = ""
            yield text[start:i]
            
    if mark:
        yield text[start:]


def parse_cif(text):
    buff = collections.deque(split_cif(text))
    data = {}
    state = 0

    while buff:
        item = buff.popleft()

        if state == 0:
            if re.match(r"_\w+", item):
                data[item] = buff.popleft()
            elif re.match(r"(loop|LOOP)_", item):
                column = []
                state = 1

        elif state == 1:
            if re.match(r"_\w+", item):
                column += [item]
                data[item] = []
            else:
                buff.appendleft(item)
                state = 2

        elif state == 2:
            if not re.match(r"(loop|LOOP)?_\w*", item):
                buff.appendleft(item)
                for tag in column:
                    data[tag] += [buff.popleft()]
            else:
                buff.appendleft(item)
                state = 0
    
    return data
    
def label2symbol(label):
    res = re.match("([A-Z][a-z]*).*", label)
    return res.group(1)
    
def read_cif(text):
    data = parse_cif(text)
    
    # atomic coordinates
    site_pos = vstack([
        array(data["_atom_site_fract_x"], dtype=float),
        array(data["_atom_site_fract_y"], dtype=float),
        array(data["_atom_site_fract_z"], dtype=float),
    ]).T
    
    # atomic labels
    site_label = data["_atom_site_label"]

    # construction of unit cell
    cos_ab = cos(pi / 180 * float(data.get("_cell_angle_gamma", 90)))
    cos_bc = cos(pi / 180 * float(data.get("_cell_angle_alpha", 90)))
    cos_ca = cos(pi / 180 * float(data.get("_cell_angle_beta", 90)))
    # unit vectors for lattice axis
    u_prim = empty([3, 3])
    # a
    u_prim[0] = [1, 0, 0]
    # b
    u_prim[1, 0] = cos_ab
    u_prim[1, 1] = sqrt(1.0 - u_prim[1, 0] ** 2)
    u_prim[1, 2] = 0
    # c
    u_prim[2, 0] = cos_ca
    u_prim[2, 1] = (cos_bc - u_prim[1, 0] * u_prim[2, 0]) / u_prim[1, 1]
    u_prim[2, 2] = sqrt(1.0 - u_prim[2, 0] ** 2 - u_prim[2, 1] ** 2)
    # latticce vector
    a_prim = empty([3, 3])
    a_prim[0] = u_prim[0] * (float(data["_cell_length_a"][0]) / au_aa)
    a_prim[1] = u_prim[1] * (float(data["_cell_length_b"][0]) / au_aa)
    a_prim[2] = u_prim[2] * (float(data["_cell_length_c"][0]) / au_aa)

    # symmetry operations
    if "_symmetry_equiv_pos_as_xyz" in data:
        symop_list = data["_symmetry_equiv_pos_as_xyz"]
    elif "_space_group_symop_operation_xyz" in data:
        symop_list = data["_space_group_symop_operation_xyz"]
    else:
        symop_list = []

    site_pos2, site_label2 = [], []
    for r, lbl in zip(site_pos, site_label):
        for symop in symop_list:
            # calculate symmetry equivalent position: r_sym
            symop2 = re.sub(r"(\d+)/", r"\1.0/", symop.lower())
            r_sym = array(eval(symop2, {"x": r[0], "y": r[1], "z": r[2]}))
            #
            if not tool.is_same_position(a_prim, r_sym, site_pos2):
                site_pos2 += [r_sym % 1]
                site_label2 += [label2symbol(lbl)]
                
    return data["_chemical_name_mineral"], a_prim, site_pos2, site_label2
