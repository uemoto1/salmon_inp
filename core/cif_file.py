#!/usr/bin/env python3
import sys
import re
import collections
import numpy as np
from . import tool
from numpy import pi, cos, sqrt, arccos
from numpy.linalg import norm

au_aa = 0.52917721067 # Bohr radius in Angstrom

class CIF:
    """Simple CIF Import/Export Class"""
    
    def __init__(self, sysname = '', a_prim=None, site_lbl=None, site_pos=None):
        """Create and Initialize CIF Object"""
        self.sysname = ''
        self.a_prim = a_prim
        self.site_pos = site_pos
        self.site_lbl = site_lbl


    def _parse_cif(self, cif_string):
        """Parse CIF Formatted File"""
        
        # tokenize the cif_string by delimiter ',",;
        buff = collections.deque()
        mark = ''
        for i, char in enumerate(cif_string):
            if not mark:
                if char in ['"', "'", ';']:
                    start, mark = i + 1, char
                elif not char.isspace():
                    start, mark = i, '*'    
                        
            elif ((mark == '*' and char.isspace()) or
                  (mark in ["'", '"', ';'] and char == mark) or
                  (mark in ["'", '"'] and char in ['\n', '\r'])):
                mark = ''
                buff += [cif_string[start:i]]
                
        if mark:
            buff += [cif_string[start:i]]
        
        # parse data structure in the CIF file
        state = 0
        data = {}
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
        

    def loads(self, cif_string):
        
        def lbl2elem(s):
            return re.match('([A-Z][a-z]*).*', s).group(1)
        
        def chem2sysname(s):
            return re.sub(r'\W', '', s)
        
        def del_e(s):
            return re.sub(r"\(\d+\)", "", s)
            
        self.data = self._parse_cif(cif_string)
        
        if '_chemical_name_systematic' in self.data:
            self.sysname = chem2sysname(self.data['_chemical_name_systematic'])
        elif '_chemical_name_common' in self.data:
            self.sysname = chem2sysname(self.data['_chemical_name_common'])
        elif '_chemical_name_mineral' in self.data:
            self.sysname = chem2sysname(self.data['_chemical_name_mineral'])
        # remove error range
        self.data['_atom_site_fract_x'] = [
            del_e(it) for it in self.data['_atom_site_fract_x']
        ]
        self.data['_atom_site_fract_y'] = [
            del_e(it) for it in self.data['_atom_site_fract_y']
        ]
        self.data['_atom_site_fract_z'] = [
            del_e(it) for it in self.data['_atom_site_fract_z']
        ]
        # construction of unit cell
        angle_ab = float(del_e(self.data.get('_cell_angle_gamma', 90)))
        angle_bc = float(del_e(self.data.get('_cell_angle_alpha', 90)))
        angle_ca = float(del_e(self.data.get('_cell_angle_beta', 90)))
        cos_ab = cos(pi / 180 * angle_ab)
        cos_bc = cos(pi / 180 * angle_bc)
        cos_ca = cos(pi / 180 * angle_ca)
        # unit vectors for lattice axis
        u_prim = np.empty([3, 3])
        # a
        u_prim[0] = [1.0, 0.0, 0.0]
        # b
        u_prim[1, 0] = cos_ab
        u_prim[1, 1] = sqrt(1.0 - u_prim[1, 0] ** 2)
        u_prim[1, 2] = 0.0
        # c
        u_prim[2, 0] = cos_ca
        u_prim[2, 1] = (cos_bc - u_prim[1, 0] * u_prim[2, 0]) / u_prim[1, 1]
        u_prim[2, 2] = sqrt(1.0 - u_prim[2, 0] ** 2 - u_prim[2, 1] ** 2)
        # lattice length
        length_a = float(del_e(self.data['_cell_length_a']))  / au_aa
        length_b = float(del_e(self.data['_cell_length_b']))  / au_aa
        length_c = float(del_e(self.data['_cell_length_c']))  / au_aa
        # latticce vector
        self.a_prim = np.empty([3, 3])
        self.a_prim[0] = u_prim[0] * length_a
        self.a_prim[1] = u_prim[1] * length_b
        self.a_prim[2] = u_prim[2] * length_c
        # atomic coordinates
        self.site_pos_orig = np.vstack([
            np.array(self.data['_atom_site_fract_x'], dtype=float),
            np.array(self.data['_atom_site_fract_y'], dtype=float),
            np.array(self.data['_atom_site_fract_z'], dtype=float),
        ]).T
        # atomic labels
        self.site_lbl_orig = [
            lbl2elem(it) for it in self.data['_atom_site_label']
        ]
        # symmetry operations
        if '_symmetry_equiv_pos_as_xyz' in self.data:
            self.symop_list = self.data['_symmetry_equiv_pos_as_xyz']
        elif '_space_group_symop_operation_xyz' in self.data:
            self.symop_list = self.data['_space_group_symop_operation_xyz']
        else:
            self.symop_list = []
        # calculate symmetry equivalent position: r_sym
        self.site_pos, self.site_lbl = [], []
        for r, lbl in zip(self.site_pos_orig, self.site_lbl_orig):
            for symop in self.symop_list:
                symop2 = re.sub(r'(\d+)/', r'\1.0/', symop.lower())
                r_sym = np.array(eval(symop2, {'x': r[0], 'y': r[1], 'z': r[2]}))
                if tool.is_individual_pos(self.a_prim, r_sym, self.site_pos):
                    self.site_pos += [r_sym % 1]
                    self.site_lbl += [lbl2elem(lbl)]        
        return
                
    
    def dumps(self):
        buff = [
            "# This file is genereated by salmon_inp tool",
            "data_%s" % self.sysname,
            "_chemical_name_mineral %s" % self.sysname,
            "_space_group_IT_number 1",
            "_symmetry_int_tables_number 1",
            "_symmetry_space_group_name_Hall 'P 1'",
            "_symmetry_space_group_name_H-M 'P 1'",
            "_symmetry_cell_setting triclinic",
            "_cell_length_a %.2f" % (norm(self.a_prim[0]) * au_aa),
            "_cell_length_b %.2f" % (norm(self.a_prim[1]) * au_aa),
            "_cell_length_c %.2f" % (norm(self.a_prim[2]) * au_aa),
            "_cell_angle_alpha %d" % (
                arccos(tool.dot_cos(self.a_prim[1], self.a_prim[2])) * 180 / pi
            ),
            "_cell_angle_beta %d" % (
                arccos(tool.dot_cos(self.a_prim[2], self.a_prim[0])) * 180 / pi
            ),
            "_cell_angle_gamma %d" % (
                arccos(tool.dot_cos(self.a_prim[0], self.a_prim[1])) * 180 / pi
            ),
            "loop_",
            "_symmetry_equiv_pos_as_xyz",
            "x,y,z",
            "loop_",
            "_atom_site_label",
            "_atom_site_fract_x",
            "_atom_site_fract_y",
            "_atom_site_fract_z",
        ]
        
        count = collections.defaultdict(int)
        for lbl, pos in zip( self.site_lbl, self.site_pos):
            count[lbl] += 1
            buff += ["%s%d %.4f %.4f %.4f" % (
                lbl, count[lbl], pos[0], pos[1], pos[2]
            )]
        
        return "\n".join(buff)
        
