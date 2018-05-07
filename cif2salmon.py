#!/usr/bin/env python3
from core import cif
from core import orthonormalize
from core import tool
from numpy.linalg import norm

tbl_fhi_abinit = {
    # 'Symbol': [Z, nelec, lloc]
    'H': [1, 1, 2],'He': [2, 2, 2],'Li': [3, 1, 2],'Be': [4, 2, 2],'B': [5, 3, 2],
    'C': [6, 4, 2],'N': [7, 5, 2],'O': [8, 6, 2],'F': [9, 7, 2],'Ne': [10, 8, 2],
    'Na': [11, 1, 2],'Mg': [12, 2, 2],'Al': [13, 3, 2],'Si': [14, 4, 2],'P': [15, 5, 2],
    'S': [16, 6, 2],'Cl': [17, 7, 2],'Ar': [18, 8, 2],'K': [19, 1, 0],'Ca': [20, 2, 0],
    'Sc': [21, 3, 1],'Ti': [22, 4, 1],'V': [23, 5, 1],'Cr': [24, 6, 1],'Mn': [25, 7, 1],
    'Fe': [26, 8, 1],'Co': [27, 9, 1],'Ni': [28, 10, 0],'Cu': [29, 11, 0],'Zn': [30, 12, 0],
    'Ga': [31, 3, 0],'Ge': [32, 4, 0],'As': [33, 5, 2],'Se': [34, 6, 2],'Br': [35, 7, 2],
    'Kr': [36, 8, 2],'Rb': [37, 1, 0],'Sr': [38, 2, 0],'Y': [39, 3, 1],'Zr': [40, 4, 1],
    'Nb': [41, 5, 1],'Mo': [42, 6, 1],'Tc': [43, 7, 1],'Ru': [44, 8, 1],'Rh': [45, 9, 1],
    'Pd': [46, 10, 1],'Ag': [47, 11, 1],'Cd': [48, 12, 1],'In': [49, 3, 0],'Sn': [50, 4, 0],
    'Sb': [51, 5, 0],'Te': [52, 6, 0],'I': [53, 7, 0],'Xe': [54, 8, 0],'Cs': [55, 1, 1],
    'Ba': [56, 2, 0],'Ce': [58, 4, 0],'Nd': [60, 6, 0],'Pm': [61, 7, 0],'Sm': [62, 8, 0],
    'Gd': [64, 10, 0],'Er': [68, 14, 0],'Tm': [69, 15, 0],'Yb': [70, 16, 0],'Lu': [71, 17, 0],
    'Hf': [72, 4, 1],'Ta': [73, 5, 1],'W': [74, 6, 3],'Re': [75, 7, 1],'Os': [76, 8, 1],
    'Ir': [77, 9, 1],'Pt': [78, 10, 1],'Au': [79, 11, 1],'Hg': [80, 12, 1],'Tl': [81, 3, 0],
    'Pb': [82, 14, 0],'Bi': [83, 5, 0],'Po': [84, 6, 0],'At': [85, 7, 0],'Rn': [86, 8, 0],
}

sysname, a_prim, site_pos, site_label = cif.read_cif(open("C.cif").read())

ns, a_orth = orthonormalize.search_orthonormal_lattice(a_prim)

site_pos, site_label = orthonormalize.generate_atomic_position(
    ns, site_pos, site_label
)

atom_red_coord = []
izatom = []
pseudo_file = []
lloc_ps = []

label_list = list(set(site_label))
for i, lbl in enumerate(label_list):
    izatom += [str(tbl_fhi_abinit[lbl][0])]
    lloc_ps += [str(tbl_fhi_abinit[lbl][2])]
    pseudo_file += ['%s.fhi' % lbl]
    
total_nelec = 0
for r, lbl in zip(site_pos, site_label):
    i = label_list.index(lbl)
    atom_red_coord += ["'%s' %.4f %.4f %.4f %d" % (
        lbl, r[0], r[1], r[2], i+1
    )]
    total_nelec += tbl_fhi_abinit[lbl][1]


with open("templates/gs_rt_response.inp") as fh:
    print(fh.read().format(
        SYSNAME = sysname,
        LX = norm(a_orth[0]),
        LY = norm(a_orth[1]),
        LZ = norm(a_orth[2]),
        NELEC = total_nelec,
        NELEM = len(label_list),
        NATOM = len(site_pos),
        IZATOM = ",".join(izatom),
        LLOC_PS = ",".join(lloc_ps),
        PSEUDO_FILE = ",".join(pseudo_file),
        ATOMIC_RED_COOR = "\n".join(atom_red_coord),
    )
)
