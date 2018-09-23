import os
import sys
import json
import numpy as np
from numpy.linalg import norm

class Salmon:
    
    
    def __init__(self, file_template, sysname, a_orth, site_lbl, site_pos, file_pptbl, np=100):
        self.sysname = sysname
        self.a_orth = a_orth
        self.site_lbl = site_lbl
        self.site_pos = site_pos
        self.pp_list = list(set(self.site_lbl))
        self.np = np
        
        if not os.path.isfile(file_pptbl):
            sys.stderr.write("# Potential '%s' is not found!" % file_pptbl)
            sys.exit(-1) 
        
        with open(file_pptbl, "r") as fh_pptbl:
            self.pptbl = json.load(fh_pptbl)
            
        if not os.path.isfile(file_template):
            sys.stderr.write("# Template '%s' is not found!" % file_pptbl)
            sys.exit(-1) 

        with open(file_template, "r") as fh_template:
            self.template = fh_template.read()
            
        self.pp_url = [
            (self.pptbl[i]['url'], self.pptbl[i]["pseudo_file"])
            for i in self.pp_list
        ]
        
    def dumps(self):
        # Search Prefered Delta
        r_cutoff = min([self.pptbl[i]["r_cutoff"] for i in self.pp_list])
        v_cutoff = (4.0 * np.pi / 3.0) * r_cutoff ** 3
        delta = np.cbrt(v_cutoff / self.np)
        
        input_data = self.template.format(
            SYSNAME = self.sysname,
            LX = norm(self.a_orth[0]),
            LY = norm(self.a_orth[1]),
            LZ = norm(self.a_orth[2]),
            NX = int(np.ceil(norm(self.a_orth[0]) / delta / 4)) * 4,
            NY = int(np.ceil(norm(self.a_orth[1]) / delta / 4)) * 4,
            NZ = int(np.ceil(norm(self.a_orth[2]) / delta / 4)) * 4,
            NATOM = len(self.site_lbl),
            NELEM = len(self.pp_list),
            IZATOM = ",".join([
                repr(self.pptbl[i]["izatom"]) for i in self.pp_list
            ]),
            LLOC_PS = ",".join([
                repr(self.pptbl[i]["lloc_ps"]) for i in self.pp_list
            ]),
            PSEUDO_FILE = ",".join([
                repr(self.pptbl[i]["pseudo_file"]) for i in self.pp_list
            ]),
            NELEC = sum([
                self.pptbl[lbl]["nelec"] for lbl in self.site_lbl
            ]),
            ATOMIC_RED_COOR = "\n".join(
                "'%s' %.4f %.4f %.4f %d" % (
                    lbl, pos[0], pos[1], pos[2], self.pp_list.index(lbl) + 1
                )
                for lbl, pos in zip(self.site_lbl, self.site_pos)
            ),
        )
        
        return input_data
            
