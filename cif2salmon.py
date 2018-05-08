#!/usr/bin/env python3
import sys
import os

from core import cif_file
from core import supercell
from core import tool
from core import salmon_file
from numpy.linalg import norm

import optparse

dir_templates = os.path.join(os.curdir, "templates")
dir_pptbl = os.path.join(os.curdir, "pptables")

def main():

    parser = optparse.OptionParser()
    parser.add_option("-p", "--pptype", dest="pptype",
                      default="abinit-fhi", type=str, help="Type of Pseudopotentials")
    parser.add_option("-t", "--template", dest="template",
                      default="gs_rt_response", help="Calculation mode")
    parser.add_option("--export-cif", dest="export_cif",
                      default="", help="Export CIF file of generated supercell")
    opts, args = parser.parse_args()


    cif = cif_file.CIF()
    cif.loads(sys.stdin.read())
    
    a_orth, site_pos, site_lbl = supercell.search_supercell(
        cif.a_prim, cif.site_pos, cif.site_lbl
    )
    
    salmon = salmon_file.Salmon(
        cif.sysname, a_orth, site_lbl, site_pos,  
        os.path.join(dir_pptbl, "%s.json" % opts.pptype),
    )
    
    
    print(salmon.dumps(
        os.path.join(dir_templates, "%s.inp" % opts.template)
    ))
    


main()
