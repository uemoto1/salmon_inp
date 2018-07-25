#!/usr/bin/env python3
import os
import sys
import optparse

from core import cif_file
from core import supercell
from core import tool
from core import salmon_file

dir_templates = os.path.join(os.curdir, "templates")
dir_pptbl = os.path.join(os.curdir, "pptables")

def main():

    parser = optparse.OptionParser(
        usage = """usage: %prog [options] < ciffile.cif > salmonfile.inp"""
    )
    parser.add_option("-p", "--pptype", dest="pptype",
                      default="abinit-fhi", type=str, help="Type of Pseudopotentials ")
    parser.add_option("-t", "--template", dest="template",
                      default="gs_rt_response", help="Calculation mode")
    parser.add_option("--np", dest="np", type=int, default=100, help="Required number of grid inside potential")
    parser.add_option("--export-cif", dest="export_cif",
                      default="", help="Export CIF file of generated supercell")
    opts, args = parser.parse_args()

    if os.path.isfile(opts.template):
        file_template = opts.template
    else:
        file_template = os.path.join(dir_templates, "%s.inp" % opts.template)
                
    cif = cif_file.CIF()
    cif.loads(sys.stdin.read())
    
    a_orth, site_pos, site_lbl = supercell.search_supercell(
        cif.a_prim, cif.site_pos, cif.site_lbl
    )
    
    salmon = salmon_file.Salmon(
        cif.sysname, a_orth, site_lbl, site_pos,  
        os.path.join(dir_pptbl, "%s.json" % opts.pptype), opts.np
    )
    

    
    sys.stdout.write(salmon.dumps(file_template))
    
    if opts.export_cif:
        cif2 = cif_file.CIF(cif.sysname, a_orth, site_lbl, site_pos)
        with open(opts.export_cif, "w") as fh_cif:
            fh_cif.write(cif2.dumps())

if __name__ == "__main__":
    main()
