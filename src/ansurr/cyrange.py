#!/usr/bin/env python

import sys
import os
import json
import subprocess

from sys import platform
from importlib import resources

from ansurr.functions import check_quiet_print

def calc_welldefined(pdb,output_dir='',quiet=False):

    if platform == "linux" or platform == "linux2":
        with resources.path("ansurr.bin", "cyrange.bin") as f:  # use this for the package
        #with resources.path("bin", "cyrange.bin") as f: 
            cyrange_file_path = f
    elif platform == "darwin":
        with resources.path("ansurr.bin", "cyrange_OSX.bin") as f:  # use this for the package
        #with resources.path("bin", "cyrange.bin") as f: 
            cyrange_file_path = f

    cyrange_output = subprocess.run([cyrange_file_path, pdb],stdout=subprocess.PIPE, universal_newlines = True, stderr=open(os.devnull, 'w'))

    cyrange = {'default chain':[]}

    for cyrange_out in cyrange_output.stdout.split('\n'):
        if "Optimal range" in cyrange_out:
            cyrange_out = cyrange_out[18:]
            gaps = int(cyrange_out.split(':')[1].split(',')[1].split()[0])
            if gaps == 0:
                if cyrange_out[0].isalpha():
                    chain = cyrange_out[0]
                    first_resi = int(cyrange_out.split(':')[0].split('..')[0][1:])
                    last_resi = int(cyrange_out.split(':')[0].split('..')[1][1:])
                    if chain not in cyrange:
                        cyrange[chain] = []
                    cyrange[chain].extend(list(range(first_resi,last_resi+1)))
                else:
                    first_resi = int(cyrange_out.split(':')[0].split('..')[0])
                    last_resi = int(cyrange_out.split(':')[0].split('..')[1])
                    cyrange['default chain'].extend(list(range(first_resi,last_resi+1)))
            else:
                for r in cyrange_out.split(':')[0].split(','):
                    r = r.strip()
                    if '..' in r:
                        if r[0].isalpha():
                            chain = r[0]
                            first_resi = int(r.split(':')[0].split('..')[0][1:])
                            last_resi = int(r.split(':')[0].split('..')[1][1:])
                            if chain not in cyrange:
                                cyrange[chain] = []
                            cyrange[chain].extend(list(range(first_resi,last_resi+1)))
                        else:
                            first_resi = int(r.split(':')[0].split('..')[0])
                            last_resi = int(r.split(':')[0].split('..')[1])
                            cyrange['default chain'].extend(list(range(first_resi,last_resi+1)))
                    else:
                        if r[0].isalpha():
                            chain = r[0]
                            cyrange[chain].extend([int(r)])
                        else:
                            cyrange['default chain'].extend([int(r)])

    chain_resi = {}
    for chain in cyrange:
        resi = len(cyrange[chain])
        if resi > 0:
            chain_resi[chain] = resi
      
    if len(chain_resi) > 0: 
        if len(chain_resi) == 1:
            check_quiet_print(quiet,' -> identified '+str(chain_resi['default chain'])+' well-defined residues')
        else:
            msg=' -> identified'
            for chain in chain_resi:
                msg+=' '+str(chain_resi[chain])+' ('+chain+'),'
            msg=msg[:-1]+" well-defined residues"
            print(msg)
        
        if not os.path.exists(output_dir+'/other_output/cyrange/'):
            os.makedirs(output_dir+'/other_output/cyrange/')
        with open(output_dir+'/other_output/cyrange/cyrange.json', 'w') as outfile:
            json.dump(cyrange, outfile)
        return True

    else:
        check_quiet_print(quiet," -> ERROR CYRANGE failed to identify well-defined residues")
        return False

def main():

    cyrange(sys.argv[1])
    sys.exit(0)
   
if __name__ == "__main__":
    main()