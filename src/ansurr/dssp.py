import sys
import os
import subprocess

from sys import platform
from importlib import resources

from ansurr.functions import check_quiet_print

def calc_secondary_structure(pdb,output_dir='',quiet=False):

    if platform == "linux" or platform == "linux2":
        with resources.path("ansurr.bin", "dssp.bin") as f:  # use this for the package
        #with resources.path("bin", "dssp.bin") as f: 
            dssp_file_path = f
    #if platform == "darwin":
    #    with resources.path("ansurr.bin", "dsspOSX.bin") as f:  # use this for the package
    #        dssp_file_path = f

    if platform != "darwin":
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        out = open(output_dir+os.path.basename(os.path.splitext(pdb)[0])+".dssp","w")

        dssp_output = subprocess.run([dssp_file_path, pdb],stdout=subprocess.PIPE, universal_newlines=True, stderr=open(os.devnull, 'w'))

        start = 0
        for line in dssp_output.stdout.split('\n'):
            if start == 1:
                try:
                    resi = str(int(line[5:10]))
                    resn = line[13]
                    ss = line[16]
                    if ss == ' ':
                        ss = '.'
                    out.write(resi+' '+resn+' '+ss+'\n')
                except:
                    pass

            elif "#  RESIDUE AA" in line:
                start = 1
        out.close()

def main():

    calc_secondary_structure(sys.argv[1])
    sys.exit(0)
   
if __name__ == "__main__":
    main()