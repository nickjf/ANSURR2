#!/usr/bin/env python

import sys
import os
import json
import argparse
import shutil
import glob
import re

from ansurr import rci_nef
from ansurr import rigidipy
from ansurr import compare
from ansurr import cyrange
from ansurr import plot2D
from ansurr import dssp
from ansurr import export_nef

from ansurr.functions import check_quiet_print
from datetime import datetime

sys.tracebacklimit = 0
ansurr_version = "2.0.47"

natsort = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)]

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pdb", type=str, help="input PDB file",required=True)
    parser.add_argument("-s", "--shifts", type=str, help="input shifts file in NEF format",required=True)
    parser.add_argument("-l", "--lig", help=" include free ligands when computing flexibility", action="store_true")
    parser.add_argument("-n", "--nonstd", help="include non-standard residues when computing flexibility", action="store_true")
    parser.add_argument("-m", "--min", help="only output the ANSURR scores in a text file", action="store_true")
    parser.add_argument("-o", "--olig", help="combine chains into a single structure when calculating flexibility", action="store_true")
    parser.add_argument("-q", "--quiet", help="suppress output to the terminal", action="store_true")
    parser.add_argument("-r", "--reref", help="re-reference shifts using PANAV before calculating RCI (requires Java)",action="store_true")
    parser.add_argument("-w", "--welldef", help="compute ANSURR scores for only well-defined residues identified by CYRANGE", action="store_true")
    parser.add_argument("-v", "--version", action='version', version='\nANSURR | Accuracy of NMR Structures Using RCI and Rigidity v'+ansurr_version+' | citation: https://doi.org/10.1038/s41467-020-20177-1')
    args = parser.parse_args()

    path_to_shifts = args.shifts
    path_to_pdb = args.pdb # need to check that file looks like a PDB file
    
    shift_file = os.path.basename(os.path.splitext(path_to_shifts)[0])
    pdb_file = os.path.basename(os.path.splitext(path_to_pdb)[0])

    now = datetime.now()
    check_quiet_print(args.quiet,"ANSURR v"+ansurr_version+" "+now.strftime("%d/%m/%Y %H:%M:%S"))
    check_quiet_print(args.quiet," -> pdb: "+path_to_pdb)
    check_quiet_print(args.quiet," -> shifts: "+path_to_shifts)

    if args.min:
        output_dir = 'ansurr_tmp'
    else:
        output_dir = os.path.join(os.getcwd(), pdb_file+'_'+shift_file)

    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)
        if not args.min:
            check_quiet_print(args.quiet," -> overwriting output directory: "+output_dir)
    else:
        if not args.min:
            check_quiet_print(args.quiet," -> output directory: "+output_dir)
    os.makedirs(output_dir)   

    check_quiet_print(args.quiet,"calculating RCI")
    set_chain_output_rci = rci_nef.calc_RCI(path_to_shifts,output_dir+'/other_output/',panav=args.reref,quiet=args.quiet)

    check_quiet_print(args.quiet,'extracting models from '+path_to_pdb)
    pdbs,args.oligomer  = rigidipy.extract_pdb(path_to_pdb,output_dir,freeligands=args.lig,nonstandardres=args.nonstd,oligomers=args.olig,quiet=args.quiet)
    if len(pdbs) > 0:
        check_quiet_print(args.quiet,' -> extracted '+str(len(pdbs))+' model(s)')
    else:
        check_quiet_print(args.quiet,' -> CRITICAL ERROR failed to extract any models, exiting')
        sys.exit(1)

    for pdb in glob.glob(output_dir+'/other_output/extracted_pdbs/*.pdb'):
        dssp.calc_secondary_structure(pdb,output_dir+'/other_output/dssp/',quiet=args.quiet)

    if args.welldef:
        args.welldef = cyrange.calc_welldefined(path_to_pdb,output_dir,args.quiet)

    check_quiet_print(args.quiet,'computing the flexibility of each model')
    chain_output_rigidity = {}
    for pdb in sorted(pdbs,key=natsort):
        check_quiet_print(args.quiet," -> "+os.path.basename(os.path.splitext(pdb)[0]),end='')
        if rigidipy.calc_rigidity(pdb,quiet=args.quiet):
            check_quiet_print(args.quiet,' DONE')
            if args.oligomer:
                rigidipy.rigid_decomp(pdb,output_dir+'/other_output/extracted_pdbs/combined/',output_dir+'/other_output/FIRST/')
                try:
                    os.remove(output_dir+'/other_output/extracted_pdbs/combined/decomp_list')  # important to delete decomp_list
                except OSError:
                    pass
            else:
                chain_output_rigidity = rigidipy.rigid_decomp(pdb,output_dir+'/other_output/extracted_pdbs/',output_dir+'/other_output/FIRST/',chain_output=chain_output_rigidity)
                try:
                    os.remove(output_dir+'/other_output/extracted_pdbs/decomp_list')  # important to delete decomp_list
                except OSError:
                    pass
        else:
            check_quiet_print(args.quiet,' FAILED')

    if len(glob.glob(output_dir+'/other_output/FIRST/*.decomp')) == 0: # check FIRST worked
        check_quiet_print(args.quiet,' -> CRITICAL ERROR FIRST failed to run for any models, exiting')
        sys.exit(1)

    if args.oligomer:
        for pdb in glob.glob(output_dir+"/other_output/FIRST/*.decomp"):
            chain_output_rigidity = rigidipy.split_decomp(pdb,pdb_file,output_dir+"/other_output/FIRST/",chain_output_rigidity)
            os.remove(pdb)

        try:
            os.remove("resi_ref.tmp")
        except:
            pass

    decomps = glob.glob(output_dir+'/other_output/FIRST/*.decomp')

    check_quiet_print(args.quiet,'computing ANSURR scores')

    nef_output = {'per_model':{},'per_residue':{}}
    for chain in chain_output_rigidity:
        for rigidity in sorted(chain_output_rigidity[chain],key=natsort):
            found = 0
            for rci_set in set_chain_output_rci:
                if chain in set_chain_output_rci[rci_set]:
                    for rci in sorted(set_chain_output_rci[rci_set][chain],key=natsort):
                        nef_output = compare.compare_rci_rigidity(rigidity,rci,pdb_chain=chain,rci_set=rci_set,rci_chain=chain,nef_output=nef_output,output_dir=output_dir,cyrange=args.welldef,quiet=args.quiet)
                        found = 1
            if found == 0:  
                for rci_set in set_chain_output_rci:
                    for rci_chain in sorted(set_chain_output_rci[rci_set],key=natsort):
                        for rci in set_chain_output_rci[rci_set][rci_chain]:
                            nef_output = compare.compare_rci_rigidity(rigidity,rci,pdb_chain=chain,rci_set=rci_set,rci_chain=rci_chain,nef_output=nef_output,output_dir=output_dir,cyrange=args.welldef,quiet=args.quiet)  
        
    if len(glob.glob(output_dir+'/ANSURR_output/out/*.out')) > 0:
        
        export_nef.export(nef_output,output_dir=output_dir+'/ANSURR_output/'+pdb_file+'_'+shift_file,ansurr_version=ansurr_version,pdb_file=os.path.abspath(path_to_pdb),shift_file=os.path.abspath(path_to_shifts),reref=args.reref,lig=args.lig,nonstd=args.nonstd,olig=args.olig)

        if os.path.exists(output_dir+'/ANSURR_output/scores.out'):
            if not args.min:
                plot2D.plot_scores(output_dir+'/ANSURR_output/',cyrange=args.welldef)
            check_quiet_print(args.quiet,"summary of ANSURR scores")
            for line in open(output_dir+'/ANSURR_output/scores.out','r'):
                check_quiet_print(args.quiet," -> "+line.strip())
    else:
        check_quiet_print(args.quiet," -> CRITICAL ERROR failed to calculate any validation scores, exiting")
        sys.exit(1)

    if args.min:
        check_quiet_print(args.quiet,'ANSURR scores saved to scores_'+pdb_file+'_'+shift_file+'.out')
        shutil.move(output_dir+'/ANSURR_output/scores.out','scores_'+pdb_file+'_'+shift_file+'.out')
        shutil.rmtree(output_dir)

    check_quiet_print(args.quiet,"")
    sys.exit(0)

      
if __name__ == "__main__":
    main()

