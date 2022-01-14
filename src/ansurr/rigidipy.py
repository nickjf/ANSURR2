#!/usr/bin/env python

import subprocess
import sys
import os
import json
import shutil
import glob

from importlib import resources
from sys import platform
from collections import Counter

from ansurr.functions import check_quiet_print

def calc_rigidity(pdb,quiet=False):

    if platform == "linux" or platform == "linux2":
        with resources.path("ansurr.bin", "calc_rigidity_gcc4.8.4.bin") as f:  # use this for the package
        #with resources.path("bin", "calc_rigidity_gcc4.8.4.bin") as f: 
            data_file_path = f
    elif platform == "darwin":
        with resources.path("ansurr.bin", "calc_rigidity_OSX.bin") as f:  # use this for the package
        #with resources.path("bin", "calc_rigidity_OSX.bin") as f:
            data_file_path = f
    elif platform == "win32":
        check_quiet_print(quiet,'ANSURR will not yet run on Windows')
        sys.exit(0)

    lib_location = str(os.path.dirname(data_file_path))+'/../'
    
    s = subprocess.run([data_file_path, pdb, lib_location],stdout=subprocess.PIPE, stderr=open(os.devnull, 'w'))
   
    if "Calculating rigidity" in str(s.stdout):
        return True
    else:
        return False

def make_monomers(pdbid,chains,model,resi_ref,output_dir):
    pdb = output_dir+'combined/'+pdbid+chains+'_'+model+'.pdb'
    prev_chain = 'XXX'
    for line in open(pdb,'r'):
        if 'ATOM' in line[:4] or 'HETATM' in line[:6]:
            chain = line[21]
            resi = int(line[22:26])
            orig_resi = resi - resi_ref[model][chain]['new_first'] + resi_ref[model][chain]['orig_first']
            if prev_chain == 'XXX':
                prev_chain = chain
                out = open(output_dir+pdbid+chain+'_'+model+'.pdb','w')
            elif chain != prev_chain:
                out.close()
                out = open(output_dir+pdbid+chain+'_'+model+'.pdb','w')
            out.write(line[0:22]+''.join([' ']*(4-len(str(orig_resi))))+str(orig_resi)+line[26:])
            prev_chain = chain
    out.close()

def parse_pdb_lines(pdb_lines,model,models,freeligands, nonstandardres, oligomers, conects, nonstandard, quiet=False):
    standard_res = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
    selected_pdb_lines = {}

    for chain in pdb_lines:
        std_res = 0
        for line in pdb_lines[chain]:
            resn = line[17:20]
            if resn in standard_res and line[:4] == 'ATOM': # helps diagnose non-standard residues from free ligands, non-standard residues are HETATMS that come with standard residues too
                std_res = 1
                break
        if std_res == 1:
            for line in pdb_lines[chain]:
                resn = line[17:20]
                if resn in standard_res:
                    if chain not in selected_pdb_lines:
                        selected_pdb_lines[chain] = [line]
                    else:
                        selected_pdb_lines[chain].append(line)
                elif nonstandardres:
                    if chain not in selected_pdb_lines:
                        selected_pdb_lines[chain] = ['ATOM  '+line[6:]]
                    else:
                        selected_pdb_lines[chain].append('ATOM  '+line[6:])
                    if resn not in nonstandard:
                        check_quiet_print(quiet," -> non-standard residue "+resn.replace(' ','')+" will be included in rigidity calculations")
                        nonstandard.append(resn)
                else:
                    if resn not in nonstandard:
                        check_quiet_print(quiet," -> found a non-standard residue ("+resn.replace(' ','')+") which are ignored by default. To include non-standard residues, re-run ANSURR with the -n flag")
                        nonstandard.append(resn)
        else:
            for line in pdb_lines[chain]:
                resn = line[17:20]
                if freeligands:
                    if chain not in selected_pdb_lines:
                        selected_pdb_lines[chain] = [line]
                    else:
                        selected_pdb_lines[chain].append(line)
                    if resn not in nonstandard:
                        check_quiet_print(quiet," -> free ligand "+resn.replace(' ','')+" will be included in rigidity calculations")
                        nonstandard.append(resn)
                else:
                    if resn not in nonstandard:
                        check_quiet_print(quiet," -> found a free ligand ("+resn.replace(' ','')+") which are ignored by default. To include free ligands, re-run ANSURR with the -l flag")
                        nonstandard.append(resn)
    for chain in selected_pdb_lines:
        if len(selected_pdb_lines[chain]) > 0:
            if model not in models:  
                models[model] = {chain:selected_pdb_lines[chain]}
            elif chain not in models[model]:
                models[model][chain] = selected_pdb_lines[chain]
            else:
                models[model][chain].extend(selected_pdb_lines[chain])
    return models

    
def extract_pdb(pdb_file, output_dir, freeligands=False, nonstandardres=False, oligomers=False, conects=False,quiet=False):

    pdb_in = open(pdb_file,'r')
    pdb = os.path.basename(os.path.splitext(pdb_file)[0])

    model = '1' # assume first model is model 1
    models = {}

    pdb_lines = {}
    chain = ''
    conect = []
    nonstandard = []

    pdbs_out = []

    conformer_chosen = ''
    conformer_rejects = []

    chain_output = {}

    for line in pdb_in:
        if line[:6] == 'EXPDTA':
            if 'NMR' not in line and 'SOLID' not in line: # some older solution NMR structures just referred to as NMR
                check_quiet_print(quiet,' -> CRITICAL ERROR this version of ANSURR is only for validating structures solved using SOLUTION NMR, exiting')
                quit()
        elif line[:5] == 'MODEL':
            model = line.split()[1]
        elif line[:4] == 'ATOM' or line[:6] == 'HETATM':
            chain = '' if line[21] == ' ' else line[21]
            conformer = '' if line[16] == ' ' else line[16]

            if conformer != '' and conformer_chosen == '':
                conformer_chosen = conformer

            if conformer == '' or conformer == conformer_chosen:
                if chain not in pdb_lines:
                    pdb_lines[chain] = [line]
                else:
                    pdb_lines[chain].append(line)
            else:
                conformer_rejects.append(conformer)

        elif 'TER' in line[:3] or 'END' in line[:3]: # parse after TER to help diagnose free ligands from non-sandard residues
            if len(pdb_lines) > 0:   
                models = parse_pdb_lines(pdb_lines,model,models,freeligands, nonstandardres, oligomers, conects,nonstandard,quiet=quiet)
                pdb_lines = {}
        elif 'CONECT'in line[:6] and conects == True:
            conect.append(line)

    if len(pdb_lines) > 0:   # catch any extra stuff in case structure doesn't end with TER or END, e.g. for CNS output
        models = parse_pdb_lines(pdb_lines,model,models,freeligands, nonstandardres, oligomers, conects,nonstandard,quiet=quiet)
        pdb_lines = {}

    for model in models: # this orders by resi number as sometimes pdbs are not ordered correctly (e.g. pdb 2lrl)
        for chain in models[model]:
            pdb_lines = models[model][chain]
            resi_pdb_lines = {}

            for line in pdb_lines:
                resi = int(line[22:26])
                if resi not in resi_pdb_lines:
                    resi_pdb_lines[resi] = [line]
                else:
                    resi_pdb_lines[resi].append(line) 

            pdb_lines = [resi_pdb_lines[i] for i in sorted(resi_pdb_lines)]
            models[model][chain] =  [item for sublist in pdb_lines for item in sublist]


    chains_done = []
    resi_ref = {}
    old_new_atom_num = {} # 
    for model in models:
        chains = ''
        if oligomers:
            if len(models[model]) > 1:
                if not os.path.exists(output_dir+'/other_output/extracted_pdbs/combined'):
                    os.makedirs(output_dir+'/other_output/extracted_pdbs/combined')
                for chain in models[model]:
                    chains += chain
                out = open(output_dir+'/other_output/extracted_pdbs/combined/'+pdb+chains+'_'+model+'.pdb','w')
                pdbs_out.append(output_dir+'/other_output/extracted_pdbs/combined/'+pdb+chains+'_'+model+'.pdb')
                count = 1
                num = 0
                resi_ref[model] = {}
                for chain in models[model]:
                    prev_resi = -99999
                    resi_ref[model][chain] = {'orig_first':int(models[model][chain][0][22:26]),'orig_last':int(models[model][chain][-1][22:26]),'new_first':'','new_last':''}
                    for l in models[model][chain]:
                        old_new_atom_num[l[6:11].strip()] = str(count)
                        resi = int(l[22:26])
                        if prev_resi != -99999:
                            num += resi - prev_resi
                        else:
                            num += 1
                        prev_resi = resi
                        if resi_ref[model][chain]['new_first'] == '':
                            resi_ref[model][chain]['new_first'] = num
                        out.write(l[0:6]+format(str(count)," >5s")+l[11:22]+''.join([' ']*(4-len(str(num))))+str(num)+l[26:].replace('\n','')+'\n') #this makes sure there is a single new line char at end (needed for FIRST to run correctly)
                        count +=1 
                    resi_ref[model][chain]['new_last'] = num
                out.write('END\n')
                for c in conect:
                    corrected_c = ''
                    for i in c.split()[1:]:
                        if i in old_new_atom_num:
                            corrected_c += ' '+old_new_atom_num[i]
                    if len(corrected_c.split()) > 1:
                        out.write('CONECT '+corrected_c+'\n')
                out.close()
                make_monomers(pdb,chains,model,resi_ref,output_dir+'/other_output/extracted_pdbs/')
                if chains not in chains_done:
                    check_quiet_print(quiet," -> chains "+chains+" combined into a single structure to calculate flexibility ")   
                    chains_done.append(chains)
            else:
                oligomers = False
                check_quiet_print(quiet,' -> found only a single chain ('+str([i for i in models[model]][0])+'), no need to combine')
        if oligomers == False:
            for chain in models[model]:
                chains += chain
                if not os.path.exists(output_dir+'/other_output/extracted_pdbs/'):
                    os.makedirs(output_dir+'/other_output/extracted_pdbs/')
                out = open(output_dir+'/other_output/extracted_pdbs/'+pdb+chain+'_'+model+'.pdb','w')
                pdbs_out.append(output_dir+'/other_output/extracted_pdbs/'+pdb+chain+'_'+model+'.pdb')
                count = 1
                for l in models[model][chain]:
                    old_new_atom_num[l[6:11].strip()] = str(count)
                    out.write(l[0:6]+format(str(count)," >5s")+l[11:])
                    count +=1 
                out.write('END\n')
                for c in conect:
                    corrected_c = ''
                    for i in c.split()[1:]:
                        if i in old_new_atom_num:
                            corrected_c += ' '+old_new_atom_num[i]
                    if len(corrected_c.split()) > 1:
                        out.write('CONECT '+corrected_c+'\n')
                out.close()
            if len(chains) > 1 and chains not in chains_done:
                check_quiet_print(quiet," -> found "+str(len(chains))+" chains, which are extracted seperately by default. To combine chains when calculating flexibility, re-run ANSURR with the -o flag")
                chains_done.append(chains)

    if conformer_chosen != '':
        err_string = " -> found alternative locations, using: "+conformer_chosen+", ignoring: "
        for c in sorted(list(set(conformer_rejects))):
            err_string+=c+', '
        check_quiet_print(quiet,err_string[:-2])

    if len(conect) > 0:
        check_quiet_print(quiet," -> read in "+str(len(conect))+" CONECT records")
    if oligomers:
        json.dump(resi_ref, open("resi_ref.tmp",'w'))

    return pdbs_out, oligomers


def rigid_decomp(path_to_pdb,path_to_decomp,output_dir,chain_output=''):

    pdb = open(path_to_pdb,'r')
    decomp = open(path_to_decomp+'decomp_list','r')

    class res(object):
        _registry = []
        def __init__(self,i,name,atom):
            self._registry.append(self)
            self.i = i
            self.name = name
            self.CA = atom
            self.clusters = []
            self.energy = 0

    # get atom numbers for each CA atom from pdb, also get chain for chain_output_rigidity
    for line in pdb:
        if line[:4] == 'ATOM' or line[:6] == 'HETATM':
            atom_number = int(line[6:11])
            atom_name = line[12:16]
            resn = line[17:20].replace(' ','')
            resi = int(line[22:26])
            chain = line[21]
            if atom_name == ' CA ': # this is occasionally a problem for ligands which have a CA atom, however, it's difficult to spot a ligand CA atom from a non-std residue CA atom, and doesn't effect scoring etc. Leave for now
                r = res(resi,resn,atom_number)

    # import decomp_list         
    data = []
    all_data = []
    energy = []
    rigid_clusters = []
    for line in decomp:
        if 'HEADER' not in line[:6]:
            if 'END' in line[:3]:
                rc = []
                data = data[:atom_number] # remove overflow
                all_data.append(data)
                cluster_count = Counter(data)
                for c in cluster_count:
                    if int(cluster_count[c]) >= 15:
                        rc.append(c)
                rigid_clusters.append(rc)
                data = []
            elif 'A:' in line[:2]:
                energy.append(float(line[10:21]))
            elif 'A:' not in line[:2] and 'B:' not in line[:2]:
                line = line.replace('\n','').split(':')
                for l in line:
                    data.append(int(l))  

    # append rigid clusters to residues
    for d in all_data:
         for r in res._registry:
            r.clusters.append(d[r.CA-1])

    # find rigid clusters
    for x in range(len(all_data)):
        for r in res._registry:
            if r.clusters[x] in rigid_clusters[x]:
                r.energy = energy[x]

    #write output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    out = open(output_dir+os.path.basename(os.path.splitext(path_to_pdb)[0])+'.decomp','w')

    if chain_output != '':
        if chain not in chain_output:
            chain_output[chain] = [output_dir+os.path.basename(os.path.splitext(path_to_pdb)[0])+'.decomp']
        else:
            chain_output[chain].append(output_dir+os.path.basename(os.path.splitext(path_to_pdb)[0])+'.decomp')

    for r in res._registry:
        out.write(str(r.i) + ' ' + r.name +' ' + str(r.energy) + '\n')
    out.close()

    if chain_output != '':
        return chain_output


def split_decomp(decomp,pdb_file,output_dir,chain_output):

    resi_ref = json.load(open("resi_ref.tmp"))

    decomp_num = sys.argv[1].split('/')[-1]
    num = decomp.split('.decomp')[0].split('_')[-1]

    chains = ''
    for chain in resi_ref[num]:
        chains += chain

    for chain in resi_ref[num]:
        out = open(output_dir+pdb_file+chain+'_'+num+'.decomp','w')
        orig_first = resi_ref[num][chain]['orig_first']
        orig_last = resi_ref[num][chain]['orig_last']
        new_first = resi_ref[num][chain]['new_first']
        new_last = resi_ref[num][chain]['new_last']
        for line in open(decomp,'r'):
            resi= int(line.split()[0]) 
            if resi >= new_first and resi <= new_last:
                out.write(str(resi - new_first + orig_first)+ ' '+line.split()[1]+' '+line.split()[2]+'\n')
            if resi == new_last:
                out.close()
                break

        if chain not in chain_output:
            chain_output[chain] = [output_dir+pdb_file+chain+'_'+num+'.decomp']
        else:
            chain_output[chain].append(output_dir+pdb_file+chain+'_'+num+'.decomp')

    return chain_output
    


