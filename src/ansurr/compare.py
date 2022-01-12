#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib
import json

from matplotlib import pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import percentileofscore
from importlib import resources

from ansurr.functions import check_quiet_print

matplotlib.use('Agg')

def compare_seq(seq1,seq2):
    score = 0
    count = 0
    for s1 in seq1:
        if s1 == seq2[count]:
            score += 1
        count += 1
    return score

def fix_missing_res(d,report_break,quiet=False):
    count = d['resi'][0]
    pos = 0
    to_delete = []
    for i in d['resi']:
        if i > count:
            if report_break == 1:
                check_quiet_print(quiet,'*WARNING break in PDB numbering at resi '+str(count-1)+'* ', end='')
                report_break = 2
            for k in d:
                if k == 'resi':
                    d[k].insert(pos,count)
                elif k == 'resn':
                    d[k].insert(pos,'XXX')
                else:
                    d[k].insert(pos,np.nan)
        elif i < count:
            if report_break == 1:
                check_quiet_print(quiet,'*WARNING unexpected PDB numbering at resi '+str(count-1)+', renumbering* ', end='')
                report_break = 2
            d['resi'][pos] = count   
        elif report_break == 2:
            report_break = 1
        count += 1
        pos += 1
    return d

def calc_score(seq1,seq2):
    score = 0
    count = 0
    for s1 in seq1:
        if s1 == seq2[count]:
            score += 1
        count += 1
    return score

def find_breaks(nums):
    ranges = sum((list(t) for t in zip(nums, nums[1:]) if t[0]+1 != t[1]), [])
    iranges = iter(nums[0:1] + ranges + nums[-1:])
    breaks = [[n,int(next(iranges))+1] for n in iranges]
    return breaks

def align(RCI,FIRST,cut_off=0.5,quiet=False):
    quit_lowseqID = 0
    if len(RCI['resn']) <= len(FIRST['resn']): # "a" is always the shorter seq
        a = RCI['resn']
        b = FIRST['resn']
    else:
        a = FIRST['resn']
        b = RCI['resn']

    max_score_N = 0
    for i in range(len(b)): # start from N
        a_test = a[-i-1:len(b)+1]
        a_test.extend([np.nan]*(len(b)-i-1))
        a_test = a_test[::-1]
        a_test.extend([np.nan]*(len(b)-len(a_test)))
        a_test = a_test[::-1]  
        score = calc_score(a_test,b)
        if score > max_score_N:
            max_score_N = score
            chosen_i = i

    max_score_C = max_score_N # then start from C and look for better fit
    for i in range(1,len(a)):
        a_test = a[:-i]
        a_test = a_test[::-1]
        a_test.extend([np.nan]*(len(b)-len(a_test)))
        a_test = a_test[::-1]   
        score = calc_score(a_test,b)
        if score > max_score_C:
            max_score_C = score
            chosen_i = i

    if max_score_N >= max_score_C:
        if len(RCI['resi']) <= len(FIRST['resi']):
            for r in RCI:  
                RCI[r] = RCI[r][-chosen_i-1:len(b)+1] 
                RCI[r].extend([np.nan]*(len(b)-chosen_i-1))
                RCI[r] = RCI[r][::-1]   
                RCI[r].extend([np.nan]*(len(b)-len(RCI[r])))
                RCI[r] = RCI[r][::-1]
        else:
            for f in FIRST:
                FIRST[f] = FIRST[f][-chosen_i-1:len(b)+1] 
                FIRST[f].extend([np.nan]*(len(b)-chosen_i-1)) 
                FIRST[f] = FIRST[f][::-1]   
                FIRST[f].extend([np.nan]*(len(b)-len(FIRST[f])))
                FIRST[f] = FIRST[f][::-1]
        num_resn = len([i for i in FIRST['resn'] if i != 'XXX'])
        if max_score_N < cut_off * num_resn:
            check_quiet_print(quiet,'sequence identity ('+str(round(100*max_score_N/num_resn,1))+'%) is below cut-off ('+str(100*cut_off)+'%), skipping')
            quit_lowseqID = 1
    else:
        if len(RCI['resi']) <= len(FIRST['resi']):
            for r in RCI:
                RCI[r] = RCI[r][:-chosen_i]
                RCI[r] = RCI[r][::-1]
                RCI[r].extend([np.nan]*(len(b)-len(RCI[r])))
                RCI[r] = RCI[r][::-1]
        else:
            for f in FIRST:
                FIRST[f] = FIRST[f][:-chosen_i]
                FIRST[f] = FIRST[f][::-1]
                FIRST[f].extend([np.nan]*(len(b)-len(FIRST[f])))
                FIRST[f] = FIRST[f][::-1]
        num_resn = len([i for i in FIRST['resn'] if i != 'XXX'])
        if max_score_C < cut_off * num_resn:
            check_quiet_print(quiet,'sequence identity ('+str(round(100*max_score_C/num_resn,1))+'%) is below cut-off ('+str(100*cut_off)+'%), skipping')
            quit_lowseqID = 1
    RCI['resi'] = FIRST['resi']
    return RCI, FIRST, quit_lowseqID

def trim(res1,res2): # lazy variable names fix this!
    s = 0
    e = 0
    sc = 0
    ec = 0
    for i in enumerate(res1):
        if str(i[1]) != 'nan': 
            if sc == 3:
                if i[1] != res2[i[0]]:
                    ec += 1
                else:
                    ec = 0
                if ec == 3:
                    e = i[0]
            else:
                if i[1] == res2[i[0]]:
                    sc += 1
                else:
                    sc = 0
                if sc == 3:
                    s = i[0]
    s = s-2
    if ec == 3:
        e = e-2
    else:
        e = len(res1) - ec
    return s,e
   
def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

def get_segments(resi):
    prev_res = min(resi) if resi else None
    segments = list()
    for number in enumerate(resi):
        if number[1] != prev_res+1:
            segments.append([number[0]])
        elif len(segments[-1]) > 1:
            segments[-1][-1] = number[0]
        else:
            segments[-1].append(number[0])
        prev_res = number[1]
    return segments

def smooth(resi,data):
    data_smoothed = []
    segs = get_segments(resi)
    for s in segs:
        if len(s) > 2:
            data_smoothed_temp = movingaverage(data[s[0]:s[1]+1],3)
            N = (data[s[0]] + data[s[0]+1]) / 2.0
            data_smoothed_temp[0] = N
            C = (data[s[1]] + data[s[1]-1]) / 2.0
            data_smoothed_temp[-1] = C
            data_smoothed.extend(data_smoothed_temp)
        elif len(s) == 2:
            data_smoothed_temp = movingaverage(data[s[0]:s[1]+1],2)
            N = (data[s[0]] + data[s[0]+1]) / 2.0
            data_smoothed_temp[0] = N
            C = (data[s[1]] + data[s[1]-1]) / 2.0
            data_smoothed_temp[-1] = C
            data_smoothed.extend(data_smoothed_temp)
        else:
            data_smoothed.append(data[s[0]])
    return data_smoothed

def calc_RMSD(RCI_score_for_corr,FIRST_score_for_corr):
    diff = []
    for i in enumerate(RCI_score_for_corr):
        r = i[1]
        f = FIRST_score_for_corr[i[0]]
        diff.append(np.square(r-f))
    return np.sqrt(np.mean(diff))

def rescale_FIRST(FIRST):
    K = 315777.09
    T = 298.15
    Hartree = 0.00038
    FIRST['score'] = [np.exp((4.2* x * Hartree)*K/T) if not np.isnan(x) else np.nan for x in FIRST['score']]
    FIRST_nan = [i for i,x in enumerate(FIRST['score']) if np.isnan(x)]
    FIRST['score'] = smooth([FIRST['resi'][i] for i,x in enumerate(FIRST['score']) if not np.isnan(x)],[i for i in FIRST['score'] if not np.isnan(i)])
    for n in FIRST_nan:
        FIRST['score'].insert(n,np.nan)
    return FIRST['score']

def rescale_RCI(RCI_score): 
    offset = 0.024
    scale = 0.2-0.024
    RCI_score = [min(max((i-offset)/scale,0.0),1.0) for i in RCI_score]
    return RCI_score

def all_same(items):
    return all(x == items[0] for x in items)

def compare_rci_rigidity(input_first,input_rci,pdb_chain,rci_set,rci_chain,nef_output,output_dir='',cyrange=False,quiet=False):
    
    # import data
    FIRST_in = open(input_first,'r')
    RCI_in = open(input_rci,'r')
    PDB_ID = os.path.basename(os.path.splitext(input_first)[0])
    SHIFT_ID= os.path.basename(os.path.splitext(input_rci)[0])

    MODEL_ID = PDB_ID.split('_')[-1]

    pdb_chain = pdb_chain if pdb_chain !=' ' else '.'

    if cyrange:
        with open(output_dir+'/other_output/cyrange/cyrange.json') as json_file:
            welldefined_residues = json.load(json_file)

    with resources.path("ansurr.lib", "benchmarks.json") as f:
        benchmark_file_path = f
    #benchmark_file_path = 'lib/benchmarks.json'

    with open(benchmark_file_path) as json_file:
        benchmarks = json.load(json_file)

    check_quiet_print(quiet," -> "+PDB_ID + '|'+ SHIFT_ID+' ',end='')

    # secondary structure
    ss_dict = {}
    try:
        secondary_structure_in = open(output_dir+'/other_output/dssp/'+PDB_ID+'.dssp','r')
        for line in secondary_structure_in:
            line = line.split()
            try:
                ss_dict[int(line[0])] = line[2]
            except:
                pass
    except:
        pass

    # read in data
    RCI = {'resi':[],'resn':[],'score':[],'shifts':[],'shift_types':[]}
    FIRST = {'resi':[],'resn':[],'score':[],'ss':[]}

    for line in RCI_in:
        line = line.split()
        RCI['resi'].append(int(line[0]))
        RCI['resn'].append(line[1])
        RCI['score'].append(float(line[2]))
        if line[1] == 'GLY' or line[1] == 'PRO':
            RCI['shifts'].append(float(line[3])/5.0)
            RCI['shift_types'].append(line[4])
        else:
            RCI['shifts'].append(float(line[3])/6.0)
            RCI['shift_types'].append(line[4])
            
    for line in FIRST_in:
        line = line.split()
        FIRST['resi'].append(int(line[0]))
        FIRST['resn'].append(line[1])
        FIRST['score'].append(float(line[2]))
        try:
            FIRST['ss'].append(ss_dict[int(line[0])])
        except:
            FIRST['ss'].append(np.nan)

    if cyrange:
        FIRST['welldefined'] = []
        if len(welldefined_residues) == 1:
            chain = 'default chain'
        else:
            chain = pdb_chain

        for resi in FIRST['resi']:
            if resi in welldefined_residues[chain]:
                FIRST['welldefined'].append('1')
            else:
                FIRST['welldefined'].append('0')
    else:
        FIRST['welldefined'] = ['0']*len(FIRST['resi'])
        
    # fix missing residues
    FIRST = fix_missing_res(FIRST,1,quiet=quiet)
    RCI = fix_missing_res(RCI,0,quiet=quiet)

    # align
    RCI,FIRST,quit_lowseqID = align(RCI,FIRST,quiet=quiet)

    if quit_lowseqID == 0:

        # trim - look for where first 3 residues and last 3 residues in aligned RCI/FIRST output are and make these the beginning/end
        start,end = trim(RCI['resn'],FIRST['resn'])
        for r in RCI:
            RCI[r] = RCI[r][start:end] 
        for f in FIRST:
            FIRST[f] = FIRST[f][start:end] 

        # rescale
        FIRST['score'] = rescale_FIRST(FIRST)
        #RCI['score'] = rescale_RCI(RCI['score']) # RCI now arrives rescaled

        # various stats
        RCI_nan = [i for i,x in enumerate(RCI['score']) if np.isnan(x)]
        FIRST_nan = [i for i,x in enumerate(FIRST['score']) if np.isnan(x)]
        RCI_noRC = [x for i,x in enumerate(RCI['score']) if i not in FIRST_nan and not np.isnan(x) and RCI['shift_types'][i] != 'RC']
        FIRST_noRC = [x for i,x in enumerate(FIRST['score']) if i not in RCI_nan and not np.isnan(x) and RCI['shift_types'][i] != 'RC']

        RMSD_noRC =  calc_RMSD(RCI_noRC,FIRST_noRC)
        RMSD_score = 100.0 - percentileofscore(benchmarks['rmsd'],RMSD_noRC)

        if all_same(FIRST_noRC) or all_same(RCI_noRC):
            spearman_noRC = np.nan
            corr_score = np.nan
            check_quiet_print(quiet,'*WARNING Spearman correlation coefficient cannot be determined, setting correlation score to NaN * ',end='')
        else:
            spearman_noRC = spearmanr(RCI_noRC,FIRST_noRC)[0]
            corr_score = percentileofscore(benchmarks['corr'],spearman_noRC)

        if cyrange:
            RCI_noRC_welldefined = [x for i,x in enumerate(RCI['score']) if i not in FIRST_nan and not np.isnan(x) and RCI['shift_types'][i] != 'RC' and FIRST['welldefined'][i]=='1']
            FIRST_noRC_welldefined = [x for i,x in enumerate(FIRST['score']) if i not in RCI_nan and not np.isnan(x) and RCI['shift_types'][i] != 'RC' and FIRST['welldefined'][i]=='1']
            
            RMSD_noRC_welldefined =  calc_RMSD(RCI_noRC_welldefined,FIRST_noRC_welldefined)
            RMSD_score_welldefined = 100.0 - percentileofscore(benchmarks['rmsd_welldefined'],RMSD_noRC_welldefined)
            #RMSD_score_welldefined = 100.0 - percentileofscore(benchmarks['rmsd_welldefined'],RMSD_noRC_welldefined)

            if all_same(FIRST_noRC_welldefined) or all_same(RCI_noRC_welldefined):
                spearman_noRC_welldefined = np.nan
                corr_score_welldefined = np.nan
                check_quiet_print(quiet,'*WARNING Spearman correlation coefficient for well-defined residues cannot be determined, setting correlation score to NaN * ',end='')
            else:
                spearman_noRC_welldefined = spearmanr(RCI_noRC_welldefined,FIRST_noRC_welldefined)[0]
                corr_score_welldefined = percentileofscore(benchmarks['corr_welldefined'],spearman_noRC_welldefined)

        av_perc_shifts = round(np.nanmean([RCI['shifts'][i[0]] for i in enumerate(RCI['resi']) if not np.isnan(i[1])])*100.0,1)

        if av_perc_shifts < 75:
            av_perc_shifts_out = str(av_perc_shifts)+' (RCI less reliable!)'
            check_quiet_print(quiet,'*WARNING chemical shift completeness (' + str(av_perc_shifts) +'%)' +' is below recommended minimum (75%), RCI may be unreliable* DONE')
        else:
            av_perc_shifts_out = str(av_perc_shifts)
            check_quiet_print(quiet,'DONE')
            
        # write output file
        if not os.path.exists(output_dir+'/ANSURR_output/'):
            os.makedirs(output_dir+'/ANSURR_output/figs')
            os.makedirs(output_dir+'/ANSURR_output/out')

        if pdb_chain not in nef_output['per_residue']:
            nef_output['per_residue'][pdb_chain] = {}
        if rci_set not in nef_output['per_residue'][pdb_chain]:
            nef_output['per_residue'][pdb_chain][rci_set] = {}
        if rci_chain not in nef_output['per_residue'][pdb_chain][rci_set] :
            nef_output['per_residue'][pdb_chain][rci_set][rci_chain] = {}
        if MODEL_ID not in nef_output['per_residue'][pdb_chain][rci_set][rci_chain]:
            nef_output['per_residue'][pdb_chain][rci_set][rci_chain][MODEL_ID] = {}

        RCI_FIRST_out = open(output_dir+'/ANSURR_output/out/'+PDB_ID+'_'+SHIFT_ID+'.out','w')
        for i in enumerate(FIRST['resi']):

            ss = str(FIRST['ss'][i[0]])
            wd = str(FIRST['welldefined'][i[0]])
            if ss == 'nan':
                ss = '.'
            if wd == 'nan':
                wd = '.'

            if FIRST['resn'][i[0]] != 'XXX':
                resi = FIRST['resi'][i[0]]
                resn = FIRST['resn'][i[0]]
                rci = float(RCI['score'][i[0]])
                rigidity = float(FIRST['score'][i[0]])
                shifts_frac = RCI['shifts'][i[0]]
                shift_types = RCI['shift_types'][i[0]]

                RCI_FIRST_out.write('{:>5}'.format(resi)+' '+'{:<4}'.format(resn)+'{:<23.20f}'.format(rci)+'{:<23.20f}'.format(rigidity)+ss+' '+wd+' '+'{:<5.2f}'.format(shifts_frac)+'{:<10}'.format(shift_types)+'\n')
                
                rci = '{0:.8f}'.format(rci) if str(rci) != 'nan' else '.'
                rigidity = '{0:.8f}'.format(rigidity) if str(rigidity) != 'nan' else '.'
                shift_types = shift_types if str(shift_types) != 'nan' else '.'

                #### NEF
                if cyrange:
                    nef_output['per_residue'][pdb_chain][rci_set][rci_chain][MODEL_ID][resi] = {'resn':resn,'rci':rci,'rigidity':rigidity,'ss':ss,'wd':bool(int(wd)),'shift_types':shift_types,'completeness':shifts_frac}
                else:
                    nef_output['per_residue'][pdb_chain][rci_set][rci_chain][MODEL_ID][resi] = {'resn':resn,'rci':rci,'rigidity':rigidity,'ss':ss,'wd':'.','shift_types':shift_types,'completeness':shifts_frac}

        RCI_FIRST_out.close()

        # append to scores.out
        scores = open(output_dir+'/ANSURR_output/scores.out','a+')
        if cyrange:
            scores.write('PDB: '+ '{:<11}'.format(PDB_ID) + ' SHIFTS: '+ '{:<11}'.format(SHIFT_ID) + ' SHIFT%: '+ '{:3}'.format(av_perc_shifts) + ' Corr: '+ '{: 5.3f}'.format(spearman_noRC) + ' CorrScore: '+ '{:4.1f}'.format(round(corr_score,1)) + ' RMSD: '+ '{:4.3f}'.format(RMSD_noRC) + ' RMSDScore: '+ '{:4.1f}'.format(round(RMSD_score,1))+ ' Corr-WD: '+ '{: 5.3f}'.format(spearman_noRC_welldefined)+' CorrScore-WD: '+ '{:4.1f}'.format(round(corr_score_welldefined,1)) + ' RMSD-WD: '+ '{:4.3f}'.format(RMSD_noRC_welldefined) + ' RMSDScore-WD: '+ '{:4.1f}'.format(round(RMSD_score_welldefined,1))+'\n')
        else:
            scores.write('PDB: '+ '{:<11}'.format(PDB_ID) + ' SHIFTS: '+ '{:<11}'.format(SHIFT_ID) + ' SHIFT%: '+ '{:3}'.format(av_perc_shifts) + ' Corr: '+ '{: 5.3f}'.format(spearman_noRC) + ' CorrScore: '+ '{:4.1f}'.format(round(corr_score,1)) + ' RMSD: '+ '{:4.3f}'.format(RMSD_noRC)+ ' RMSDScore: '+ '{:4.1f}'.format(round(RMSD_score,1))+ '\n')
        scores.close()

        if pdb_chain not in nef_output['per_model']:
            nef_output['per_model'][pdb_chain] = {}
        if rci_set not in nef_output['per_model'][pdb_chain]:
            nef_output['per_model'][pdb_chain][rci_set] = {}
        if rci_chain not in nef_output['per_model'][pdb_chain][rci_set] :
            nef_output['per_model'][pdb_chain][rci_set][rci_chain] = {}

        if cyrange:
            nef_output['per_model'][pdb_chain][rci_set][rci_chain][MODEL_ID] = {'shift_completeness':'{0:.1f}'.format(av_perc_shifts),'corr':'{0:.4f}'.format(spearman_noRC),'rmsd':'{0:.4f}'.format(RMSD_noRC),'corr_score':'{0:.1f}'.format(corr_score),'rmsd_score':'{0:.1f}'.format(RMSD_score),'corr_wd':'{0:.4f}'.format(spearman_noRC_welldefined),'rmsd_wd':'{0:.4f}'.format(RMSD_noRC_welldefined),'corr_score_wd':'{0:.1f}'.format(corr_score_welldefined),'rmsd_score_wd':'{0:.1f}'.format(RMSD_score_welldefined)}
        else:
            nef_output['per_model'][pdb_chain][rci_set][rci_chain][MODEL_ID] = {'shift_completeness':'{0:.1f}'.format(av_perc_shifts),'corr':'{0:.4f}'.format(spearman_noRC),'rmsd':'{0:.4f}'.format(RMSD_noRC),'corr_score':'{0:.1f}'.format(corr_score),'rmsd_score':'{0:.1f}'.format(RMSD_score),'corr_wd':'.','rmsd_wd':'.','corr_score_wd':'.','rmsd_score_wd':'.'}

        # plot figs
        #plot(RCI,FIRST)

        plt.figure()

        helices_x = [FIRST['resi'][i[0]] for i in enumerate(FIRST['ss']) if i[1] == 'H'] # FIRST resi can be nan!
        helices_y = [1.04]*len(helices_x)
        strands_x = [FIRST['resi'][i[0]] for i in enumerate(FIRST['ss']) if i[1] == 'E']
        strands_y = [1.04]*len(strands_x)

        RC_residues_x = [FIRST['resi'][i[0]] for i in enumerate(FIRST['resi']) if RCI['shift_types'][i[0]] == 'RC']
        if len(RC_residues_x) > 0:
            RC_residues_y = [1.02]*len(RC_residues_x)
            plt.scatter(RC_residues_x,RC_residues_y,color='black',marker="x",s=5)

        if cyrange:
            welldefined_x = [FIRST['resi'][i[0]] for i in enumerate(FIRST['welldefined']) if i[1] == '1'] # FIRST resi can be nan!
            welldefined_y = [1.06]*len(welldefined_x)
            illdefined_x = [FIRST['resi'][i[0]] for i in enumerate(FIRST['welldefined']) if i[1] == '0'] # FIRST resi can be nan!
            illdefined_y = [1.06]*len(illdefined_x)
            plt.scatter(welldefined_x,welldefined_y,color='limegreen',s=2)
            plt.scatter(illdefined_x,illdefined_y,color='grey',s=2)

        # get range over which to plot
        resi_range = []
        for i in enumerate(FIRST['resi']):
            if str(FIRST['score'][i[0]]) != 'nan' and str(RCI['score'][i[0]]) != 'nan':
                resi_range.append(i[1])

        plt.plot(RCI['resi'],RCI['score'])
        plt.plot(FIRST['resi'],FIRST['score'])
        plt.scatter(helices_x,helices_y,color='red',s=2)
        plt.scatter(strands_x,strands_y,color='blue',s=2)
        plt.xlabel('residue number',size=10)
        plt.ylabel('flexibility',size=10)
        plt.ylim(-0.025,1.075)
        plt.xlim(min(resi_range)-5,max(resi_range)+5)
        if cyrange:
            plt.title('structure: ' + PDB_ID + ' shifts: ' + SHIFT_ID + ' shift%: '+av_perc_shifts_out+'\n All residues | Corr: '+"{:.3f}".format(spearman_noRC)+' CorrScore: ' + "{:.1f}".format(corr_score)+' RMSD: '+"{:.3f}".format(RMSD_noRC)+' RMSDScore: ' + "{:.1f}".format(RMSD_score)+'\n Well-defined | Corr: '+"{:.3f}".format(spearman_noRC_welldefined)+' CorrScore: ' + "{:.1f}".format(corr_score_welldefined)+' RMSD: '+"{:.3f}".format(RMSD_noRC_welldefined)+' RMSDScore: ' + "{:.1f}".format(RMSD_score_welldefined),fontsize=10)
        else:
            plt.title('structure: ' + PDB_ID + ' shifts: ' + SHIFT_ID + ' shift%: '+av_perc_shifts_out+'\n Corr: '+"{:.3f}".format(spearman_noRC)+' CorrScore: ' + "{:.1f}".format(corr_score)+' RMSD: '+"{:.3f}".format(RMSD_noRC)+' RMSDScore: ' + "{:.1f}".format(RMSD_score),fontsize=10)
        plt.savefig(output_dir+'/ANSURR_output/figs/'+PDB_ID+'_'+SHIFT_ID+'.png',dpi=300,bbox_inches='tight')

        plt.close()

        return nef_output
    else:
        return nef_output

def main():

    #compare_rci_rigidity(sys.argv[1],sys.argv[2])
    sys.exit(0)
   
if __name__ == "__main__":
    main()
