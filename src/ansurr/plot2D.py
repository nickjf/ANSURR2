import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import numpy as np

from ansurr import adjust_text as at

def plot_scores(output_dir='',cyrange=False):
	data = {}
	for line in open(output_dir+'scores.out','r'):
		if line[:4] == 'PDB:':
			line = line.split()
			pdb = line[1][:[pos for pos, char in enumerate(line[1]) if char == '_'][-1]] # accounts for "_" appearing in structure name
			model = int(line[1].split('_')[-1])
			shifts = line[3]
			shift_perc = float(line[5])
			corr_score = float(line[9])
			rmsd_score = float(line[13])
			if cyrange:
				corr_score_wd = float(line[17])
				rmsd_score_wd = float(line[21])
		
			if line[9] != 'nan':
				if pdb not in data:
					data[pdb] = {}
		
				if shifts not in data[pdb]:
					if cyrange:
						if str(corr_score_wd) != 'nan':
							data[pdb][shifts] = {'model':[model],'rmsd':[rmsd_score],'corr':[corr_score],'rmsd_welldefined':[rmsd_score_wd],'corr_welldefined':[corr_score_wd],'shift_perc':[shift_perc]}
					else:
						if str(corr_score) != 'nan':
							data[pdb][shifts] = {'model':[model],'rmsd':[rmsd_score],'corr':[corr_score],'shift_perc':[shift_perc]}

				else:
					if str(corr_score) != 'nan':
						data[pdb][shifts]['model'].append(model)
						data[pdb][shifts]['rmsd'].append(rmsd_score)
						data[pdb][shifts]['corr'].append(corr_score)
						data[pdb][shifts]['shift_perc'].append(shift_perc)
					if cyrange:
						if str(corr_score_wd) != 'nan':
							data[pdb][shifts]['rmsd_welldefined'].append(rmsd_score_wd)
							data[pdb][shifts]['corr_welldefined'].append(corr_score_wd)

				
	for pdb in data:
		for shifts in data[pdb]:
			plt.cla()
			plt.scatter(data[pdb][shifts]['rmsd'], data[pdb][shifts]['corr'], c='black',s=8)

			labels = [plt.text(data[pdb][shifts]['rmsd'][i], data[pdb][shifts]['corr'][i], data[pdb][shifts]['model'][i],
							   ha='center', va='center', color='blue',size=6) for i in range(len(data[pdb][shifts]['rmsd']))]

			plt.ylabel('correlation score',size=12)
			plt.xlabel('RMSD score',size=12)
			plt.axis('scaled')  
			plt.xlim(-5,105)
			plt.ylim(-5,105)
			shift_perc = int(round(np.nanmean(data[pdb][shifts]['shift_perc']),0)) # just in case models don't have the same shift completeness - very unlikely
			if shift_perc < 75:
				shift_perc_out = str(shift_perc) + ' (unreliable!)'
			else:
				shift_perc_out = str(shift_perc)

			plt.title('Structure: '+pdb+'\nShifts: '+shifts+' Shift%: '+shift_perc_out)
			at.adjust_text(labels,expand_text=(1.75, 1.75), expand_points=(1.75, 1.75),arrowprops=dict(arrowstyle='->', color='red', alpha=0.6, linewidth=0.5))
			plt.tight_layout()
			plt.savefig(output_dir+pdb+'_'+shifts+'.png',dpi=300)


			if cyrange:
				plt.cla()
				plt.scatter(data[pdb][shifts]['rmsd_welldefined'], data[pdb][shifts]['corr_welldefined'], c='black',s=8)
				labels_cyrange = [plt.text(data[pdb][shifts]['rmsd_welldefined'][i], data[pdb][shifts]['corr_welldefined'][i], data[pdb][shifts]['model'][i],
							   ha='center', va='center', color='blue',size=6) for i in range(len(data[pdb][shifts]['rmsd_welldefined']))]

				plt.ylabel('correlation score',size=12)
				plt.xlabel('RMSD score',size=12)
				plt.axis('scaled')  
				plt.xlim(-5,105)
				plt.ylim(-5,105)
				shift_perc = int(round(np.nanmean(data[pdb][shifts]['shift_perc']),0)) # have sep shift completness for well defined region?
				if shift_perc < 75:
					shift_perc_out = str(shift_perc) + ' (unreliable!)'
				else:
					shift_perc_out = str(shift_perc)

				plt.title('Structure: '+pdb+' (WD)\nShifts: '+shifts+' Shift%: '+shift_perc_out)
				at.adjust_text(labels_cyrange,expand_text=(1.75, 1.75), expand_points=(1.75, 1.75),arrowprops=dict(arrowstyle='->', color='red', alpha=0.6, linewidth=0.5))
				plt.tight_layout()
				plt.savefig(output_dir+pdb+'_'+shifts+'_welldefined.png',dpi=300)

def main():

    plot_scores(sys.argv[1])
    sys.exit(0)
   
if __name__ == "__main__":
    main()