import sys

import pynmrstar

#from ansurr.functions import check_quiet_print

def export(nef_output, output_dir='',ansurr_version='.',pdb_file='.',shift_file='.',reref='.',lig='.',nonstd='.',olig='.',quiet=False):

    ansurr_per_model_lp = pynmrstar.Loop.from_scratch()
    ansurr_per_model_lp.add_tag(['_ansurr_per_model.model_code','_ansurr_per_model.pdb_chain_code','_ansurr_per_model.shifts_set_code','_ansurr_per_model.shifts_chain_code','_ansurr_per_model.backbone_shift_completeness', '_ansurr_per_model.correlation', '_ansurr_per_model.correlation_score','_ansurr_per_model.rmsd','_ansurr_per_model.rmsd_score', '_ansurr_per_model.correlation_well_defined', '_ansurr_per_model.correlation_score_well_defined','_ansurr_per_model.rmsd_well_defined','_ansurr_per_model.rmsd_score_well_defined'])

    for chain in nef_output['per_model']:
        for rci_set in nef_output['per_model'][chain]:
            for rci_chain in nef_output['per_model'][chain][rci_set]:
                for model in nef_output['per_model'][chain][rci_set][rci_chain]:
                    ansurr_per_model_lp.add_data([model, chain, rci_set, rci_chain,nef_output['per_model'][chain][rci_set][rci_chain][model]['shift_completeness'], nef_output['per_model'][chain][rci_set][rci_chain][model]['corr'],nef_output['per_model'][chain][rci_set][rci_chain][model]['corr_score'],nef_output['per_model'][chain][rci_set][rci_chain][model]['rmsd'],nef_output['per_model'][chain][rci_set][rci_chain][model]['rmsd_score'],nef_output['per_model'][chain][rci_set][rci_chain][model]['corr_wd'],nef_output['per_model'][chain][rci_set][rci_chain][model]['corr_score_wd'],nef_output['per_model'][chain][rci_set][rci_chain][model]['rmsd_wd'],nef_output['per_model'][chain][rci_set][rci_chain][model]['rmsd_score_wd']])


    #print(ansurr_per_model_lp)

    ansurr_per_residue_lp = pynmrstar.Loop.from_scratch()
    ansurr_per_residue_lp.add_tag(['_ansurr_per_residue.model_code','_ansurr_per_residue.chain_code','_ansurr_per_residue.shifts_set_code','_ansurr_per_residue.shifts_chain_code','_ansurr_per_residue.sequence_code', '_ansurr_per_residue.residue_name', '_ansurr_per_residue.rci', '_ansurr_per_residue.rigidity','_ansurr_per_residue.secondary_structure','_ansurr_per_residue.welldefined','_ansurr_per_residue.shift_types'])

    for chain in nef_output['per_residue']:
        for rci_set in nef_output['per_residue'][chain]:
            for rci_chain in nef_output['per_residue'][chain][rci_set]:
                for model in nef_output['per_residue'][chain][rci_set][rci_chain]:
                    for residue in nef_output['per_residue'][chain][rci_set][rci_chain][model]:
                        ansurr_per_residue_lp.add_data([model, chain, rci_set, rci_chain,residue, nef_output['per_residue'][chain][rci_set][rci_chain][model][residue]['resn'],nef_output['per_residue'][chain][rci_set][rci_chain][model][residue]['rci'],nef_output['per_residue'][chain][rci_set][rci_chain][model][residue]['rigidity'],nef_output['per_residue'][chain][rci_set][rci_chain][model][residue]['ss'],nef_output['per_residue'][chain][rci_set][rci_chain][model][residue]['wd'],nef_output['per_residue'][chain][rci_set][rci_chain][model][residue]['shift_types']])

    #print(ansurr_per_residue_lp)

    ansurr_sf = pynmrstar.Saveframe.from_scratch("ansurr", "ansurr")
          
            
    ansurr_sf.add_tag("sf_category", "ansurr")
    ansurr_sf.add_tag("sf_framecode", "ansurr")

    ansurr_sf.add_tag("_ansurr.ansurr_version", "ansurr")
    ansurr_sf.add_tag("_ansurr.panav_version", "ansurr")
    ansurr_sf.add_tag("_ansurr.dssp_version", "ansurr")
    ansurr_sf.add_tag("_ansurr.cyrange_version", "ansurr")

    ansurr_sf.add_tag("_ansurr.input_pdb_file", "ansurr")
    ansurr_sf.add_tag("_ansurr.input_shift_file", "ansurr")
    ansurr_sf.add_tag("_ansurr.rereference_shifts", "ansurr")
    ansurr_sf.add_tag("_ansurr.include_ligands", "ansurr")
    ansurr_sf.add_tag("_ansurr.include_nonstd_res", "ansurr")
    ansurr_sf.add_tag("_ansurr.oligomer", "ansurr")

    ansurr_sf['ansurr_version'] = ansurr_version
    ansurr_sf['panav_version'] = "2.1"
    ansurr_sf['dssp_version'] = "3.0.0"
    ansurr_sf['cyrange_version'] = "2.0"


    ansurr_sf['input_pdb_file'] = pdb_file
    ansurr_sf['input_shift_file'] = shift_file

    ansurr_sf['rereference_shifts'] = reref
    ansurr_sf['include_ligands'] = lig
    ansurr_sf['include_nonstd_res'] = nonstd
    ansurr_sf['oligomer'] = olig

    # Now add the loop we created before
    ansurr_sf.add_loop(ansurr_per_model_lp)
    ansurr_sf.add_loop(ansurr_per_residue_lp)
    
    # Now write out our saveframe to a file. Optionally specify format="json" to write in JSON format.
    ansurr_sf.write_to_file(output_dir+"_ansurr.nef")
    ansurr_sf.write_to_file(output_dir+"_ansurr.json", format_="json")


def main():
    sys.exit(0)
   
if __name__ == "__main__":
    main()