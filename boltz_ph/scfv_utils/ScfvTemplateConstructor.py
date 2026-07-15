import os
import numpy as np
import biotite
import gemmi
import json
import biotite.sequence as seq
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
import py2Dmol
import sys
# Class Imports
from AnalyzeFAB import AnalyzeFAB
from mpnn_scorer import MPNNScorer
from StrucTools import *

class ScfvTemplateConstructor:
    def __init__(self, struc_fab_path: str, struc_fab_target_path: str, scheme: str = 'martin', verbose = False, threshold = 4):
        self.struc_fab_path = struc_fab_path
        self.struc_fab_target_path = struc_fab_target_path
        self.structure = extract_atom_array(struc_file_path = struc_fab_target_path, ca_only = False)
        self.fab_analysis = AnalyzeFAB(scheme = scheme, verbose = verbose)
        self.fixed_residues = self.fab_analysis.analyze_fab(struc_fab_path, struc_fab_target_path, threshold = threshold)
        self.fab_dict = self.fab_analysis.fab_dict
    
    def visualize_structure(self, pdb_file_path):
        """ Visualizes the PDB structure using Sokrypton py2Dmol """
        
        viewer = py2Dmol.view()
        viewer.add_pdb(pdb_file_path)
        return viewer.show()
    
    def annotate_fab(self):
        """ Annotates fab using the AnalyzeFAB Class """
        
        annotated_paired_seq = self.fab_dict['annotated_fab']
        seq_dict = {'heavy': annotated_paired_seq['heavy']['seq'], 'light': annotated_paired_seq['light']['seq']}

        return seq_dict, annotated_paired_seq
    
    def create_fixed_designable_variable_array(self, annotated_paired_seq: dict, linker_length: int = 20, cdr_extend: int = 2):
        """ Create a list of the linker length with fixed, variable, designable residues
            class_mapping = {'variable': 0, 'fixed': 1, 'designable': 2}
            Approach:
                1. Define list of 0s of desired scfv length
                2. Fill in at specified positions (designable or fixed) with the appropriate value: 1 or 2
                3. Remaining values are by default variable residues and already have value of 0
                4. Return np array of class_mapping values
            Args:
                annotated_paired_seq (dict): Dictionary of annotated paired sequence
                linker_length (int): Length of linker
                cdr_extend (int): Number of residues to extend CDR Definition by on both sides of CDR region
                risk (bool): Allows for mutation of potential Vernier Zone Residues
            Returns:
                scfv_mapping (np.ndarray): Numpy array of class_mapping values
                paratope_residues_str (str): String of 1-indexed residue positions in paratope separated by commas
                cys_indices_one_indexed (list): List of 1-indexed residue positions of cysteine residues in antibody framework regions 
        """

        if annotated_paired_seq['orientation'] == "VH-VL":
            first_chain = "heavy"
            sec_chain = "light"
        else:
            first_chain = "light"
            sec_chain = "heavy"
        
        annotated_chains = {'heavy': annotated_paired_seq['heavy'], 'light': annotated_paired_seq['light']}
        first_chain_length = len(annotated_chains[first_chain]['seq'])
        sec_chain_length = len(annotated_chains[sec_chain]['seq'])
        scfv_length = first_chain_length + sec_chain_length + linker_length
        
        cys_indices_full = []
        paratope_residues_zero_indexed = np.zeros(scfv_length, dtype=int)
        for chain, chain_dict in annotated_chains.items():
    
            # Determine residue offset based on chain order
            if chain == first_chain:
                residue_offset = 0  
            else:
                residue_offset = first_chain_length + linker_length  # Account for linker length and first chain length
    
            # For each region in the chain
            for region, region_index_dict in chain_dict['region_loc_dict'].items():
        
                # Get adjusted start and end indices for the region. Adjustment accounts for linker and position of heavy/light in scfv
                start_adj = region_index_dict['start'] + residue_offset
                end_adj = region_index_dict['end'] + residue_offset 
    
                # Assign to appropriate category in indices_dict
                if "fmwk" in region:
                    # Check whether Cysteine is in fmwk region -> must mutate  out -> index should be part of designable residues
                    seq_fmwk = chain_dict['region_seqs_dict'][region]
                    # 1. Find all indices where 'C' appears
                    cys_indices_list = [i + start_adj for i, char in enumerate(seq_fmwk) if char == "C"]
                    cys_indices_full.extend(cys_indices_list)
                    
                elif "cdr" in region:
                    
                    # Extend CDR definition by cdr_extend on both sides (Accounts for discrepancies between antibody annotation schemes on CDR definitions)
                    safe_start = max(start_adj - cdr_extend, 0) # Ensure start does not go below 0
                    safe_end = min(end_adj + cdr_extend + 1, scfv_length) # Ensure end doesn't exceed scfv length. + 1 to include end index 
                    
                    # Define paratope based on the extended start and end indices 
                    paratope_residues_zero_indexed[safe_start: safe_end] = 1
                        

        # Convert Paratope Residues Indices to a 1-indexed string format for Protein Hunter input (separated by commas)
        paratope_indices_one_indexed = [index + 1 for index, val in enumerate(paratope_residues_zero_indexed) if val == 1]
        #print("Paratope Indices 1-indexed: ", paratope_indices_one_indexed)
        paratope_indices_str = ','.join([str(i) for i in paratope_indices_one_indexed])
        #print("Paratope Indices String: ", paratope_indices_str)
        # Convert Cys indices to 1-indexed list
        cys_indices_one_indexed = [val + 1 for val in cys_indices_full]

        return paratope_indices_str, cys_indices_one_indexed
    
    def create_final_fixed_designable_dict(self, cys_indices_one_pos: list,  linker_length: int):
        """ Creates final fixed and designable dict for Protein Hunter input based on mpnn_probs & upper-core mask """
    
        # Extract information from the fab_analysis object
        residue_mask_dictionary = self.fab_dict['residue_masks']
        annotated_fab = self.fab_dict['annotated_fab']
        mask_fixed = residue_mask_dictionary['fixed']
        orientation = annotated_fab['orientation']
        seq_len_heavy = len(annotated_fab['heavy']['seq'])
        seq_len_light = len(annotated_fab['light']['seq'])

        # Add Linker residues as designable residues
        # Linker exists between heavy and light chains but need to know order
        if orientation == "VH-VL":
            linker_start = seq_len_heavy
            seq_input = annotated_fab['heavy']['seq'] + ("X" * linker_length) + annotated_fab['light']['seq']
        elif orientation == "VL-VH":
            linker_start = seq_len_light
            seq_input = annotated_fab['light']['seq'] + ("X" * linker_length) + annotated_fab['heavy']['seq']
        else:
            raise ValueError("Invalid orientation. Must be either 'VH-VL' or 'VL-VH.")
        
        # 0. Create list of linker 1-indexed res pos to ensure always designable when sampling for additional fixed residues
        linker_1_pos = [index + 1 for index, val in enumerate(seq_input) if val == "X"]
        # 1. Create the MPNN with Linker and Mask with Linker arrays
        mask_linker = np.zeros(linker_length)
        mask_fixed_linker = np.concatenate((
            mask_fixed[:linker_start], 
            mask_linker, 
            mask_fixed[linker_start:]
        ))
        

        # 2. Loop through each residue in the fab and check if it is in the set of fixed residues
        # If it is a fixed residue, it's design class is fixed. If not, design class is designable
        
        # Define dictionaries to store fixed and designable residues indices.
        final_fixed_designable_dict = {'fixed': [], 'designable': [], 'linker': linker_1_pos}
        
        for index, val in enumerate(mask_fixed_linker):
            res_1pos = index + 1
            
            # If residue is not in set of fixed residues, check if mpnn_probs_linker is above threshold
            if not val:
                design_cat = 'designable'
            
            # If residue is in set of fixed residues, it is fixed
            elif val:
                design_cat = 'fixed'
                
            # If residue is a cysteine in framework region, define class as 'designable'
            # Possible for cysteines in framework to be in fixed residues, but cannot allow them to remain fixed
            if res_1pos in cys_indices_one_pos:
                design_cat = 'designable'
                mask_fixed_linker[index] = False # Guarantee that cysteines in framework are designable in the fixed mask aligned to user-specified linker length
            
            # Fix the first residue
            if res_1pos == 1:
                design_cat = 'fixed'

            # Add residue to final_fixed_designable_dict
            final_fixed_designable_dict[design_cat].append(res_1pos)
        
        self.fab_dict['residue_masks']['fixed_with_linker'] = mask_fixed_linker

        return final_fixed_designable_dict, seq_input
    
    def create_sc_pdb_file(self, pdb_save_path, linker_length=20):
        """ 
        Writes a single chain of a PDB file using Biotite AtomArrays.
        Compatible with both PDB and CIF inputs.
        """
        # Create a working copy to prevent mutating the stored structure
        structure = self.structure.copy()

        num_chains = len(struc.get_chains(structure))
        
        # 1. Separate Chains A and B
        if num_chains == 2:
            chain_a = structure[structure.chain_id == 'A']
            chain_b = structure[structure.chain_id == 'B']
        elif num_chains == 3:
            chain_a = structure[structure.chain_id == 'A']
            chain_b = structure[structure.chain_id == 'B']
            chain_c = structure[structure.chain_id == 'C']
        else:
            raise ValueError("Input structure must contain chains 'A' and 'B'.")

        # 2. Determine offset
        # Get the highest residue number in Chain A
        last_residue_a = np.max(chain_a.res_id)
        
        # Calculate the shift required for Chain B
        # The offset ensures B starts after A + Linker
        # Original logic: new_B_res_num = old_B_res_num + (End_A + Linker)
        # Note: This adds the offset to existing numbers. If B starts at 1, it works. 
        offset = last_residue_a + linker_length
        
        # 3. Apply transformation to Chain B
        chain_b.res_id += offset
        chain_b.chain_id[:] = 'A' # Rename chain B to A

        if num_chains == 3:
            chain_c.chain_id[:] = 'B' # Rename chain C to B
        
        # 4. Merge chains
        # Biotite automatically handles atom concatenation
        if num_chains == 2:
            combined_structure = chain_a + chain_b
        elif num_chains == 3:
            combined_structure = chain_a + chain_b + chain_c
        
        # 5. Write to PDB
        # Biotite PDBFile handles the formatting (columns, spacing) automatically
        pdb_file = pdb.PDBFile()
        pdb_file.set_structure(combined_structure)
        pdb_file.write(pdb_save_path)
        
        print(f"Written to {pdb_save_path}")

    def convert_pdb_to_cif(self, input_pdb_path):
        """Adapts the logic from Boltz's parse_pdb to convert a PDB file to mmCIF."""
        
        print(f"Reading PDB: {input_pdb_path}")
        output_cif_path = input_pdb_path.replace('.pdb', '.cif')
        
        # 1. Read the structure using Gemmi
        structure = gemmi.read_structure(input_pdb_path)
        structure.setup_entities()
        
        # 2. Apply the subchain renaming logic (Copied from Boltz source)
        # This ensures chains are correctly named for mmCIF format (e.g., handling multiple segments)
        subchain_counts, subchain_renaming = {}, {}
        for chain in structure[0]:
            subchain_counts[chain.name] = 0
            for res in chain:
                if res.subchain not in subchain_renaming:
                    subchain_renaming[res.subchain] = chain.name + str(subchain_counts[chain.name] + 1)
                    subchain_renaming[res.subchain] = str(subchain_counts[chain.name] + 1) # Simplified renaming
                    subchain_counts[chain.name] += 1
                res.subchain = subchain_renaming[res.subchain]
            
        # Update entities with new subchain names
        for entity in structure.entities:
            entity.subchains = [subchain_renaming.get(subchain, subchain) for subchain in entity.subchains]

        # 3. Create mmCIF document and write to file
        doc = structure.make_mmcif_document()
        doc.write_file(output_cif_path)
    
        print(f"✅ Converted to CIF: {output_cif_path}")
        return output_cif_path
    
    def create_mpnn_per_residue_bias(self, cys_residue_numbers: list, output_json_path: str = "/tmp/mpnn_per_residue_bias.json"):
        """ Creates the JSON file for the per-residue bias for the cysteines involved in the disulfide bond in original Antibody Framework Regions
            Cysteines typically located in positions: 23 and 87.
            Per-residue bias will ensure cysteines get mutated to common hydrophobic residues typically used for this specific scenario
            Goal: Create JSON file and save to /tmp/mpnn_per_residue_bias.json
            Args:
                cys_residue_numbers (list): List of residue numbers for cysteines involved in disulfide bond
                output_json_path (str): Path to save the JSON file
            Return:
                output_json_path (str): Path to the saved JSON file
        """
        # 1. Setup your staples (High positive bias = forced selection)
        staple_residues = ["V", "I", "L", "A", "F"]
        staple_bias = 100.0  # Strong Positive bias imposes MPNN to select between residues

        # 2. Create the dictionary
        # Format: {"H23": {"V": 100.0, "I": 100.0, ...}, "H88": {...}}
        bias_dict = {} 
        for pdb_idx in cys_residue_numbers:
            # Construct key: pdb chain + residue number, e.g., "H23", "A88", etc. (typical PDB notation)
            key = f"A{pdb_idx}"
            
            # Initialize position
            bias_dict[key] = {}
            # Add bias for each hydrophobic staple
            for aa in staple_residues:
                bias_dict[key][aa] = staple_bias

        # 3. Save to a temporary JSON file
        with open(output_json_path, "w") as f:
            json.dump(bias_dict, f, indent=2)
        print(f"Bias file created at {output_json_path} for indices: {list(bias_dict.keys())}")
        return output_json_path
    
    def create_protein_hunter_inputs(self, pdb_save_path:str, linker_length: int =20, cdr_extend: int = 2):
        """ Creates the single chain PDB file and extract seq input for Protein Hunter inputs """
        
        # 1. Annotate fab
        seq_dict, annotated_paired_seq = self.annotate_fab()
       
        # 2. Create mapping of fixed/designable/variable indices
        paratope_indices, cys_indices = self.create_fixed_designable_variable_array(annotated_paired_seq,
                                                                                    linker_length=linker_length, 
                                                                                    cdr_extend=cdr_extend)
        # 3. Create Final Fixed and Designable Indices Dict & seq_input 
        final_fixed_designable_dict, seq_input = self.create_final_fixed_designable_dict(cys_indices_one_pos= cys_indices,
                                                                                         linker_length= linker_length)
        # Add cys indices to dictionary to ensure always sampled for MPNN
        final_fixed_designable_dict['cys'] = cys_indices
        # 4. Create single chain PDB file with linker
        self.create_sc_pdb_file(pdb_save_path= pdb_save_path, linker_length=linker_length)
        
        # 5. Convert to mmCIF format
        output_cif_path = self.convert_pdb_to_cif(input_pdb_path=pdb_save_path)

        # 7. Create JSON file for per-residue bias. Need to mutate cysteine residues to common hydrophobic residues
        json_path = self.create_mpnn_per_residue_bias(cys_residue_numbers=cys_indices)

        print("Number of Fixed Residues:", len(final_fixed_designable_dict['fixed']))
        print("Number of Designable Residues:", len(final_fixed_designable_dict['designable']))
        print("✅ Created Protein Hunter Inputs")
        return seq_dict, seq_input, output_cif_path, final_fixed_designable_dict, paratope_indices, json_path