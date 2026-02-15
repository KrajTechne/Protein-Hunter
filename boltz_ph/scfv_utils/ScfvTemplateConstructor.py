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
# Class Imports
from AnalyzeFAB import AnalyzeFAB
from mpnn_scorer import MPNNScorer

class ScfvTemplateConstructor:
    def __init__(self, struc_fab_path: str, struc_fab_target_path: str, scheme: str = 'martin'):
        self.struc_fab_path = struc_fab_path
        self.struc_fab_target_path = struc_fab_target_path
        self.atom_lines = self._read_pdb_file()
        self.fab_analysis = AnalyzeFAB(scheme = scheme)
        self.fixed_residues = self.fab_analysis.analyze_fab(struc_fab_path, struc_fab_target_path)
        self.fab_dict = self.fab_analysis.fab_dict
    
    def _read_pdb_file(self):
        """ Reads the PDB file and returns the atom lines """
        with open(self.struc_fab_path, 'r') as f:
            lines = f.readlines()
        
        atom_lines = [line for line in lines if line.startswith("ATOM") or line.startswith("HETATM")]
        return atom_lines
    
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
    
    def extract_fab_mpnn_scores(self, fab_seq: str, model_type:str = "soluble_mpnn"):
        """ Extracts MPNN scores for each residue in the Fab sequence """
        
        mpnn_scorer = MPNNScorer()
        
        def mpnn_score(model_type:str = "soluble_mpnn"):
            path_output_folder = "/Volumes/sandbox/denovotrial/mpnn_scoring/cd3_fab_target/"
            results = mpnn_scorer.score(pdb_path=self.struc_fab_target_path,
                                out_folder= path_output_folder,
                                model_type= model_type,
                                chains_to_score=["A", "B"],
                                single_aa_score = True,
                                use_sequence = True,
                                fixed_residues= None,  
                                use_atom_context=False, 
                                verbose=False,
                                batch_size = 1,
                                number_of_batches = 10)
            return results
        
        scores = mpnn_score(model_type= model_type)
        scores_fab = mpnn_scorer.extract_probs_fab(fab_seq = fab_seq, scores = scores)
        return scores_fab
    
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
    
    def create_final_fixed_designable_dict(self, cys_indices_one_pos: list,  linker_length: int, mpnn_probs: list, prob_fixed_threshold: float = 0.25):
        """ Creates final fixed and designable dict for Protein Hunter input based on mpnn_probs & upper-core mask """
    
        # Extract information from the fab_analysis object
        residue_mask_dictionary = self.fab_dict['residue_masks']
        annotated_fab = self.fab_dict['annotated_fab']
        mask_upper_core = residue_mask_dictionary['upper_core']
        orientation = annotated_fab['orientation']
        seq_len_heavy = len(annotated_fab['heavy']['seq'])
        seq_len_light = len(annotated_fab['light']['seq'])
    
        # Define dictionaries to store fixed and designable residues indices
        final_fixed_designable_dict = {'fixed': [], 'designable': []}
        res_designable_dict = {}

        # Add Linker residues as designable residues
        # Linker exists between heavy and light chains but need to know order
        if orientation == "VH-VL":
            linker_start = seq_len_heavy
        elif orientation == "VL-VH":
            linker_start = seq_len_light
        else:
            raise ValueError("Invalid orientation. Must be either 'VH-VL' or 'VL-VH.")

        # 1. Create the MPNN with Linker and Mask with Linker arrays
        mask_linker = np.zeros(linker_length)
        mask_upper_core_linker = np.concatenate((
            mask_upper_core[:linker_start], 
            mask_linker, 
            mask_upper_core[linker_start:]
        ))
        mpnn_probs_linker = np.concatenate((
            mpnn_probs[:linker_start], 
            mask_linker, 
            mpnn_probs[linker_start:]
        ))

        # 2. Loop through each residue in the fab and check if it is in the upper core
        # If in upper core, it is a fixed residue
        # If not in upper core and mpnn_prob > prob_fixed_threshold, it is fixed
        # If not in upper core and mpnn_prob < prob_fixed_threshold, it is designable
        for index, val in enumerate(mask_upper_core_linker):
            res_1pos = index + 1
            # If residue is in upper core, keep fixed
            if val:
                design_cat = 'fixed'
            # If residue is not in upper core, check if mpnn_probs_linker is above threshold
            else:
                # if mpnn prob of current amino acid is above threshold, keep fixed
                if mpnn_probs_linker[index] > prob_fixed_threshold:
                    design_cat = 'fixed'
                else:
                    design_cat = 'designable'
            
            # If residue is a cysteine in framework region, define class as 'designable'
            if res_1pos in cys_indices_one_pos:
                design_cat = 'designable'

            # Add residue to final_fixed_designable_dict
            res_designable_dict[res_1pos] = design_cat
            final_fixed_designable_dict[design_cat].append(res_1pos)

        return final_fixed_designable_dict, res_designable_dict
    
    def create_sc_pdb_file(self, output_file_path, res_designable_dict, linker_length=20):
        """ Writes a single chain of a PDB file to a new file, 
            with an optional linker of specified linker length added between first chain's last residue and second chain's first residue."""

        pdb_chain_index = 21
        pdb_residue_index = slice(22, 26) # PDB residue number columns. Gets all residue numbers including those >= 1000, but less than 10000
        chain_A_lines = [lines for lines in self.atom_lines if lines.startswith('ATOM') and lines[pdb_chain_index] == 'A']
        chain_B_lines = [lines for lines in self.atom_lines if lines.startswith('ATOM') and lines[pdb_chain_index] == 'B']
        final_chain_A_residue_index = int(chain_A_lines[-1][pdb_residue_index].strip())

        # Define linker residue numbering (Cleaner to define start and end separately)
        linker_start = final_chain_A_residue_index + 1
        linker_end = final_chain_A_residue_index + linker_length
        chain_B_offset = linker_end

        new_chain_B_starting_residue_index = final_chain_A_residue_index + linker_length # Want to add in linker of specified size
        # Challenge: How to ensure groups of lines associated with same residue are modified together?
        # One way: Use a dictionary to map residue numbers to their corresponding lines
        residue_to_lines = {}
        # Chain A: Extract lines for each residue in chain A. Preserving Chain A numbering
        for line in chain_A_lines:
            residue_number = int(line[pdb_residue_index].strip())
            residue_line_list = residue_to_lines.get(residue_number, [])
            residue_line_list.append(line)
            residue_to_lines[residue_number] = residue_line_list
        
        # Linker: Create fake linker residues (no atom lines, just for numbering)
        for residue_number in range(linker_start, linker_end + 1): # +1 to include end residue
            residue_to_lines[residue_number] = [] # No atom lines for linker residues
        
        # Chain B: Extract lines for each residue in chain B, but adjusting residue numbering to account for linker
        for line in chain_B_lines:
            residue_number = int(line[pdb_residue_index].strip())
            chain_linker_adjusted_residue_number = residue_number + chain_B_offset
            residue_line_list = residue_to_lines.get(chain_linker_adjusted_residue_number, [])
            residue_line_list.append(line)
            residue_to_lines[chain_linker_adjusted_residue_number] = residue_line_list


        # Add in lines based on whether designable or fixed. Output in single chain and with second chain's residue numbering adjusted by linker_length
        # Create seq input string for Protein Hunter
        pdb_lines = []
        seq_input = ""
        for residue_number, lines in residue_to_lines.items():
            designability_attribute = res_designable_dict[residue_number]
            
            # Preserve residue lines as is, but adjust chain to 'A' and residue numbering (regardless of fixed or not)
            for line in lines:
                new_line = line[:pdb_chain_index] + 'A' + line[pdb_chain_index+1:pdb_residue_index.start] + str(residue_number).rjust(4) + line[pdb_residue_index.stop:]
                pdb_lines.append(new_line)
            
            if designability_attribute == 'fixed':
                # Extract amino acid for seq input
                line = lines[0] # Take first line of residue to extract amino acid
                aa3 = line[17:20].strip() # Extract 3-letter amino acid code
                aa1 = seq.ProteinSequence().convert_letter_3to1(aa3)
                seq_input += aa1
            else: # designable, skip residue lines and add 'X' to seq input as model will redesign them
                seq_input += "X"

        # Write to output PDB file       
        with open(output_file_path, 'w') as file:
            for line in pdb_lines:
                file.write(line)
            file.write('END\n')
        print(f"Written to {output_file_path}")

        return seq_input
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
    
    def create_protein_hunter_inputs(self, output_file_path, linker_length=20):
        """ Creates the single chain PDB file and extract seq input for Protein Hunter inputs """
        
        # 1. Annotate fab
        seq_dict, annotated_paired_seq = self.annotate_fab()
        paired_seq = annotated_paired_seq['seq']
        
        # 2. Extract fab mpnn scores
        mpnn_scores_fab = self.extract_fab_mpnn_scores(fab_seq= paired_seq, model_type= "soluble_mpnn")
        # 3. Create mapping of fixed/designable/variable indices
        paratope_indices, cys_indices = self.create_fixed_designable_variable_array(annotated_paired_seq,
                                                                                    linker_length=linker_length, 
                                                                                    cdr_extend=2)
        # 4. Create Final Fixed and Designable Dict 
        final_fixed_designable_dict, sc_res_designable_dict = self.create_final_fixed_designable_dict(cys_indices_one_pos= cys_indices,
                                                                                                      linker_length= linker_length,
                                                                                                      mpnn_probs= mpnn_scores_fab)
        # 5. Create single chain PDB file with linker
        seq_input = self.create_sc_pdb_file(output_file_path=output_file_path, res_designable_dict= sc_res_designable_dict, linker_length=linker_length)
        
        # 6. Convert to mmCIF format
        output_cif_path = self.convert_pdb_to_cif(input_pdb_path=output_file_path)

        # 7. Create JSON file for per-residue bias. Need to mutate cysteine residues to common hydrophobic residues
        json_path = self.create_mpnn_per_residue_bias(cys_residue_numbers=cys_indices)

        print("Number of Fixed Residues:", len(final_fixed_designable_dict['fixed']))
        print("Number of Designable Residues:", len(final_fixed_designable_dict['designable']))
        print("✅ Created Protein Hunter Inputs")
        return seq_dict, paired_seq, seq_input, output_cif_path, final_fixed_designable_dict, paratope_indices, json_path