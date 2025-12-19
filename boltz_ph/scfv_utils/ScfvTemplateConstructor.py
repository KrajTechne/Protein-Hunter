from itertools import chain
import os
import numpy as np
import biotite
import gemmi
import biotite.sequence as seq
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
import py2Dmol
# Class Imports
print(os.getcwd())
from boltz_ph.scfv_utils.Scfv import Scfv

class ScfvTemplateConstructor:
    def __init__(self, fv_pdb_path: str, scheme: str = 'martin', target_name: str = 'target', orientation: str = 'VH-VL'):
        self.fv_pdb_path = fv_pdb_path
        self.atom_lines = self._read_pdb_file()
        self.scfv_annotator = Scfv(scheme = scheme)
        self.target_name = target_name
        self.orientation = orientation
    
    def _read_pdb_file(self):
        """ Reads the PDB file and returns the atom lines """
        with open(self.fv_pdb_path, 'r') as f:
            lines = f.readlines()
        
        atom_lines = [line for line in lines if line.startswith("ATOM") or line.startswith("HETATM")]
        return atom_lines
    
    def visualize_structure(self, pdb_file_path):
        """ Visualizes the PDB structure using Sokrypton py2Dmol """
        
        viewer = py2Dmol.view()
        viewer.add_pdb(pdb_file_path)
        return viewer.show()
    
    def extract_chain_sequences(self):
        """ Extract sequences for chains A and B from the PDB file """
        seq_obj = seq.ProteinSequence()
        
        # Account for possibility that input path may not be PDB. Could be PDB or mmCIF
        if self.fv_pdb_path.endswith('.pdb'):
            struc_file = pdb.PDBFile.read(self.fv_pdb_path)
            struc_atom_array = pdb.get_structure(struc_file, model = 1) # Get first model (Assumption: only one model present)
        elif self.fv_pdb_path.endswith('.cif'):
            struc_file = pdbx.CIFFile.read(self.fv_pdb_path)
            struc_atom_array = pdbx.get_structure(struc_file, model = 1) #
    
        # Extract sequences for each chain from structure atom array
        chain_ids = np.unique(struc_atom_array.chain_id)
        seq_dict = {}
        paired_seq = ""
        for chain_id in chain_ids:
            chain_atoms = struc_atom_array[struc_atom_array.chain_id == chain_id]
            aa_pos, aa3_list = struc.get_residues(chain_atoms)
            chain_aa1_seq = ''.join([seq_obj.convert_letter_3to1(aa3) for aa3 in aa3_list])
            seq_dict[chain_id] = chain_aa1_seq
            paired_seq += chain_aa1_seq
        
        return seq_dict, paired_seq
    
    def annotate_paired_sequence(self, paired_seq):
        """ Annotates the paired sequence using Scfv Annotator Class """
        
        # Add {target_id: seq} to annotator
        linker_dict = {self.target_name: ""}
        orientation_dict = {self.target_name: self.orientation}
        self.scfv_annotator.update_seqs(seq_id = self.target_name, seq = paired_seq)
        annotated_paired_seq  = self.scfv_annotator.annotate_seqs(linker_dict, orientation_dict, target_dict={}, generate_motif_commands= False)
        return annotated_paired_seq
    
    def create_fixed_designable_variable_indices(self, annotated_paired_seq: dict, linker_length: int =20, cdr_extend: int = 3):
        """ Creates fixed, designable, and variable residues indices based on annotated paired sequence """

        if self.orientation == "VH-VL":
            first_chain = "heavy"
            sec_chain = "light"
        else:
            first_chain = "light"
            sec_chain = "heavy"
        
        target_annotation = annotated_paired_seq[self.target_name]
        target_annotation_chains = {'heavy': target_annotation['heavy'], 'light': target_annotation['light']}
        first_chain_length = len(target_annotation[first_chain]['seq'])
        sec_chain_length = len(target_annotation[sec_chain]['seq'])
        print(f"Heavy Chain Regions: {target_annotation['heavy']['region_seqs_dict']}")
        print(f"Light Chain Regions: {target_annotation['light']['region_seqs_dict']}")
        #print(f"Seq Lengths: First Chain {first_chain_length}, Second Chain {sec_chain_length}, Linker Length: {linker_length}")
        #print(f"Total Length: {first_chain_length + sec_chain_length + linker_length}")

        indices_dict = {'fixed': [], 'designable': [], 'variable': []}
        paratope_indices = []

        for chain, chain_dict in target_annotation_chains.items():
    
            # Determine residue offset based on chain order
            if chain == first_chain:
                residue_offset = 1 # Residue numbering starts at 1
            else:
                residue_offset = first_chain_length + linker_length + 1 # +1 to account for 1-based indexing. Account for linker length and first chain length
    
            # For each region in the chain
            for region, region_index_dict in chain_dict['region_loc_dict'].items():
        
                # Get zero-based indices for the region
                start = region_index_dict['start']
                end = region_index_dict['end'] 
        
                # Create list of zero-based indices and adjust by residue offset to make them 1-based and account for previous chain and linker
                zero_indices_list = list(range(start, end + 1)) # +1 to include end index
                adjusted_indices_list = [i + residue_offset for i in zero_indices_list]

                # Assign to appropriate category in indices_dict
                if "fmwk" in region:    
                    indices_dict['variable'].extend(adjusted_indices_list)
                elif "cdr" in region:
                    indices_dict['fixed'].extend(adjusted_indices_list)

                    # Extract residues potentially responsible for binding (extend CDRs by cdr_extend on both sides, within chain limits)
                    # Aim: Capture potential paratope residues and pass to Protein Hunter to ensure contact check only occurs with these residues and target epitope
                    extended_start = max(start - cdr_extend, 0) # Ensure start does not go below 0
                    extended_end = min(end + cdr_extend, len(chain_dict['seq']) - 1) # Ensure end does not exceed chain length - 1
                    extended_zero_indices_list = list(range(extended_start, extended_end + 1))
                    extended_adjusted_indices_list = [i + residue_offset for i in extended_zero_indices_list]
                    paratope_indices.extend(extended_adjusted_indices_list)
        
        # For Linker region
        indices_dict['designable'] = list(range(first_chain_length + 1, first_chain_length + linker_length + 1))

        #print("Fixed indices:", indices_dict['fixed'])
        #print("Designable indices:", indices_dict['designable'])
        #print("Variable indices:", indices_dict['variable'])

        #print(f"Number of Fixed Residues: {len(indices_dict['fixed'])}")
        #print(f"Number of Designable Residues: {len(indices_dict['designable'])}")
        #print(f"Number of Variable Residues: {len(indices_dict['variable'])}")
        return indices_dict, paratope_indices
    
    def create_final_fixed_designable_dict(self, fixed_designable_variable_indices: dict, probs_dict: dict):
        """ Creates final fixed and designable dict for Protein Hunter input based on fixed, designable, and variable indices """
        
        # Establish sequence length
        seq_len = max(fixed_designable_variable_indices['fixed'] + fixed_designable_variable_indices['designable'] + fixed_designable_variable_indices['variable'])
        #print(f"Length of seq: {len(fixed_designable_variable_indices['fixed'] + fixed_designable_variable_indices['designable'] + fixed_designable_variable_indices['variable'])}")
        #print("Sequence Length:", seq_len)
        
        # Create initial dict to store residue categories: fixed and designable 
        final_fixed_designable_dict = {'fixed': [], 'designable': []}
        res_designable_dict = {}
        for res_1pos in range(1, seq_len + 1):

            # Determine category of residue
            if res_1pos in fixed_designable_variable_indices['fixed']:
                category = 'fixed'
            elif res_1pos in fixed_designable_variable_indices['designable']:
                category = 'designable'
            else:
                category = 'variable'
            
            # Assign residue to final dict based on category and associated probability
            prob = probs_dict[category]
            random_prob = np.random.rand() # Random number between 0 and 1. Outputting probs of fixed/designable based on this value
            if random_prob >= prob:
                final_fixed_designable_dict['fixed'].append(res_1pos)
                res_designable_dict[res_1pos] = 'fixed'
            else: # Assign to designable random prob 
                final_fixed_designable_dict['designable'].append(res_1pos)
                res_designable_dict[res_1pos] = 'designable'
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
            if designability_attribute == 'fixed':
                # Preserve residue lines as is, but adjust chain to 'A' and residue numbering (if chosen to be fixed and thus, included in template PDB)
                for line in lines:
                    new_line = line[:pdb_chain_index] + 'A' + line[pdb_chain_index+1:pdb_residue_index.start] + str(residue_number).rjust(4) + line[pdb_residue_index.stop:]
                    pdb_lines.append(new_line)
                
                # Extract amino acid for seq input
                line = lines[0] # Take first line of residue to extract amino acid
                aa3 = line[17:20].strip() # Extract 3-letter amino acid code
                aa1 = seq.ProteinSequence().convert_letter_3to1(aa3)
                seq_input += aa1
            else: # designable, skip residue lines and add 'X' to seq input as model will redesign them
                # Skip adding lines to pdb_lines
                seq_input += "X"
                continue # Skip designable residues as model will redesign them

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
    
    def create_protein_hunter_inputs(self, output_file_path, linker_length=20, cdr_extend: int = 3, probs_dict: dict = {'fixed': 0.0, 'variable': 0.7, 'designable': 1.01}):
        """ Creates the single chain PDB file and extract seq input for Protein Hunter inputs """
        
        # 1. Extract chain sequences
        seq_dict, paired_seq = self.extract_chain_sequences()
        # 2. Annotate paired sequence
        annotated_paired_seq = self.annotate_paired_sequence(paired_seq)
        # 3. Create fixed/designable indices
        fixed_designable_variable_indices, paratope_indices = self.create_fixed_designable_variable_indices(annotated_paired_seq, linker_length=linker_length, 
                                                                                                            cdr_extend=cdr_extend)
        # 4. Create Final Fixed and Designable Dict 
        final_fixed_designable_dict, sc_res_designable_dict = self.create_final_fixed_designable_dict(fixed_designable_variable_indices, probs_dict=probs_dict)
        # Create single chain PDB file with linker
        seq_input = self.create_sc_pdb_file(output_file_path=output_file_path, res_designable_dict= sc_res_designable_dict, linker_length=linker_length)
        
        # Convert to mmCIF format
        output_cif_path = self.convert_pdb_to_cif(input_pdb_path=output_file_path)

        print("Number of Fixed Residues:", len(final_fixed_designable_dict['fixed']))
        print("Number of Designable Residues:", len(final_fixed_designable_dict['designable']))
        print("✅ Created Protein Hunter Inputs")
        return seq_dict, paired_seq, seq_input, output_cif_path, final_fixed_designable_dict, paratope_indices