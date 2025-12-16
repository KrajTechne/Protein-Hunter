import os
from pydoc import doc
import numpy as np
import biotite
import gemmi
import biotite.sequence as seq
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
import py2Dmol

class ScfvTemplateConstructor:
    def __init__(self, fv_pdb_path: str):
        self.fv_pdb_path = fv_pdb_path
        self.atom_lines = self._read_pdb_file()
    
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
    
    def create_sc_pdb_file(self, output_file_path, linker_length=20):
        """ Writes a single chain of a PDB file to a new file, 
            with an optional linker of specified linker length added between first chain's last residue and second chain's first residue."""

        pdb_chain_index = 21
        pdb_residue_index = slice(23, 26)
        chain_A_lines = [lines for lines in self.atom_lines if lines.startswith('ATOM') and lines[pdb_chain_index] == 'A']
        chain_B_lines = [lines for lines in self.atom_lines if lines.startswith('ATOM') and lines[pdb_chain_index] == 'B']
        final_chain_A_residue_index = int(chain_A_lines[-1][pdb_residue_index].strip())

        new_chain_B_starting_residue_index = final_chain_A_residue_index + linker_length # Want to add in linker of specified size
        # Challenge: How to ensure groups of lines associated with same residue are modified together?
        # One way: Use a dictionary to map residue numbers to their corresponding lines
        residue_to_lines = {}
        for line in chain_B_lines:
            residue_number = int(line[pdb_residue_index].strip())
            residue_line_list = residue_to_lines.get(residue_number, [])
            residue_line_list.append(line)
            residue_to_lines[residue_number] = residue_line_list
    
        # Adjusting residue numbers and chain IDs for chain B -> want single chain and offset residue number by linker_length
        residue_lines_with_linker = []
        for residue_number, lines in residue_to_lines.items():
            new_residue_number = new_chain_B_starting_residue_index + residue_number
            for line in lines:
                new_line = line[:pdb_chain_index] + 'A' + line[pdb_chain_index+1:pdb_residue_index.start] + str(new_residue_number).rjust(3) + line[pdb_residue_index.stop:]
                residue_lines_with_linker.append(new_line)
    
        with open(output_file_path, 'w') as file:
            for line in chain_A_lines:
                file.write(line)
            for line in residue_lines_with_linker:
                file.write(line)
            file.write('END\n')
        print(f"Written to {output_file_path}")
    
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
    
    def create_protein_hunter_inputs(self, output_file_path, linker_length=20):
        """ Creates the single chain PDB file and extract seq input for Protein Hunter inputs """
        
        # Create single chain PDB file with linker
        self.create_sc_pdb_file(output_file_path=output_file_path, linker_length=linker_length)
        
        # Convert to mmCIF format
        output_cif_path = self.convert_pdb_to_cif(input_pdb_path=output_file_path)

        # Assuming Redesign area only restricted to linker region for now
        seq_dict, paired_seq = self.extract_chain_sequences()
        protein_hunter_mutated_seq = "X" * linker_length
        seq_input = protein_hunter_mutated_seq.join(chain_seq for chain_seq in seq_dict.values())
        return seq_dict, paired_seq, seq_input, output_cif_path