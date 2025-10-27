import os
import sys
import argparse
from pathlib import Path
from protein_hunter_utils import *
import pandas as pd
import yaml
import torch
import random
import numpy as np
import copy

def parse_args():
    parser = argparse.ArgumentParser(
        description="Boltz protein design with cycle optimization"
    )
    parser.add_argument("--gpu_id", default=0, type=int)
    parser.add_argument("--grad_enabled", action="store_true", default=False)
    parser.add_argument("--name", default="target_name_is_missing", type=str)
    parser.add_argument(
        "--mode", default="binder", choices=["binder", "unconditional"], type=str
    )
    parser.add_argument("--num_designs", default=50, type=int)
    parser.add_argument("--num_cycles", default=5, type=int)
    parser.add_argument("--binder_chain", default="A", type=str)
    parser.add_argument("--min_design_protein_length", default=100, type=int)
    parser.add_argument("--max_design_protein_length", default=150, type=int)
    parser.add_argument("--protein_ids", default="B", type=str)
    parser.add_argument(
        "--protein_seqs",
        default="",
        type=str,
    )
    parser.add_argument("--protein_msas", default="empty", type=str)
    parser.add_argument("--cyclics", default="", type=str)
    parser.add_argument("--ligand_id", default="B", type=str)
    parser.add_argument("--ligand_smiles", default="", type=str)
    parser.add_argument("--ligand_ccd", default="", type=str)
    parser.add_argument(
        "--nucleic_type", default="dna", choices=["dna", "rna"], type=str
    )
    parser.add_argument("--nucleic_id", default="B", type=str)
    parser.add_argument("--nucleic_seq", default="", type=str)
    parser.add_argument(
        "--template_path", default="", type=str
    )  # can be "2VSM", or path(s) to .cif/.pdb, multiple allowed separated by comma
    parser.add_argument("--mediator_chain", default="", type=str)
    parser.add_argument("--template_chain_id", default="", type=str)
    parser.add_argument("--no_potentials", action="store_true")
    parser.add_argument("--diffuse_steps", default=200, type=int)
    parser.add_argument("--recycling_steps", default=3, type=int)
    parser.add_argument("--boltz_model_version", default="boltz2", type=str)
    parser.add_argument(
        "--boltz_model_path",
        default="~/.boltz/boltz2_conf.ckpt",
        type=str,
    )
    parser.add_argument(
        "--ccd_path", default="~/.boltz/mols", type=str
    )
    parser.add_argument("--randomly_kill_helix_feature", action="store_true")
    parser.add_argument("--negative_helix_constant", default=0.2, type=float)
    parser.add_argument("--logmd", action="store_true")
    parser.add_argument("--save_dir", default="", type=str)
    parser.add_argument("--add_constraints", action="store_true")
    parser.add_argument("--contact_residues", default="", type=str)
    parser.add_argument("--omit_AA", default="C", type=str)
    parser.add_argument("--all_X", action="store_true", default=False)
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Plot cycles figs per run (requires matplotlib)",
    )

    # NEW: Add constraint_target_chain argument
    parser.add_argument(
        "--constraint_target_chain",
        default="B",
        type=str,
        help="Target chain for constraints and contact calculations",
    )

    # NEW: Add no_contact_filter argument
    parser.add_argument(
        "--no_contact_filter",
        action="store_true",
        help="Do not filter or restart for unbound contact residues at cycle 0",
    )
    parser.add_argument("--max_contact_filter_retries", default=6, type=int)
    parser.add_argument("--contact_cutoff", default=15.0, type=float)

    parser.add_argument("--alphafold_dir", default="~/alphafold3", type=str)
    parser.add_argument("--af3_docker_name", default="alphafold3_yc", type=str)
    parser.add_argument(
        "--af3_database_settings", default="~/alphafold3/alphafold3_data_save", type=str
    )
    parser.add_argument(
        "--hmmer_path",
        default="~/.conda/envs/alphafold3_venv",
        type=str,
    )
    parser.add_argument("--use_msa_for_af3", action="store_true")
    parser.add_argument(
        "--work_dir", default="", type=str
    )

    # temp and bias params
    parser.add_argument("--temperature_start", default=0.05, type=float)
    parser.add_argument("--temperature_end", default=0.001, type=float)
    parser.add_argument("--alanine_bias_start", default=-0.5, type=float)
    parser.add_argument("--alanine_bias_end", default=-0.2, type=float)
    parser.add_argument("--alanine_bias", action="store_true")

    parser.add_argument("--high_iptm_threshold", default=0.8, type=float)

    return parser.parse_args()


def main():
    args = parse_args()
    no_potentials = args.no_potentials
    grad_enabled = args.grad_enabled
    diffuse_steps = args.diffuse_steps
    recycling_steps = args.recycling_steps
    boltz_model_path = args.boltz_model_path
    ccd_path = Path(args.ccd_path)
    ccd_lib = load_canonicals(os.path.expanduser(str(ccd_path)))
    logmd = args.logmd

    predict_args = {
        "recycling_steps": recycling_steps,
        "sampling_steps": diffuse_steps,
        "diffusion_samples": 1,
        "write_confidence_summary": True,
        "write_full_pae": False,
        "write_full_pde": False,
        "max_parallel_samples": 1,
    }

    device = f"cuda:{args.gpu_id}" if torch.cuda.is_available() else "cpu"
    print("Using device:", device)
    boltz_model = get_boltz_model(
        checkpoint=boltz_model_path,
        predict_args=predict_args,
        device=device,
        model_version=args.boltz_model_version,
        no_potentials=no_potentials,
        grad_enabled=grad_enabled,
    )
    designer = LigandMPNNWrapper(
        "./LigandMPNN/run.py"
    )


    # ========== Configuration ==========
    mode = args.mode
    name = args.name
    save_dir = (
        args.save_dir
        if args.save_dir
        else f"./results/{name}"
    )
    protein_hunter_save_dir = f"{save_dir}/0_protein_hunter_design"
    os.makedirs(save_dir, exist_ok=True)
    os.makedirs(protein_hunter_save_dir, exist_ok=True)
    num_designs = args.num_designs
    num_cycles = args.num_cycles
    binder_chain = args.binder_chain
    min_design_protein_length = args.min_design_protein_length
    max_design_protein_length = args.max_design_protein_length

    protein_ids = args.protein_ids
    protein_seqs = args.protein_seqs
    protein_msas = args.protein_msas
    cyclics = args.cyclics
    ligand_id = args.ligand_id
    ligand_smiles = args.ligand_smiles
    ligand_ccd = args.ligand_ccd
    nucleic_type = args.nucleic_type
    nucleic_id = args.nucleic_id
    nucleic_seq = args.nucleic_seq

    def smart_split(s):
        if "," in s:
            return [x.strip() for x in s.split(",")]
        elif ":" in s:
            return [x.strip() for x in s.split(":")]
        else:
            return [x.strip() for x in s.split()] if s else []

    if args.template_path:
        template_path_list = smart_split(args.template_path)
        template_chain_id_list = smart_split(args.template_chain_id) if args.template_chain_id else []
        template_files = [get_cif(tp) for tp in template_path_list]
        print("Parsed template(s):", template_files)
    else:
        template_files = []
        template_chain_id_list = []
    constraint_target_chain = (
        args.constraint_target_chain
    )  # NEW: target chain for constraints and contact calculations

    if ":" not in protein_ids and "," not in protein_ids:
        protein_ids_list = [protein_ids.strip()] if protein_ids.strip() else []
        protein_seqs_list = [protein_seqs.strip()] if protein_seqs.strip() else []
        protein_msas_list = [protein_msas.strip()] if protein_msas.strip() else []
        cyclics_list = [cyclics.strip()] if cyclics.strip() else []
    else:
        protein_ids_list = smart_split(protein_ids)
        protein_seqs_list = smart_split(protein_seqs)
        protein_msas_list = (
            smart_split(protein_msas) if protein_msas else [""] * len(protein_ids_list)
        )
        cyclics_list = (
            smart_split(cyclics) if cyclics else ["False"] * len(protein_ids_list)
        )
    max_len = max(
        len(protein_ids_list),
        len(protein_seqs_list),
        len(protein_msas_list),
        len(cyclics_list),
    )
    while len(protein_ids_list) < max_len:
        protein_ids_list.append("")
    while len(protein_seqs_list) < max_len:
        protein_seqs_list.append("")
    while len(protein_msas_list) < max_len:
        protein_msas_list.append("")
    while len(cyclics_list) < max_len:
        cyclics_list.append("False")

    if mode == "unconditional":
        data = {
            "sequences": [{"protein": {"id": ["A"], "sequence": "X", "msa": "empty"}}]
        }
    else:
        sequences = []

        # Step 1: For each entry, group by sequence.
        from collections import defaultdict

        seq_to_indices = defaultdict(list)
        for idx, seq in enumerate(protein_seqs_list):
            if seq:
                seq_to_indices[seq].append(idx)

        # Step 2: Assign a representative index for each unique sequence,
        msa_final_for_seq = {}
        for seq, idx_list in seq_to_indices.items():
            chosen_msa = None
            # Find first non-empty MSA ('' means to generate; 'empty' is explicit single sequence; else path)
            for i in idx_list:
                val = protein_msas_list[i]
                if val and val != "":
                    chosen_msa = val
                    break
            # If no non-empty found, keep empty string ("") to indicate generation, else assign whichever
            if chosen_msa is None:
                chosen_msa = ""
            msa_final_for_seq[seq] = chosen_msa

        seq_to_final_msa = {}
        for seq, chosen_msa in msa_final_for_seq.items():
            idx0 = seq_to_indices[seq][0]
            pid0 = protein_ids_list[idx0] if protein_ids_list[idx0] else f"CHAIN_{idx0}"
            if chosen_msa == "":
                print("Processing MSA for", pid0)
                print("pseq", seq)
                print("pmsa", chosen_msa)
                msa_value = process_msa(pid0, seq, Path(protein_hunter_save_dir))
                seq_to_final_msa[seq] = str(msa_value)
            elif chosen_msa == "empty":
                seq_to_final_msa[seq] = "empty"
            else:
                seq_to_final_msa[seq] = chosen_msa

        # (2) Build the sequences list, using the canonical MSA for each sequence (all identical sequences share same MSA)
        for idx, (pid, seq, cyc) in enumerate(
                zip(protein_ids_list, protein_seqs_list, cyclics_list)
        ):
            if not pid or not seq:
                continue
            final_msa = seq_to_final_msa.get(seq, "empty")
            cyc_val = cyc.lower() in ["true", "1", "yes"]
            sequences.append(
                {
                    "protein": {
                        "id": [pid],
                        "sequence": seq,
                        "msa": final_msa,
                        "cyclic": cyc_val,
                    }
                }
            )

        # Add binder chain after all targets.
        sequences.append(
            {
                "protein": {
                    "id": [binder_chain],
                    "sequence": "X",
                    "msa": "empty",
                    "cyclic": False,
                }
            }
        )

        # Support both small-molecule ligand and nucleic acid together, always use LigandMPNN
        if ligand_smiles:
            sequences.append({"ligand": {"id": [ligand_id], "smiles": ligand_smiles}})
        elif ligand_ccd:
            sequences.append({"ligand": {"id": [ligand_id], "ccd": ligand_ccd}})
        if nucleic_seq:
            sequences.append(
                {nucleic_type: {"id": [nucleic_id], "sequence": nucleic_seq}}
            )

        # ===================== MULTIPLE TEMPLATES SUPPORT =====================
        templates = []
        # If multiple templates provided, add all as template blocks
        for i, template_file in enumerate(template_files):
            t_block = {}
            if template_file.endswith(".cif"):
                t_block["cif"] = template_file
            else:
                t_block["pdb"] = template_file
            # assign chain id from template_chain_id(s) if available and matches index
            if template_chain_id_list and i < len(template_chain_id_list):
                t_block["chain_id"] = template_chain_id_list[i]
            elif args.template_chain_id and not template_chain_id_list:
                t_block["chain_id"] = args.template_chain_id
            templates.append(t_block)
        data = {"sequences": sequences}
        if templates:
            data["templates"] = templates
        # ======================================================================

    data["sequences"] = sorted(
        data["sequences"], key=lambda entry: list(entry.values())[0]["id"][0]
    )

    if args.add_constraints:
        residues = args.contact_residues.split(",")
        contacts = []
        # Only use the constraint_target_chain as target for constraints regardless of protein_ids_list
        contacts.extend(
            [
                [constraint_target_chain, int(res)]
                for res in residues
                if res.strip() != ""
            ]
        )
        constraints = {
            "pocket": {
                "binder": binder_chain,
                "contacts": contacts,
            }
        }
        data["constraints"] = [constraints]
        pocket_conditioning = True
    else:
        pocket_conditioning = False

    print("✅ Configuration ready.")
    print("Mode:", mode)
    print("Name:", name)
    print("Save directory:", protein_hunter_save_dir)
    print("Data dictionary:\n", data)

    def sample_seq(length, all_X=False):
        if all_X:
            amino_acids = "X" * 200
        else:
            amino_acids = "ADEFGHIKLMNQRSTVWY" + "X" * 18
        return "".join(random.choices(amino_acids, k=length))

    def shallow_copy_tensor_dict(d):
        if isinstance(d, dict):
            return {k: shallow_copy_tensor_dict(v) for k, v in d.items()}
        elif isinstance(d, list):
            return [shallow_copy_tensor_dict(x) for x in d]
        elif isinstance(d, torch.Tensor):
            return d.detach().clone()
        else:
            return d

    # ========== Helper: check contacts after structure prediction ==========
    clean_memory()
    any_ligand_or_nucleic = ligand_smiles or ligand_ccd or nucleic_seq
    if any_ligand_or_nucleic:
        model_type = "ligand_mpnn"
    else:
        model_type = "soluble_mpnn"
    print("model_type", model_type)

    alanine_bias = args.alanine_bias
    temperature_start = args.temperature_start
    temperature_end = args.temperature_end
    alanine_bias_start = args.alanine_bias_start
    alanine_bias_end = args.alanine_bias_end

    all_run_metrics = []

    for design_id in range(num_designs):
        run_id = str(design_id)
        run_save_dir = os.path.join(protein_hunter_save_dir, f"run_{run_id}")
        os.makedirs(run_save_dir, exist_ok=True)

        data_cp = copy.deepcopy(data)
        print(f"\n=== Starting Design Run {run_id} ===")
        best_iptm = float("-inf")
        best_structure = None
        best_output = None
        best_pdb_filename = None
        best_cycle_idx = -1
        best_alanine_percentage = None
        run_metrics = {"run_id": run_id}
        binder_length = random.randint(
            min_design_protein_length, max_design_protein_length
        )
        if mode == "unconditional":
            data_cp["sequences"][0]["protein"]["sequence"] = sample_seq(
                binder_length, all_X=args.all_X
            )
        else:
            for seq_entry in data_cp["sequences"]:
                if (
                    "protein" in seq_entry
                    and binder_chain in seq_entry["protein"]["id"]
                ):
                    seq_entry["protein"]["sequence"] = sample_seq(
                        binder_length, all_X=args.all_X
                    )
                    print(seq_entry["protein"]["sequence"])
                    break
        print(f"Binder initial sequence length: {binder_length}")

        # ========== Cycle 0 structure prediction, with contact filtering check ==========
        contact_filter_attempt = 0
        max_contact_filter_retries = args.max_contact_filter_retries
        cycle0_contacts_verified = False
        initial_contact_residues = args.contact_residues.strip()
        no_contact_filter = getattr(args, "no_contact_filter", False)
        # Use constraint_target_chain instead of guessing target chain for contacts
        target_chain_for_contacts = constraint_target_chain
        while True:
            output, structure = run_prediction(
                data_cp,
                binder_chain,
                randomly_kill_helix_feature=args.randomly_kill_helix_feature,
                negative_helix_constant=args.negative_helix_constant,
                boltz_model=boltz_model,
                ccd_lib=ccd_lib,
                ccd_path=ccd_path,
                logmd=logmd,
                device=device,
                boltz_model_version=args.boltz_model_version,
                pocket_conditioning=pocket_conditioning,
            )
            structure.atoms["coords"] = (
                output["coords"][0]
                .detach()
                .cpu()
                .numpy()[: structure.atoms["coords"].shape[0], :]
            )
            pdb_filename = f"{run_save_dir}/{name}_run_{run_id}_predicted_cycle_0.pdb"
            plddts = output["plddt"].detach().cpu().numpy()[0]
            save_pdb(structure, output["coords"], plddts, pdb_filename)

            # If contact_residues specified, check if binder contacts them after cycle 0 structure prediction
            contact_check_okay = True
            if initial_contact_residues and not no_contact_filter:
                print(
                    f"Checking if binder binds target contact residues ({initial_contact_residues}) on {target_chain_for_contacts} after run_prediction cycle 0, attempt {contact_filter_attempt + 1}"
                )
                try:
                    binds = binder_binds_contacts(
                        pdb_filename,
                        binder_chain,
                        target_chain_for_contacts,
                        initial_contact_residues,
                        cutoff=args.contact_cutoff,
                    )
                    if not binds:
                        print(
                            f"❌ Binder does NOT contact required residues {initial_contact_residues} after structure prediction. Retrying run_prediction..."
                        )
                        contact_check_okay = False
                except Exception as e:
                    print(f"WARNING: Could not perform binder-contact check: {e}")
                    contact_check_okay = True  # fail open if prody/numpy not installed
            else:
                binds = True
            if contact_check_okay:
                break
            else:
                # Redo structure prediction with a new random X sequence if contact is NOT satisfied
                contact_filter_attempt += 1
                if contact_filter_attempt >= max_contact_filter_retries:
                    print(
                        "WARNING: Maximum retries for contact filtering reached. Proceeding anyway."
                    )
                    break
                print("Resampling initial sequence and predicting again.")
                # Resample initial sequence to new random sequence of same binder_length
                if mode == "unconditional":
                    data_cp["sequences"][0]["protein"]["sequence"] = sample_seq(
                        binder_length, all_X=args.all_X
                    )
                else:
                    for seq_entry in data_cp["sequences"]:
                        if (
                            "protein" in seq_entry
                            and binder_chain in seq_entry["protein"]["id"]
                        ):
                            seq_entry["protein"]["sequence"] = sample_seq(
                                binder_length, all_X=args.all_X 
                            )
                            print(
                                f"New resampled binder sequence: {seq_entry['protein']['sequence']}"
                            )
                            break
                clean_memory()

        (
            iptm_values,
            plddt_values,
            iplddt_values,
            cycle_numbers,
            alanine_counts,
            pair_iptm_1_avg,
            pair_iptm_2_avg,
        ) = (
            [],
            [],
            [],
            [],
            [],
            [],
            [],
        )
        logits_list, seq_str_list = [], []
        binder_chain_idx = chain_to_number[binder_chain]
        pair_chains = output['pair_chains_iptm']

        # Calculate mean pairwise ipTM for binder_chain vs all other chains (see @file_context_0)
        if len(pair_chains) > 1:
            # Calculate the average of the pairwise iptm numbers (not their sum)
            values = [
                (pair_chains[binder_chain_idx][i].detach().cpu().numpy() + pair_chains[i][binder_chain_idx].detach().cpu().numpy()) / 2.0
                for i in range(len(pair_chains)) if i != binder_chain_idx
            ]
            cycle_0_iptm = float(np.mean(values) if values else 0.0)
            print("cycle_0_iptm", cycle_0_iptm)
        else:
            cycle_0_iptm = 0.0
        run_metrics["cycle_0_iptm"] = float(cycle_0_iptm)

        run_metrics["cycle_0_plddt"] = float(
            output.get("complex_plddt", torch.tensor([0.0])).detach().cpu().numpy()[0]
        )
        run_metrics["cycle_0_iplddt"] = float(
            output.get("complex_iplddt", torch.tensor([0.0])).detach().cpu().numpy()[0]
        )
        run_metrics["cycle_0_alanine"] = 0  # no alanine before design

        for cycle in range(num_cycles):
            print(f"\n--- Run {run_id}, Cycle {cycle + 1} ---")
            alpha = (
                alanine_bias_start
                - (cycle / (num_cycles - 1)) * (alanine_bias_start - alanine_bias_end)
                if num_cycles > 1
                else alanine_bias_start
            )
            temperature = (
                temperature_start
                - (cycle / (num_cycles - 1)) * (temperature_start - temperature_end)
                if num_cycles > 1
                else temperature_start
            )

            print(f"Cycle {cycle + 1}: binder length = {binder_length}")

            design_kwargs = {
                "pdb_file": pdb_filename,
                "temperature": temperature,
            }
 
            design_kwargs["chains_to_design"] = binder_chain

            if cycle == 0:
                design_kwargs["omit_AA"] = f"{args.omit_AA},P"
            else:
                design_kwargs["omit_AA"] = args.omit_AA

            if alanine_bias:
                design_kwargs["bias_AA"] = f"A:{alpha}"
            seq_str, logits = design_sequence(designer, model_type, **design_kwargs)
            seq = seq_str.split(":")[chain_to_number[binder_chain]]
            assert len(seq) == binder_length, (
                f"Sequence length mismatch: {len(seq)} != {binder_length}"
            )
            # Calculate alanine percentage for the designed sequence
            alanine_count = seq.count("A")
            alanine_percentage = (
                alanine_count / binder_length if binder_length != 0 else 0.0
            )
            for seq_entry in data_cp["sequences"]:
                if (
                    "protein" in seq_entry
                    and binder_chain in seq_entry["protein"]["id"]
                ):
                    seq_entry["protein"]["sequence"] = seq
                    break

            alanine_counts.append(alanine_count)
            output, structure = run_prediction(
                data_cp,
                binder_chain,
                seq=seq,
                randomly_kill_helix_feature=args.randomly_kill_helix_feature
                if cycle == 0
                else False,
                negative_helix_constant=args.negative_helix_constant
                if cycle == 0
                else 0.0,
                boltz_model=boltz_model,
                ccd_lib=ccd_lib,
                ccd_path=ccd_path,
                logmd=False,
                device=device,
            )
            # current_iptm = (
            #     output.get("iptm", torch.tensor([0.0])).detach().cpu().numpy()[0]
            # )
            current_chain_idx = chain_to_number[binder_chain]
            pair_chains = output['pair_chains_iptm']
            if len(pair_chains) > 1:
                values = [
                    (pair_chains[current_chain_idx][i].detach().cpu().numpy() + pair_chains[i][current_chain_idx].detach().cpu().numpy()) / 2.0
                    for i in range(len(pair_chains)) if i != current_chain_idx
                ]
                current_iptm = float(np.mean(values) if values else 0.0)
            else:
                current_iptm = 0.0
            # Only consider best structure if alanine percentage <= 20%
            if alanine_percentage <= 0.20 and current_iptm > best_iptm:
                best_iptm = current_iptm
                best_structure = copy.deepcopy(structure)
                best_output = shallow_copy_tensor_dict(output)
                best_pdb_filename = (
                    f"{run_save_dir}/{name}_run_{run_id}_best_structure.pdb"
                )
                best_plddts = best_output["plddt"].detach().cpu().numpy()[0]
                save_pdb(
                    best_structure,
                    best_output["coords"],
                    best_plddts,
                    best_pdb_filename,
                )
                best_cycle_idx = cycle + 1
                best_alanine_percentage = alanine_percentage


            iptm_values.append(current_iptm)
            curr_plddt = float(
                output.get("complex_plddt", torch.tensor([0.0]))
                .detach()
                .cpu()
                .numpy()[0]
            )
            curr_iplddt = float(
                output.get("complex_iplddt", torch.tensor([0.0]))
                .detach()
                .cpu()
                .numpy()[0]
            )
            plddt_values.append(curr_plddt)
            iplddt_values.append(curr_iplddt)
            cycle_numbers.append(cycle + 1)

            run_metrics[f"cycle_{cycle + 1}_iptm"] = float(current_iptm)
            run_metrics[f"cycle_{cycle + 1}_plddt"] = float(curr_plddt)
            run_metrics[f"cycle_{cycle + 1}_iplddt"] = float(curr_iplddt)
            run_metrics[f"cycle_{cycle + 1}_alanine"] = int(alanine_count)

            print(
                f"Run {run_id}, Cycle {cycle + 1} sequence (len={binder_length}): {seq}"
            )
            print(
                f"ipTM: {iptm_values[-1]:.2f} pLDDT: {plddt_values[-1]:.2f} iPLDDT: {iplddt_values[-1]:.2f} Alanine count: {alanine_counts[-1]}"
            )
            pdb_filename = f"{run_save_dir}/{name}_run_{run_id}_predicted_cycle_{cycle + 1}.pdb"
            plddts = output["plddt"].detach().cpu().numpy()[0]
            save_pdb(structure, output["coords"], plddts, pdb_filename)
            clean_memory()

            last_cycle_iptm = run_metrics.get(f"cycle_{cycle + 1}_iptm", 0.0)
            # Only allow high_iptm designs with <= 20% alanine
            if alanine_percentage <= 0.20:
                high_iptm = last_cycle_iptm > args.high_iptm_threshold
            else:
                high_iptm = False

            # ===================== BEGIN PATCH FOR CONTACT CHECK WHEN SAVING HIGH IPTM YAML =====================
            save_yaml_this_design = high_iptm
            contact_residues_to_check = args.contact_residues.strip()
            yaml_contact_binding_passed = True
            if high_iptm and contact_residues_to_check:
                target_chain_for_contacts = constraint_target_chain
                this_cycle_pdb_filename = f"{run_save_dir}/{name}_run_{run_id}_predicted_cycle_{cycle + 1}.pdb"
                try:
                    contact_binds = binder_binds_contacts(
                        this_cycle_pdb_filename,
                        binder_chain,
                        target_chain_for_contacts,
                        contact_residues_to_check,
                        cutoff=args.contact_cutoff,
                    )
                    if not contact_binds:
                        print(
                            f"⛔️ Not saving YAML: binder does not contact required residues ({contact_residues_to_check}) (run={run_id}, cycle={cycle + 1})"
                        )
                        save_yaml_this_design = False
                        yaml_contact_binding_passed = False
                    else:
                        print(
                            f"✅ binder_binds_contacts PASSED when saving YAML for run={run_id} cycle={cycle + 1}"
                        )
                except Exception as e:
                    print(
                        f"WARNING: Exception occurred during binder_binds_contacts for high iptm save: {e}"
                    )
                    # Failsafe: proceed with saving
                    save_yaml_this_design = True

            if save_yaml_this_design:
                high_iptm_yaml_dir = os.path.join(save_dir, "high_iptm_yaml")
                os.makedirs(high_iptm_yaml_dir, exist_ok=True)
                yaml_filename = os.path.join(
                    high_iptm_yaml_dir,
                    f"{name}_run_{run_id}_cycle_{cycle + 1}_output.yaml",
                )
                with open(yaml_filename, "w") as f:
                    yaml.dump(data_cp, f, default_flow_style=False)
                if not yaml_contact_binding_passed and contact_residues_to_check:
                    print(
                        "⚠️ Saved YAML despite contact_binds_contacts warning (see above)."
                    )
                print(
                    f"✅ Saved run {run_id} cycle {cycle + 1} YAML to {yaml_filename} (iptm={last_cycle_iptm:.3f}, alanine%={alanine_percentage * 100:.1f})"
                )
            else:
                if last_cycle_iptm > args.high_iptm_threshold:
                    print(
                        f"⛔️ Skipped saving to high iptm: Run {run_id} Cycle {cycle + 1} (alanine%={alanine_percentage * 100:.1f})"
                    )

        # Visualization/metrics only for best designs with <= 20% alanine
        if (
            best_structure is not None
            and best_output is not None
            and best_pdb_filename is not None
            and best_alanine_percentage is not None
            and best_alanine_percentage <= 0.20
        ):
            print(
                f"\n🔝 Visualizing best structure for run {run_id} (Highest ipTM: {best_iptm:.3f}, alanine%={best_alanine_percentage * 100:.1f})\n→ {best_pdb_filename}"
            )
            plot_from_pdb(best_pdb_filename)
        else:
            print(
                f"\nNo structure was generated for run {run_id} (no eligible best design with <= 20% alanine)."
            )

        if best_alanine_percentage is not None and best_alanine_percentage <= 0.20:
            run_metrics["best_iptm"] = float(best_iptm)
            run_metrics["best_cycle"] = best_cycle_idx
            run_metrics["best_plddt"] = float(
                best_output.get("complex_plddt", torch.tensor([0.0]))
                .detach()
                .cpu()
                .numpy()[0]
            )
        else:
            run_metrics["best_iptm"] = float("nan")
            run_metrics["best_cycle"] = None
            run_metrics["best_plddt"] = float("nan")
        all_run_metrics.append(run_metrics)

        # Plot per-run figure to run subfolder
        if args.plot:
            import matplotlib.pyplot as plt

            fig, axs = plt.subplots(1, 4, figsize=(12, 3))
            colors = ["#9B59B6", "#E94560", "#FF7F11", "#2ECC71"]
            metrics = [
                (
                    "iPTM",
                    [
                        run_metrics.get(f"cycle_{i}_iptm", float("nan"))
                        for i in range(num_cycles + 1)
                    ],
                    0,
                    1,
                    "{:.3f}",
                ),
                (
                    "pLDDT",
                    [
                        run_metrics.get(f"cycle_{i}_plddt", float("nan"))
                        for i in range(num_cycles + 1)
                    ],
                    0,
                    1,
                    "{:.1f}",
                ),
                (
                    "iPLDDT",
                    [
                        run_metrics.get(f"cycle_{i}_iplddt", float("nan"))
                        for i in range(num_cycles + 1)
                    ],
                    0,
                    1,
                    "{:.1f}",
                ),
                (
                    "Alanine Count",
                    [
                        run_metrics.get(f"cycle_{i}_alanine", float("nan"))
                        for i in range(num_cycles + 1)
                    ],
                    0,
                    max(
                        [
                            run_metrics.get(f"cycle_{i}_alanine", 0)
                            for i in range(num_cycles + 1)
                        ]
                    )
                    + 2
                    if alanine_counts
                    else 10,
                    "{}",
                ),
            ]
            design_cycles = list(range(num_cycles + 1))
            for ax, (label, values, ymin, ymax, fmt), color in zip(
                axs, metrics, colors
            ):
                ax.plot(
                    design_cycles,
                    values,
                    "o-",
                    color=color,
                    linewidth=2,
                    markersize=6,
                    markerfacecolor="white",
                    markeredgewidth=2,
                )
                ax.set(
                    xlabel="Design Iteration",
                    ylabel=label,
                    title=f"{label} (Run {run_id})",
                    xticks=design_cycles,
                    ylim=(ymin, ymax),
                )
                ax.grid(True, alpha=0.3, linestyle="--")
                ax.spines[["top", "right"]].set_visible(False)
                for x, y in zip(design_cycles, values):
                    ax.annotate(
                        fmt.format(y) if not pd.isnull(y) else "",
                        (x, y),
                        textcoords="offset points",
                        xytext=(0, 8),
                        ha="center",
                        fontsize=9,
                    )
            plt.tight_layout()
            plt.savefig(
                f"{run_save_dir}/{name}_run_{run_id}_design_cycle_results.png", dpi=300
            )
            plt.close()

    # ===== Save all metrics to a single CSV =====
    columns = ["run_id"]
    for i in range(num_cycles + 1):
        columns.extend(
            [
                f"cycle_{i}_iptm",
                f"cycle_{i}_plddt",
                f"cycle_{i}_iplddt",
                f"cycle_{i}_alanine",
            ]
        )
    columns.extend(["best_iptm", "best_cycle"])

    summary_csv = os.path.join(save_dir, "summary_all_runs.csv")
    df = pd.DataFrame(all_run_metrics)
    for col in columns:
        if col not in df.columns:
            df[col] = float("nan")
    df = df[columns]
    df.to_csv(summary_csv, index=False)
    print(f"\n✅ All run/cycle metrics saved to {summary_csv}")

    # Run downstream validation: alphafold_step and rosetta_step
    alphafold_dir = os.path.expanduser(args.alphafold_dir)
    af3_docker_name = args.af3_docker_name
    af3_database_settings = os.path.expanduser(args.af3_database_settings)
    hmmer_path = os.path.expanduser(args.hmmer_path)
    work_dir = os.path.expanduser(args.work_dir) or os.getcwd()
    binder_id = binder_chain
    target_type = "protein" if protein_ids_list else "nucleic" if nucleic_type else "small_molecule" if ligand_smiles else "metal" if ligand_ccd else "protein"
    success_dir = f"{save_dir}/1_af3_rosetta_validation"
    af_pdb_dir = f"{success_dir}/03_af_pdb_success"
    af_pdb_dir_apo = f"{success_dir}/03_af_pdb_apo"
    high_iptm_yaml_dir = (
        os.path.join(save_dir, "high_iptm_yaml")
        if os.path.exists(os.path.join(save_dir, "high_iptm_yaml"))
        else save_dir
    )

    _, _, _, _ = run_alphafold_step(
        high_iptm_yaml_dir,
        alphafold_dir,
        af3_docker_name,
        af3_database_settings,
        hmmer_path,
        success_dir,
        work_dir,
        binder_id=binder_id,
        gpu_id=args.gpu_id,
        high_iptm=True,
        use_msa_for_af3=args.use_msa_for_af3,
    )

    run_rosetta_step(
        success_dir,
        af_pdb_dir,
        af_pdb_dir_apo,
        binder_id=binder_id,
        target_type=target_type,
    )

if __name__ == "__main__":
    main()