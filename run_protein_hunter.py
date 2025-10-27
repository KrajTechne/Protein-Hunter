import subprocess

cmd = [
    "python", "boltz/protein_hunter.py",
    "--num_designs", "3",
    "--num_cycles", "7",
    "--protein_seqs", "AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNAPYAAALE",
    "--protein_ids", "B",
    "--protein_msas", "",
    "--gpu_id", "2",
    "--name", "PDL1_mix_aa",
    "--min_design_protein_length", "90",
    "--max_design_protein_length", "150",
    "--high_iptm_threshold", "0.7",
    "--use_msa_for_af3",
    "--plot"
]
subprocess.run(cmd, check=True)
