import argparse

from core import ProteinHunter_Chai


def parse_args():
    parser = argparse.ArgumentParser(description="ProteinHunter input options")

    # Job and sequence settings
    parser.add_argument("--jobname", type=str, default="test", help="Job name")
    parser.add_argument("--length", type=int, default=150, help="Protein length")
    parser.add_argument(
        "--percent_X", type=int, default=50, choices=[0, 50, 100], help="Percent X"
    )
    parser.add_argument("--seq", type=str, default="", help="Input sequence (optional)")
    parser.add_argument(
        "--target_seq",
        type=str,
        default="",
        help="Target sequence or SMILES (optional)",
    )
    parser.add_argument(
        "--cyclic",
        action="store_true",
        default=False,
        help="Enable cyclic topology",
    )

    # Optimization options
    parser.add_argument(
        "--n_trials",
        type=int,
        default=1,
        help="Number of parallel trials",
    )
    parser.add_argument(
        "--n_cycles",
        type=int,
        default=5,
        help="Number of optimization cycles",
    )
    parser.add_argument(
        "--n_recycles",
        type=int,
        default=3,
        help="Number of trunk recycles",
    )
    parser.add_argument(
        "--n_diff_steps",
        type=int,
        default=200,
        help="Diffusion steps",
    )
    parser.add_argument(
        "--hysteresis_mode",
        type=str,
        default="esm",
        choices=["templates", "esm", "partial_diffusion", "none"],
        help="Hysteresis/partial diffusion mode",
    )
    parser.add_argument(
        "--repredict",
        action="store_true",
        default=True,
        help="Re-predict final best without templates",
    )

    # MPNN options
    parser.add_argument("--omit_aa", type=str, default="", help="AA types to omit")
    parser.add_argument(
        "--bias_aa",
        type=str,
        default="",
        help="AA types to bias e.g. 'A:-2.0,P:-1.0,C:-0.5'",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=0.01,
        help="MPNN sampling temperature",
    )
    parser.add_argument(
        "--scale_temp_by_plddt",
        action="store_true",
        default=False,
        help="Scale temperature by pLDDT",
    )

    # Visual options
    parser.add_argument(
        "--show_visual",
        action="store_true",
        default=False,
        help="Show interactive visualization during run",
    )
    parser.add_argument(
        "--render_freq",
        type=int,
        default=100,
        help="Visualization refresh frequency",
    )
    parser.add_argument("--gpu_id", type=int, default=0, help="GPU ID")

    return parser.parse_args()


def main():
    args = parse_args()
    # Instantiate the main class and run the pipeline
    protein_hunter = ProteinHunter_Chai(args)
    protein_hunter.run_pipeline()


if __name__ == "__main__":
    main()
