# Demo dataset

This directory contains a small simulated dataset for demonstrating input table formats and downstream script usage.

The demo is **not** intended to reproduce manuscript results. It is provided only to satisfy code-review and checklist requirements for a small dataset that can be used to inspect and test the downstream workflow.

Files:

- `demo_sample_metadata.tsv`: sample IDs and regions
- `demo_environmental_metadata.tsv`: small environmental metadata table
- `demo_mag_metadata.tsv`: MAG-level metadata
- `demo_mag_presence.tsv`: sample × MAG presence matrix
- `demo_arg_class_matrix.tsv`: MAG × ARG class count matrix
- `demo_vf_module_matrix.tsv`: MAG × virulence-associated module count matrix
- `demo_tree.nwk`: small example MAG phylogeny
- `run_demo.R`: lightweight demonstration script

Run:

```bash
Rscript demo/run_demo.R
```

Expected outputs:

- `results/demo_outputs/demo_sample_burdens.tsv`
- `results/demo_outputs/demo_model_summary.tsv`
- `results/demo_outputs/demo_readme.txt`

Expected runtime: less than one minute on a normal desktop computer.
