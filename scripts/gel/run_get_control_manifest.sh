#!/usr/bin/env bash
# run_get_control_manifest.sh
# Run get_control_manifest.py on the GEL HPC.
#
# Setup (run once before using this script):
#   ssh rdawes@doublehelix.helix.prod.aws.gel.ac
#   cd /re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature
#
# Then run: bash scripts/gel/run_get_control_manifest.sh

/resources/conda/miniconda3/envs/py3pypirev2/bin/python -u scripts/gel/get_control_manifest.py
