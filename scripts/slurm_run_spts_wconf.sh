#!/bin/bash
#
# Runs an SPTS analysis sweep over every *.conf in data/confs/ against a
# single .cxi file, one SLURM array task per config (run_spts_wconf.py
# itself is single-process, so no MPI needed here -- see run_spts.py -m
# if you want an individual analysis to run under MPI instead).
#
# Update --array below to match the number of *.conf files in CONF_DIR
# (currently 8: data/confs/spts{25,5}t{20,15,10,5}.conf).
#
# Usage: sbatch slurm_run_spts_wconf.sh <cxi_file>
#   e.g. sbatch slurm_run_spts_wconf.sh data01271.cxi
#
# Set the name of the job
#SBATCH --job-name=spts_sweep
#
# One task per array index, one CPU per task
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#
# The partition to use.
#SBATCH -p regular
#
# One array task per config file (indices 0-7 for 8 configs)
#SBATCH --array=0-7
#
# Log files, one pair per array task
#SBATCH --output=.slurm-output/spts_sweep_%A_%a.out
#SBATCH --error=.slurm-output/spts_sweep_%A_%a.err

CONF_DIR=/home/pat/spts-ana/data/confs
CXI_FILE=$1

if [ -z "$CXI_FILE" ]; then
    echo "Usage: sbatch slurm_run_spts_wconf.sh <cxi_file>" >&2
    exit 1
fi

CONFS=("$CONF_DIR"/*.conf)
CONF="${CONFS[$SLURM_ARRAY_TASK_ID]}"

echo "Task $SLURM_ARRAY_TASK_ID: running SPTS with config $CONF on $CXI_FILE"
run_spts_wconf.py "$CONF" "$CXI_FILE"
