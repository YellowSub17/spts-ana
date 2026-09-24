#!/bin/bash
#
# Runs an SPTS analysis sweep over every *.conf in data/confs/ against a
# single .cxi file, one SLURM array task per config. Each array task runs
# its analysis under MPI (run_spts_wconf.py -m), using the tasks requested
# below as MPI ranks -- rank 0 is a dedicated writer with no worker role,
# so --ntasks must be >=2 for any detection work to actually happen.
#
# Update --array below to match the number of *.conf files in CONF_DIR
# (currently 8: data/confs/spts{25,5}t{20,15,10,5}.conf). Note total cores
# requested = --ntasks x number of array tasks running concurrently, so
# tune --ntasks to your cluster's budget (8 ranks x 8 configs = 64 cores
# if the whole array runs at once).
#
# Usage: sbatch slurm_run_spts_wconf.sh <cxi_file>
#   e.g. sbatch slurm_run_spts_wconf.sh data01271.cxi
#
# Set the name of the job
#SBATCH --job-name=spts_sweep
#
# MPI ranks per array task (rank 0 = writer, the rest = workers)
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#
# The partition to use.
#SBATCH -p fast
#
# One array task per config file (indices 0-7 for 8 configs)
#SBATCH --array=0-7
#
# Log files, one pair per array task
#SBATCH --output=.slurm/spts_sweep_%A_%a.out
#SBATCH --error=.slurm/spts_sweep_%A_%a.err

CONF_DIR=/home/pat/spts-ana/data/confs
CXI_FILE=$1


source /home/pat/miniconda3/etc/profile.d/conda.sh
conda activate spts
export PATH="$CONDA_PREFIX/bin:$PATH"   # defensive, in case something on PATH still comes first
hash -r

echo $(which python)


echo $(date)


if [ -z "$CXI_FILE" ]; then
    echo "Usage: sbatch slurm_run_spts_wconf.sh <cxi_file>" >&2
    exit 1
fi

CONFS=("$CONF_DIR"/*.conf)
CONF="${CONFS[$SLURM_ARRAY_TASK_ID]}"

echo "Task $SLURM_ARRAY_TASK_ID: running SPTS with config $CONF on $CXI_FILE (MPI, $SLURM_NTASKS ranks)"
run_spts_wconf.py "$CONF" "$CXI_FILE" -m

echo $(date)
