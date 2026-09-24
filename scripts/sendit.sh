

for i in $(seq 1245 1260); do
    if [[ "$i" -eq 1251 || "$i" -eq 1260 || "$i" -eq 1262 || "$i" -eq 1270 ]]; then
        continue
    fi

    sbatch slurm_run_spts_wconf.sh data0${i}.cxi
done
