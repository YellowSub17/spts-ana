from pathlib import Path

# This locates the root relative to THIS config file
PROJ_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJ_ROOT / "data" 


def generate_full_paths(filenames, d):
    return list(map(lambda fname: f'{DATA_DIR}/newdata/{fname[:-4]}_analysis_d{d}/spts.cxi', filenames))



ps_data_ranges ={
    'ps20nm' :[f'data{x:05}.cxd' for x in range(1271, 1281)],
    'ps30nm' :[f'data{x:05}.cxd' for x in range(1263, 1270)],
    'ps40nm' : [f'data{x:05}.cxd' for x in range(1252, 1260)],
    'ps50nm' : [f'data{x:05}.cxd' for x in range(1245, 1251)],
    'groel':[f'data{x:05}.cxd' for x in range(1282, 1297)],
    'ferri':[f'data{x:05}.cxd' for x in range(1298, 1313)],
}



