import sys
# from imsrg_toolkit.imsrg import Imsrg
from imsrg_toolkit.kshell_utils import KshellWavefunctionScript, KshellDensityScript, KshellToolkit
from imsrg_toolkit.settings import username, ROOT_DIR
import numpy as np
import pandas as pd
import os
from pathlib import Path

def getNucl(Nucl, A):
  return f'{Nucl}{A}'


##########PARAMETERS TO CHANGE BEFORE RUN###################
kshell_log_path = f"/work/submit/{username}/results/kshell_log/outputs/"
kshell_error_path = f"/work/submit/{username}/results/kshell_log/errors/"
A =  24
Nucleus = "Al"
Nucl = getNucl(Nucleus, A)
e = 4
E = 12
hw = 16
opnames = ['M1']
opfiles = [['/work/submit/abelley/operators/M1_2BC_bare_hw16_emax12_e2max24.me2j.gz',"M1_2BC"]]

vs = 'sd-shell'
state = "4+1"

fn_snt = f"/home/submit/abelley/results/{Nucl}/magic/sd-shell_magic_{Nucl}_magnus_e{e}_E{E}_hw{hw}.snt"
# fn_snt = f"/work/submit/abelley/kshell-20230714-ver4/snt/usdb.snt"

kshell_params = {}
kshell_params['run_cmd'] = """\
  mpirun -np $SLURM_NTASKS """
# kshell_params['run_cmd'] = """\
#    srun apptainer exec \\
#      --bind /home/submit \\
#      --bind /work/submit \\
#      --bind /scratch/submit \\
#      --bind /ceph/submit \\
#      /work/submit/abelley/work/kshell/kshell.sif"""
###########################################################
# If the paths dont exist, create them
Path(kshell_log_path).mkdir(parents=True, exist_ok=True)
Path(kshell_error_path).mkdir(parents=True, exist_ok=True)
kshell_params['header'] = f"""#!/bin/bash
#SBATCH --job-name=test_kshell_Al24_emax{e}_magic_%j
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=0-01:00
# ulimit -s unlimited
module load mpi"""
# kshell_params['header'] = f"""#!/bin/bash
# #SBATCH --job-name=test_kshell_Al24_emax{e}_magic_%j
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=10
# #SBATCH --output=/work/submit/abelley/results/kshell_log/outputs/{Nucl}_emax{e}_magic_%j.txt
# #SBATCH --error=/work/submit/abelley/results/kshell_log/errors/{Nucl}_emax{e}_magic_%j.txt
# #SBATCH --time=10:00
# """
# kshell_params['header'] = f"""#!/bin/bash
# #SBATCH --job-name=test_kshell_Al24_emax{e}_magic_%j
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=10
# #SBATCH --time=10:00
# """


output_dir = f"/home/submit/abelley/results/{Nucl}/magic/"
filebase = f"{vs}_magic_{Nucl}_magnus_e{e}_E{E}_hw{hw}"
fn_ops = [f"{output_dir}{filebase}_{op}.snt" for op in opnames]
# tmp = [f"{output_dir}{filebase}_{op[1]}.snt" for op in opfiles]
# fn_ops.extend(tmp)


header_expvals = f"""#SBATCH --output={kshell_log_path}/testmpi_eval_%j.txt
#SBATCH --error={kshell_error_path}/testmpi_eval_%j.txt"""

kshell = KshellToolkit(fn_snt, Nucl, [state, state], Nucl_daughter=Nucl, submit_cmd="sbatch",  **kshell_params)
kshell.submit_all(f"{output_dir}/{filebase}_M1_test.csv", fn_ops = fn_ops, ops_rankJ=[1,1], verbose=True,  gen_partition=True, header=header_expvals)

