import sys
# from imsrg_toolkit.imsrg import Imsrg
from imsrg_toolkit.kshell_utils import KshellWavefunctionScript, KshellDensityScript, KshellToolkit
from imsrg_toolkit.settings import username, ROOT_DIR
import numpy as np
import pandas as pd
import os
from pathlib import Path
import glob

def getNucl(Nucl, A):
  return f'{Nucl}{A}'


##########PARAMETERS TO CHANGE BEFORE RUN###################
kshell_log_path = f"/home/belleya/projects/rrg-holt/belleya/results/kshell_log/outputs/"
kshell_error_path = f"/home/belleya/projects/rrg-holt/belleya/results/kshell_log/errors/"
A =  133
Nucleus = "Cs"
Nucl = getNucl(Nucleus, A)
hw = 16
opnames = ['Rp2']


vs = 'jj55'
state = "1+1"

snt_files = glob.glob(f"/home/belleya/projects/rrg-holt/belleya/results/Cs133/magic/jj55_magic_Cs133_magnus_e12_E24_hw16_BetaCM1.0.snt")
for snt in snt_files:
  snt_list = snt.split("/")[-1].split("_")
  SampleID = snt_list[2]
  emax = snt_list[5][1:]
  E3max = snt_list[6][1:]
  # if int(emax) < 6: continue
  fn_snt = snt
  kshell_params = {}
  kshell_params['run_cmd'] = """\
mpirun -np $SLURM_NTASKS"""
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
#SBATCH --job-name=kshell_{Nucl}_emax{imsrg_params['emax']}_N3LO_%j
#SBATCH --nodes=30
#SBATCH --ntasks-per-node=1  
#SBATCH --mem=5000M
#SBATCH --cpus-per-task=48
#SBATCH --output={kshell_log_path}/kshell_{Nucl}_emax{imsrg_params['emax']}_{imsrg_params['label']}_%j.txt
#SBATCH --error={kshell_error_path}/kshell_{Nucl}_emax{imsrg_params['emax']}_{imsrg_params['label']}_%j.txt
#SBATCH --time=2:00:00
#SBATCH --account=rrg-holt
# ulimit -s unlimited
module load openmpi """



  output_dir = f"/home/submit/abelley/results/{Nucl}/SampleDelta/"
  filebase = f"{vs}_SampleDelta_{SampleID}_{Nucl}_magnus_e{emax}_E{E3max}_hw{hw}"
  fn_ops = [f"{output_dir}{filebase}_{op}.snt" for op in opnames]


  header_expvals = f"""#SBATCH --output={kshell_log_path}/{Nucl}_Sample_{SampleID}_e{emax}_E{E3max}_%j.txt
  #SBATCH --error={kshell_error_path}/{Nucl}_Sample_{SampleID}_e{emax}_E{E3max}_%j.txt"""

  kshell = KshellToolkit(fn_snt, Nucl, [state, state], Nucl_daughter=Nucl, submit_cmd="sbatch",  **kshell_params)
  kshell.submit_all(f"{output_dir}/{filebase}_Rp2_1+.csv", fn_ops = fn_ops, ops_rankJ=[0], verbose=True,  gen_partition=True, header=header_expvals)

