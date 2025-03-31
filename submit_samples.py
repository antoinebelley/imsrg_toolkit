import sys
# from imsrg_toolkit.imsrg import Imsrg
from imsrg_toolkit.kshell_utils import KshellWavefunctionScript, KshellDensityScript, KshellToolkit
from imsrg_toolkit.utils import Utils
from imsrg_toolkit.settings import username, ROOT_DIR
import numpy as np
import pandas as pd
import os
from pathlib import PATH

#for O14
#emax=4 can request 2min with 256M and it is more than enough
#emax=6 can request 15min with 2G
#emax=8 can request 1:00 with 16G
#emax=10 can request 6:00 with 50G


##########PARAMETERS TO CHANGE BEFORE RUN###################
emax = [4,]
time = ["00:10:00"]
memory = ['1G',]
imsrg_log_path = f"/work/submit/{username}/results/imsrg_log/outputs/"
imsrg_error_path = f"/work/submit/{username}/results/imsrg_log/errors/"
kshell_log_path = f"/work/submit/{username}/results/kshell_log/outputs/"
kshell_error_path = f"/work/submit/{username}/results/kshell_log/errors/"
mass =  [16]
Nucleus = "O"

vs = 'sd-shell'
state = "+1"
num_samples = 1
###########################################################
# If the paths dont exist, create them
Path(kshell_log_path).mkdir(parents=True, exist_ok=True)
Path(kshell_error_path).mkdir(parents=True, exist_ok=True)
Path(imsrg_log_path).mkdir(parents=True, exist_ok=True)
Path(imsrg_error_path).mkdir(parents=True, exist_ok=True)


def getNucl(Nucl, A):
  return f'{Nucl}{A}'


LECs = ['Ct1S0pp','Ct1S0np','Ct1S0nn','Ct3S1','C1S0','C3P0','C1P1','C3P1','C3S1','CE1','C3P2','c1','c2','c3','c4','cD','cE']
df = pd.read_csv(f"{ROOT_DIR}/data/8000Samples.txt")

index = np.array(df.index)
rng = np.random.default_rng(seed=42)
index = rng.choice(index, num_samples, replace=False, shuffle=False)

for A in mass:
  for e, t, m in zip(emax,time, memory):
    Nucl = getNucl(Nucleus, A)
    if 3*e < 28:
      E3max = 3*e
    else:
      E3max = 28
    imsrg_params = {}
    imsrg_params['emax'] = e
    imsrg_params['E3max'] = E3max
    imsrg_params['hw'] = 10
    imsrg_params['A'] = A
    imsrg_params['opnames'] = ['Rp2']#, 'M1']
    imsrg_params['ref'] = Nucl
    imsrg_params['valence_space'] = vs # this is just a label when custom_valence_space is set
    imsrg_params['label'] = 'SampleDelta'
    imsrg_params['run_cmd'] = """\
  srun apptainer exec \\
    --bind /home/submit \\
    --bind /work/submit \\
    --bind /scratch/submit \\
    --bind /ceph/submit \\
    /work/submit/abelley/pyimsrg.sif """

    kshell_params = {}
    kshell_params['run_cmd'] = """\
  srun apptainer exec \\
    --bind /home/submit \\
    --bind /work/submit \\
    --bind /scratch/submit \\
    --bind /ceph/submit \\
    /work/submit/abelley/work/kshell/kshell.sif """
    

    for i in index:
      sample = df.iloc[i]
      SampleID = int(sample["SampleID"])
      weights = list(sample[LECs])

      imsrg_params['header'] = f"""#!/bin/bash
#SBATCH --job-name={imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output={imsrg_log_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_%j.txt
#SBATCH --error={imsrg_error_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_%j.txt
#SBATCH --time={t}
#SBATCH --mem={m}

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=24
"""
      kshell_params['header'] = f"""#!/bin/bash
#SBATCH --job-name=kshell_{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output={kshell_log_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}%j.txt
#SBATCH --error={kshell_error_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}%j.txt
#SBATCH --time=10:00 """



      header_expvals = f"""#SBATCH --output={kshell_log_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_eval_%j.txt
#SBATCH --error={kshell_error_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_eval_%j.txt"""

      imsrg_submit = Utils(Nucl, [state, state], imsrg_params, kshell_params, SampleID=SampleID)
      imsrg_submit.submit_all_combine_delta(weights, SampleID, f"{imsrg_submit.output_dir}/{imsrg_submit.filebase}_R2p.csv", header_expvals = header_expvals, verbose=True)
