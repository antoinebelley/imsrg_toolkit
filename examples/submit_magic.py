import sys
# from imsrg_toolkit.imsrg import Imsrg
from imsrg_toolkit.kshell_utils import KshellWavefunctionScript, KshellDensityScript, KshellToolkit
from imsrg_toolkit.utils import Utils
import numpy as np
import pandas as pd

#for O14
#emax=4 can request 10min with 1G and it is more than enough
#emax=6 can request 15min with 4G
#emax=8 can request 2:00 with 20G
#emax=10 can request 10:00:00 with 100G

opfiles_path = "/work/submit/abelley/operators"

##########PARAMETERS TO CHANGE BEFORE RUN###################
# emax = [4,6,8,10]
# time = ["00:10:00", "00:30:00", "02:00:00","08:00:00"]
# memory = ['10G', "10G", "20G", "100G"]
# mass =  [22,23,24,25,27]
# Nucleus = "Al"
# vs = 'sd-shell'
emax = [10]
time = ["08:00:00"]
memory = ['50G']
mass =  [37]
Nucleus = "Ca"
vs = 'sd-shell'
state = "1.5+1"
file2b = "TwBME-HO_NN-only_N3LO_EM500_srg1.80_hw16_emax18_e2max36.me2j.gz"
file3b = "NO2B_half_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms16_32_28.stream.bin"
opnames = ['M1']
opfiles = [['/work/submit/abelley/operators/M1_2BC_bare_hw16_emax12_e2max24.me2j.gz',"M1_2BC"]]
###########################################################


def getNucl(Nucl, A):
  return f'{Nucl}{A}'

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
    imsrg_params['hw'] = 16
    imsrg_params['A'] = A
    imsrg_params['opnames'] = opnames
    imsrg_params['opfiles'] = opfiles
    imsrg_params['ref'] = Nucl
    imsrg_params['valence_space'] = vs # this is just a label when custom_valence_space is set
    imsrg_params['label'] = 'magic'
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
    kshell_params['header'] = f"""#!/bin/bash
#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=/work/submit/abelley/results/kshell_log/outputs/{imsrg_params['ref']}_emax{imsrg_params['emax']}_magic_%j.txt
#SBATCH --error=/work/submit/abelley/results/kshell_log/errors/{imsrg_params['ref']}_emax{imsrg_params['emax']}_magic_%j.txt
#SBATCH --time=10:00 """


    imsrg_params['header'] = f"""#!/bin/bash
#SBATCH --job-name={imsrg_params['ref']}_emax{imsrg_params['emax']}_magic_%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=/work/submit/abelley/results/imsrg_log/outputs/{imsrg_params['ref']}_emax{imsrg_params['emax']}_magic_%j.txt
#SBATCH --error=/work/submit/abelley/results/imsrg_log/errors/{imsrg_params['ref']}_emax{imsrg_params['emax']}_magic_%j.txt
#SBATCH --time={t}
#SBATCH --mem={m}

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=24
  """



    header_expvals = f"""#SBATCH --output=/work/submit/abelley/results/kshell_log/outputs/{imsrg_params['ref']}_emax{imsrg_params['emax']}_magic_eval_%j.txt
#SBATCH --error=/work/submit/abelley/results/kshell_log/errors/{imsrg_params['ref']}_emax{imsrg_params['emax']}_magic_eval_%j.txt"""

    imsrg_submit = Utils(Nucl, [state, state], imsrg_params, kshell_params)
    imsrg_submit.submit_all(file2b, file3b, f"{imsrg_submit.output_dir}/{imsrg_submit.filebase}_M1.csv", header_expvals = header_expvals, ops_rankJ = [1, 1])
    
  # count +=1 