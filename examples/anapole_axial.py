import sys
import os
FILE2_THIS_FILE = os.path.abspath(__file__)
sys.path.append(os.path.dirname(os.path.dirname(FILE2_THIS_FILE)))
from imsrg_toolkit.kshell_utils import KshellWavefunctionScript, KshellDensityScript, KshellToolkit
from imsrg_toolkit.utils import Utils
from imsrg_toolkit.settings import username
import numpy as np
import pandas as pd

#for O14
#emax=4 can request 10min with 1G and it is more than enough
#emax=6 can request 15min with 4G
#emax=8 can request 2:00 with 20G
#emax=10 can request 10:00:00 with 100G

VPNC_path = "/home/belleya/projects/rrg-holt/belleya/operators/"
imsrg_log_path = f"/home/belleya/projects/rrg-holt/belleya/results/imsrg_log/outputs/"
imsrg_error_path = f"/home/belleya/projects/rrg-holt/belleya/results/imsrg_log/errors/"
kshell_log_path = f"/home/belleya/projects/rrg-holt/belleya/results/kshell_log/outputs/"
kshell_error_path = f"/home/belleya/projects/rrg-holt/belleya/results/kshell_log/errors/"

##########PARAMETERS TO CHANGE BEFORE RUN###################
# emax = [1,2,3,4,5,6,7,8,9,10]
# time = ["00:10:00", "00:10:00","00:10:00","00:10:00","00:10:00","00:30:00","01:00:00", "02:00:00","06:00:00","08:00:00"]
# memory = ['10G','10G','10G','10G','10G','10G', "10G", '20G', "20G",'60G',"100G"]
# vs = 'p-shell'
# mass =  [9, 13, 14, 15]
# Nuclei = ["Be", "C", "N", "N"]
# states = ["1.5-1", "0.5-1", "1+1", "0.5-1"]
# emax = [2,3,4,5,6,7,8,9,10]
# time = ["00:10:00","00:10:00","00:10:00","00:10:00","00:30:00","01:00:00", "02:00:00","06:00:00","08:00:00"]
# memory = ['10G','10G','10G','10G','10G', "10G", '20G', "20G",'60G',"100G"]
# vs = 'sd-shell'
# mass =  [25, 29]
# Nuclei = ["Mg", "Si"]
# states = ["2.5+1", "0.5+1"]
# mass =  [29]
# Nuclei = ["Si"]
# states = ["0.5+1"]

# mass =  [113,115]
# vs = 'jj54'
# customvs = "Ni92,p1p1,p1p3,p0f5,p0g9,n2s1,n1d3,n0h11,n1f7"
# Nuclei = ["In","In"]
# states = ["4.5+1","4.5+1"]
mass =  [133]
vs = 'jj55'
customvs = None
Nuclei = ["Cs"]
states = ["3.5+1"]
emax = [12]
# time = ["3:00:00"]
# memory = ['200G']
emax = [6,8,10,12]
time = ['01:00:00','6:00:00', '24:00:00', "3-00:00:00"]
memory = ['5G','20G','60G','200G']

label = "N3LO_LNL"
file2b = "TwBME-HO_NN-only_N3LO_EM500_srg2.0_hw16_emax16_e2max32.me2j.gz"
file3b = "NO2B_ThBME_srg2.0_ramp46_N3LO_EM500_JJmax13_c1_-0.81_c3_-3.2_c4_5.4_cD_0.7_cE_-0.06_LNL2_650_500_IS_hw16from30_ms18_36_24.stream.bin"
# label = 'magic'
# file2b = "TwBME-HO_NN-only_N3LO_EM500_srg1.80_hw16_emax18_e2max36.me2j.gz"
# file3b = "NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_24.stream.bin"

###########################################################


def getNucl(Nucl, A):
  return f'{Nucl}{A}'
for Nucleus, A, state in zip(Nuclei, mass, states):
  for e, t, mem in zip(emax,time, memory):
    Nucl = getNucl(Nucleus, A)
    if 3*e < 24:
      E3max = 3*e
    else:
      E3max = 24
    for BCM in [1.0]:
      print(f"Running {Nucl} emax={e} E3max={E3max} BCM={BCM}")
      imsrg_params = {}
      imsrg_params['scratch_directory'] = f"/home/belleya/scratch/imsrg_work/"
      imsrg_params['output_directory_base'] = f"/home/belleya/projects/rrg-holt/belleya/results/"
      imsrg_params['file2b_directory'] = '/home/belleya/projects/def-holt/shared/me2j/'
      imsrg_params['emax'] = e
      imsrg_params['E3max'] = E3max
      imsrg_params['hw'] = 16
      imsrg_params['A'] = A
      imsrg_params['ref'] = Nucl
      imsrg_params['opnames'] = ['AnapoleAxial']
      imsrg_params['valence_space'] = vs # this is just a label when custom_valence_space is set
      if customvs:
        imsrg_params['custom_valence_space'] = customvs
      # imsrg_params['custom_valence_space'] = "Ni92,p1p1,p1p3,p0f5,p0g9,n2s1,n1d3,n0h11,n1f7"
      imsrg_params['label'] = label
      imsrg_params['BetaCM'] = BCM
      imsrg_params['hwBetaCM'] = 16
      imsrg_params['file2e1max'] = 16
      imsrg_params['file2e2max'] = 32
      imsrg_params['file2lmax'] = 16
      imsrg_params['file3e1max'] = 18
      imsrg_params['file3e2max'] = 36
      imsrg_params['file3e3max'] = 24
      imsrg_params['file3_precision'] = 'full'
      imsrg_params['run_cmd'] = """\
srun apptainer exec \\
  --bind /home/belleya/projects \
  --bind /home/belleya/scratch \
  /home/belleya/projects/rrg-holt/belleya/imsrg/imsrg.sif"""

      kshell_params = {}
      kshell_params['output_directory'] = f"/home/belleya/projects/rrg-holt/belleya/results/kshell/"
      kshell_params['scratch_directory'] = f"/home/belleya/scratch/kshell_anapole/"
      kshell_params['run_cmd'] = """srun"""
      kshell_params['header'] = (
      "#!/bin/bash\n"
      "module load openmpi\n"
  )


      imsrg_params['header'] = f"""#!/bin/bash
#SBATCH --job-name={imsrg_params['ref']}_emax{imsrg_params['emax']}_{imsrg_params['label']}_%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output={imsrg_log_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_{imsrg_params['label']}_%j.txt
#SBATCH --error={imsrg_error_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_{imsrg_params['label']}_%j.txt
#SBATCH --time={t}
#SBATCH --account=rrg-holt
#SBATCH --mem={mem}

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=24
"""

      kshell_params['header'] = f"""#!/bin/bash
#SBATCH --job-name=kshell_{Nucl}_emax{imsrg_params['emax']}_N3LO_%j
#SBATCH --nodes=30
#SBATCH --ntasks-per-node=1  
#SBATCH --mem=5000M
#SBATCH --cpus-per-task=48
#SBATCH --output={kshell_log_path}/kshell_{Nucl}_emax{imsrg_params['emax']}_{imsrg_params['label']}_%j.txt
#SBATCH --error={kshell_error_path}/kshell_{Nucl}_emax{imsrg_params['emax']}_{imsrg_params['label']}_%j.txt
#SBATCH --time=3:00:00
#SBATCH --account=rrg-holt
# ulimit -s unlimited
module load openmpi """


      header_expvals = f"""#SBATCH --output={kshell_log_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_{imsrg_params['label']}_eval_%j.txt
#SBATCH --error={kshell_error_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_{imsrg_params['label']}_eval_%j.txt
#SBATCH --account=rrg-holt"""

      imsrg_submit = Utils(Nucl, [state, state], imsrg_params, kshell_params)
      imsrg_submit.submit_all(file2b, file3b, f"{imsrg_submit.output_dir}/{imsrg_submit.filebase}_anapole_axial.csv", header_expvals = header_expvals, verbose=True, ops_rankJ=[1])
    # # count +=1 