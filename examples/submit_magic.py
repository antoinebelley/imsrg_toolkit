import sys
import sys
import os
FILE2_THIS_FILE = os.path.abspath(__file__)
sys.path.append(os.path.dirname(os.path.dirname(FILE2_THIS_FILE)))
# from imsrg_toolkit.imsrg import Imsrg
# from examples.submit_kshell_only import SampleID
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

opfiles_path = "/work/submit/abelley/operators"
imsrg_log_path = f"/home/belleya/projects/rrg-holt/belleya/results/imsrg_log/outputs/"
imsrg_error_path = f"/home/belleya/projects/rrg-holt/belleya/results/imsrg_log/errors/"
kshell_log_path = f"/home/belleya/projects/rrg-holt/belleya/results/kshell_log/outputs/"
kshell_error_path = f"/home/belleya/projects/rrg-holt/belleya/results/kshell_log/errors/"

##########PARAMETERS TO CHANGE BEFORE RUN###################
# emax = [4,6,8,10]
# time = ["00:10:00", "01:00:00", "03:00:00","08:00:00"]
# memory = ['10G', "10G", "30G","100G"]

# emax = [6,8,10]
# time = ["01:00:00", "03:00:00","08:00:00"]
# memory = ["10G", "30G","100G"]

# emax = [10]
# time = ["10:00:00"]
# memory = ["100G"]

# Nucleus = "Al"
# vs = 'sd-shell'
# mass = [22,23,24,25,26,27]
# states = ["4+1","2.5+1","4+1","2.5+1","5+1","2.5+1"]

# Nucleus = 'Mg'
# states = ["0+1", "0.5+1", "0+1", "1.5+1", "0+1", "0.5+1", "0+1", "1.5-1", "0+1"]
# As = [26,27,28,29,30,31,32,33,34]

emax = [6]
time = ["01:30:00"]
memory = ['10G']

Nucleus = 'Mg'
states = ["1.5-1"]#, "0+1", "0.5+1", "0+1", "1.5+1", "0+1", "0.5+1", "0+1", "1.5-1", "0+1"]
As = [33]#26,27,28,29,30,31,32,33,34]

# Nucleus = 'Si'
# states = ["0+1", "3.5-1"]
# As = [34,35]

# vs = '0hw-shell'

# Nucleus = 'C'
# states = ["0+1"]
# As = [12]
# Label = "magic"
# file2b = "TwBME-HO_NN-only_N3LO_EM500_srg1.80_hw16_emax18_e2max36.me2j.gz"
# file3b = "NO2B_half_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms16_32_28.stream.bin"
# Label = "DeltaGo"
# file2b = "TwBME-HO_NN-only_DNNLOgo_bare_hw12_emax18_e2max36.me2j.gz"
# file3b = "NO2B_half_ThBME_DNNLOgo_3NFJmax15_IS_hw12_ms16_32_28.stream.bin"
Label = "Texas"
file2b = "TwBME-HO_NN-only_N3LO_Texas_bare_hw16_emax18_e2max36.me2j.gz"
file3b = "NO2B_half_ThBME_Texas_3NFJmax15_IS_hw16_ms18_36_28.stream.bin"
# opnames = ['Eccentricity_2_0']
opnames = ['Rp2']
# opfiles = [['/work/submit/abelley/operators/M1_2BC_bare_hw16_emax12_e2max24.me2j.gz',"M1_2BC"]]
###########################################################


def getNucl(Nucl, A):
  return f'{Nucl}{A}'

for A, state in zip(As,states):
  Nucl = getNucl(Nucleus, A)
  for e, t, mem in zip(emax,time, memory):
    if 3*e < 28:
      E3max = 3*e
    else:
      E3max = 28
    imsrg_params = {}
    # imsrg_params['approx'] = 'imsrg3f2'
    imsrg_params['scratch_directory'] = f"/home/belleya/scratch/imsrg_work/"
    imsrg_params['output_directory_base'] = f"/home/belleya/projects/rrg-holt/belleya/results/"
    imsrg_params['file2b_directory'] = '/home/belleya/projects/def-holt/shared/me2j/'
    imsrg_params['file2e1max'] = 18
    imsrg_params['file2e2max'] = 36
    imsrg_params['file2lmax'] = 18
    if Label == "Texas":
      imsrg_params['file2b_directory'] = '/home/belleya/projects/rrg-holt/shared/me_texas/'
      imsrg_params['file3b_directory'] = '/home/belleya/projects/rrg-holt/shared/me_texas/'
      imsrg_params['file3e1max'] = 18
      imsrg_params['file3e2max'] = 36
    imsrg_params['emax'] = e
    imsrg_params['E3max'] = E3max
    imsrg_params['hw'] = 16
    if Label == "DeltaGo":
      imsrg_params['hw'] = 12
    imsrg_params['A'] = A
    imsrg_params['opnames'] = opnames
    # imsrg_params['opfiles'] = opfilqes
    imsrg_params['ref'] = Nucl
    # imsrg_params['valence_space'] = vs # this is just a label when custom_valence_space is set
    imsrg_params['valence_space'] = 'PsdNsdfp-shell' # this is just a label when custom_valence_space is set
    imsrg_params['custom_valence_space'] = "O16,p0d5,p0d3,p1s1,n0d5,n0d3,n1s1,n0f7,n1p3"
    imsrg_params['BetaCM'] = 3
    imsrg_params['denominator_delta'] = 10
    imsrg_params['denominator_delta_orbit'] = 'all'
    imsrg_params['label'] = Label
    imsrg_params['run_cmd'] = """\
srun apptainer exec \\
  --bind /home/belleya/projects \
  --bind /home/belleya/scratch \
  /home/belleya/projects/rrg-holt/belleya/imsrg/imsrg.sif"""

    kshell_params = {}
    kshell_params['output_directory'] = f"/home/belleya/projects/rrg-holt/belleya/results/kshell/"
    kshell_params['scratch_directory'] = f"/home/belleya/scratch/kshell_magnesium/"
    kshell_params['run_cmd'] = """srun"""
    kshell_params['header'] = (
    "#!/bin/bash\n"
    "module load openmpi\n"
)

#     kshell_params['header'] = f"""#!/bin/bash
# #SBATCH --job-name=kshell_{Nucl}_emax{emax}_magic_%j
# #SBATCH --nodes=4
# #SBATCH --ntasks-per-node=1
# #SBATCH --ntasks=4
# #SBATCH --cpus-per-task=2
# #SBATCH --mem-per-cpu=100
# #SBATCH --time=0-01:00
# #SBATCH --output={kshell_log_path}/{Nucl}_magic_e{emax}_E{E3max}_%j.txt
# #SBATCH --error={kshell_error_path}/{Nucl}_magic_e{emax}_E{E3max}_%j.txt
# # ulimit -s unlimited
# module load mpi"""
    kshell_params['header'] = f"""#!/bin/bash
#SBATCH --job-name=kshell_{Nucl}_emax{imsrg_params['emax']}_magic_%j
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1  
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=/work/submit/abelley/results/kshell_log/outputs/{imsrg_params['ref']}_emax{imsrg_params['emax']}_magic_%j.txt
#SBATCH --error=/work/submit/abelley/results/kshell_log/errors/{imsrg_params['ref']}_emax{imsrg_params['emax']}_magic_%j.txt
#SBATCH --time=30:00
# ulimit -s unlimited
module load mpi """


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
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=1  
#SBATCH --mem=2000M
#SBATCH --cpus-per-task=48
#SBATCH --output={kshell_log_path}/kshell_{Nucl}_emax{imsrg_params['emax']}_{imsrg_params['label']}_%j.txt
#SBATCH --error={kshell_error_path}/kshell_{Nucl}_emax{imsrg_params['emax']}_{imsrg_params['label']}_%j.txt
#SBATCH --time=2:00:00
#SBATCH --account=rrg-holt
# ulimit -s unlimited
module load openmpi """



    header_expvals = f"""#SBATCH --output={kshell_log_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_{imsrg_params['label']}_eval_%j.txt
#SBATCH --error={kshell_error_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_{imsrg_params['label']}_eval_%j.txt
#SBATCH --account=rrg-holt"""
    imsrg_submit = Utils(Nucl, [state, state], imsrg_params, kshell_params)
    imsrg_submit.submit_all(file2b, file3b, f"{imsrg_submit.output_dir}/{imsrg_submit.filebase}_radii.csv", header_expvals = header_expvals, verbose=True)
    # fn_ops = imsrg_submit.gen_oplist()
    # imsrg_submit.kshell.submit_all(f"{imsrg_submit.output_dir}/{imsrg_submit.filebase}_radii.csv", fn_ops, header = header_expvals, verbose=True)
  # count +=1 