import sys
import os

FILE2_THIS_FILE = os.path.abspath(__file__)
sys.path.append(os.path.dirname(os.path.dirname(FILE2_THIS_FILE)))
from imsrg_toolkit.utils import Utils
from imsrg_toolkit.settings import username, ROOT_DIR
import numpy as np
import pandas as pd
from pathlib import Path

# Use 'python3 submit_multiple_samples.py --dry-run' to generate everything
# without submitting
# DRY_RUN = '--dry-run' in sys.argv
# # Limit how many array tasks run simultaneously (None = no limit)
# MAX_CONCURRENT = None

LECs = ['Ct1S0pp','Ct1S0np','Ct1S0nn','Ct3S1','C1S0','C3P0','C1P1','C3P1','C3S1','CE1','C3P2','c1','c2','c3','c4','cD','cE']
# df = pd.read_csv(f"{ROOT_DIR}/data/8000Samples.txt")
df = pd.read_csv(f"{ROOT_DIR}/data/34Samples.csv")
# samples_run = [500, 1429, 2125,2411,2618,3105,3260,3373,3448,3813,3895,4117]
samples_run = [493,500,606,774,1141,1429,2125,2411,2618,3260,3319,3373,3448,3621,3813,4117]
# index = np.array([1728])
# sample1 = [0.356199,-0.355142,-0.349985,-0.248157,2.78763,0.433097,-0.086705,-1.299329,0.788,0.520089,-0.846265,-0.736674,-0.606916,-0.47898,0.930314,-0.588286,-0.068944]
# sample13 = [-0.349296,-0.34783,-0.345087,-0.251179,2.775803,0.51126,-0.014905,-0.995633,0.932483,0.358223,-0.787001,-0.707744,-0.441933,-0.706367,0.998484,0.493866,0.372671]
# samples = [sample1, sample13]
# Nucl = "Al22"
# emax = [4,6,8,10]
# time = ["00:10:00","00:30:00", "02:00:00","08:00:00"]
# memory = ['10G', "10G", '20G',"100G"]
# emax = [12,14]
# time = ["1-00:00:00","2-12:00:00"]
# memory = ['200G','502G']
emax = [14]
time = ["3-00:00:00"]
memory = ['502G']
# emax = [6]
# time = ["00:30:00"]
# memory = ['10G']
# emax = [14]
# time = ["3-00:00:00"]
# memory = ['502G']
# emax = [12]
# time = ["1-00:00:00"]
# memory = ['180G']
samples_numbers = [34]
# As = [22,23,24,25,26,27]
# states = ['4+1',"2.5+1", "4+1", "2.5+1", "5+1", "2.5+1"]
# As = [26,27]
states = ["0+1"]
As = [136]

imsrg_log_path = f"/home/belleya/projects/rrg-holt/belleya/results/imsrg_log/outputs/"
imsrg_error_path = f"/home/belleya/projects/rrg-holt/belleya/results/imsrg_log/errors/"
kshell_log_path = f"/home/belleya/projects/rrg-holt/belleya/results/kshell_log/outputs/"
kshell_error_path = f"/home/belleya/projects/rrg-holt/belleya/results/kshell_log/errors/"


Path(kshell_log_path).mkdir(parents=True, exist_ok=True)
Path(kshell_error_path).mkdir(parents=True, exist_ok=True)
Path(imsrg_log_path).mkdir(parents=True, exist_ok=True)
Path(imsrg_error_path).mkdir(parents=True, exist_ok=True)

for A, state in zip(As,states):
  Nucl = f"Xe{A}" 
  Nucl_daughter = f"Ba{A}"
  # Nucl = f"Ge{A}" 
  # Nucl_daughter = f"Se{A}"
  # Nucl = f"Te{A}" 
  # Nucl_daughter = f"Xe{A}"
  imsrg_params = {}
  imsrg_params['imsrg_scripts_directory'] = f"/home/belleya/scratch/imsrg_work/"
  imsrg_params['output_directory_base'] = f"/home/belleya/projects/rrg-holt/belleya/results/"
  imsrg_params['E3max'] = 28
  imsrg_params['hw'] = 10
  # imsrg_params['BetaCM'] = 4
  imsrg_params['A'] = A
#   imsrg_params['opnames'] = ['Rp2']
  imsrg_params['opnames_decay'] = ['M0nu_GT_13.06_none','M0nu_F_13.06_none','M0nu_T_13.06_none','M0nu_C_394_4']
  # imsrg_params['opnames_decay'] = ['M0nu_GT_9.61_none','M0nu_F_9.61_none','M0nu_T_9.61_none','M0nu_C_394_4']
  # imsrg_params['opnames_decay'] = ['M0nu_GT_12.76_none','M0nu_F_12.76_none','M0nu_T_12.76_none','M0nu_C_394_4']
  # imsrg_params['opfiles'] = opfiles = [['/work/submit/abelley/operators/M1_2BC_bare_hw10_emax12_e2max24.me2j.gz',"M1_2BC"]]
  imsrg_params['ref'] = Nucl
#   imsrg_params['valence_space'] = 'p-shell' # this is just a label when custom_valence_space is set
  # imsrg_params['valence_space'] = 'pf5g9' # this is just a label when custom_valence_space is set
  # imsrg_params['custom_valence_space'] = "Ni56,p0f5,n0f5,p1p3,n1p3,p1p1,n1p1,p0g9,n0g9"
  imsrg_params['valence_space'] = 'sdg' # this is just a label when custom_valence_space is set
  imsrg_params['custom_valence_space'] = 'Sn100,p0g7,n0g7,p1d3,n1d3,p1d5,n1d5,p2s1,n2s1,p0h11,n0h11'
  imsrg_params['label'] = 'Sample34'
  imsrg_params['denominator_delta'] = 10
  imsrg_params['denominator_delta_orbit'] = 'all'
  imsrg_params['run_cmd'] = """\
srun apptainer exec \\
  --bind /home/belleya/projects \
  --bind /home/belleya/scratch \
  /home/belleya/projects/rrg-holt/belleya/imsrg/imsrg.sif"""
  kshell_params = {}
  kshell_params['output_directory'] = f"/home/belleya/projects/rrg-holt/belleya/results/kshell/"
  kshell_params['scratch_directory'] = f"/home/belleya/scratch/kshell_works/"
  kshell_params['run_cmd'] = """srun"""
  kshell_params['header'] = (
      "#!/bin/bash\n"
      "module load openmpi\n"
  )

  for e, t, mem in zip(emax, time, memory):
    imsrg_params['emax'] = e 
    index = np.array(df.index)
    # rng = np.random.default_rng(seed=42)
    # index = rng.choice(index, num_samples, replace=False, shuffle=False)

    for i in index:
      sample = df.iloc[i]
      SampleID = int(sample["SampleID"])
      if SampleID not in samples_run: 
        continue
      weights = list(sample[LECs])
      if e==14:
        imsrg_params['scratch_directory'] = f'/home/belleya/scratch/imsrg_{Nucl}_Sample_{SampleID}_e{e}/'

    # for i, sample in enumerate(samples):
    #   sample = df.iloc[i]
    #   SampleID = i
    #   weights = sample

      imsrg_params['header'] = f"""#!/bin/bash
#SBATCH --job-name={imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output={imsrg_log_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_%j.txt
#SBATCH --error={imsrg_error_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_%j.txt
#SBATCH --time={t}
#SBATCH --account=rrg-holt
#SBATCH --mem={mem}

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=24
"""
      kshell_params['header'] = f"""#!/bin/bash
#SBATCH --job-name=kshell_{Nucl}_emax{imsrg_params['emax']}_Sample{SampleID}_%j
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1  
#SBATCH --mem-per-cpu=1000M
#SBATCH --cpus-per-task=24
#SBATCH --output={kshell_log_path}/kshell_{Nucl}_emax{imsrg_params['emax']}_Sample{SampleID}_%j.txt
#SBATCH --error={kshell_error_path}/kshell_{Nucl}_emax{imsrg_params['emax']}_{SampleID}_%j.txt
#SBATCH --time=0:10:00
#SBATCH --account=rrg-holt
# ulimit -s unlimited
module load openmpi """
      kshell_params['header_daughter'] = f"""#!/bin/bash
#SBATCH --job-name=kshell_{Nucl_daughter}_emax{imsrg_params['emax']}_Sample{SampleID}_%j
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1  
#SBATCH --mem-per-cpu=2000M
#SBATCH --cpus-per-task=24
#SBATCH --output={kshell_log_path}/kshell_{Nucl_daughter}_emax{imsrg_params['emax']}_Sample{SampleID}_%j.txt
#SBATCH --error={kshell_error_path}/kshell_{Nucl_daughter}_emax{imsrg_params['emax']}_{SampleID}_%j.txt
#SBATCH --time=1:00:00
#SBATCH --account=rrg-holt
# ulimit -s unlimited
module load openmpi """
      kshell_params['header_density'] = f"""#!/bin/bash
#SBATCH --job-name=density_{Nucl}_{Nucl_daughter}_emax{imsrg_params['emax']}_Sample{SampleID}_%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
#SBATCH --mem=20G
#SBATCH --time=0:05:00
#SBATCH --account=rrg-holt
# ulimit -s unlimited
module load openmpi """

      imsrg_params['SampleID'] =  SampleID
      imsrg_params['LECs'] =  weights
      header_expvals = f"""#SBATCH --output={kshell_log_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_eval_%j.txt
#SBATCH --error={kshell_error_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_eval_%j.txt
#SBATCH --account=rrg-holt"""

      imsrg_submit = Utils(Nucl, [state,state], imsrg_params, kshell_params, Nucl_daughter=Nucl_daughter)
      imsrg_submit.submit_all_combine_delta(f"{imsrg_submit.output_dir}/{imsrg_submit.filebase}_decay_test.csv", 
                                            header_expvals = header_expvals, ops_rankZ_decay=[2,2,2,2], verbose=True)
  
      # exit(0)
      # fn_ops = [f"{imsrg_submit.output_dir}{imsrg_submit.filebase}_{op}.snt" for op in imsrg_submit.opnames_decay]
      # tmp = [f"{imsrg_submit.output_dir}{imsrg_submit.filebase}_{op[1]}.snt" for op in imsrg_submit.opfiles]
      # fn_ops.extend(tmp)
      # try:
      #   # imsrg_submit.kshell.submit_all(f"{imsrg_submit.output_dir}/{imsrg_submit.filebase}_ops.csv", fn_ops_decay=fn_ops,  ops_rankZ_decay=[2,2,2,2],  verbose=True, header=header_expvals, gen_partition=True)
      #   imsrg_submit.kshell.submit_from_density(f"{imsrg_submit.output_dir}/{imsrg_submit.filebase}_ops.csv", fn_ops_decay=fn_ops,  ops_rankZ_decay=[2,2,2,2],  verbose=True, header=header_expvals, gen_partition=True)

      # except FileNotFoundError:
      #   print(f"Could find sample {SampleID} snt file.")
      #   continue
      
