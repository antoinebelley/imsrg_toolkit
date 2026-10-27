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
df = pd.read_csv(f"{ROOT_DIR}/data/8000Samples.txt")


# Nucl = "Al22"
emax = [4,6,8,10]
time = ["00:10:00","00:30:00", "02:00:00","08:00:00"]
memory = ['10G', "10G", '20G',"100G"]
# emax = [4]
# time = ["00:10:00"]
# memory = ['10G']
samples_numbers = [400,200,100,50]
# As = [22,23,24,25,26,27]
# states = ['4+1',"2.5+1", "4+1", "2.5+1", "5+1", "2.5+1"]
# As = [26,27]
# states = ["0+1", "0.5+1", "0+1", "1.5+1", "0+1", "0.5+1", "0+1", "1.5-1", "0+1"]
# As = [26,27,28,29,30,31,32,33,34]
# states = ["1.5+1"]
# As = [29]
states = ["0+1"]
As = [32]

imsrg_log_path = f"/home/belleya/projects/rrg-holt/belleya/results/imsrg_log/outputs/"
imsrg_error_path = f"/home/belleya/projects/rrg-holt/belleya/results/imsrg_log/errors/"
kshell_log_path = f"/home/belleya/projects/rrg-holt/belleya/results/kshell_log/outputs/"
kshell_error_path = f"/home/belleya/projects/rrg-holt/belleya/results/kshell_log/errors/"

Path(kshell_log_path).mkdir(parents=True, exist_ok=True)
Path(kshell_error_path).mkdir(parents=True, exist_ok=True)
Path(imsrg_log_path).mkdir(parents=True, exist_ok=True)
Path(imsrg_error_path).mkdir(parents=True, exist_ok=True)

for A, state in zip(As,states):
  Nucl = f"Mg{A}" 
#   Nucl_daughter = f"Be{A}"
  imsrg_params = {}
  imsrg_params['imsrg_scripts_directory'] = f"/home/belleya/scratch/imsrg_work/"
  imsrg_params['output_directory_base'] = f"/home/belleya/projects/rrg-holt/belleya/results/"
  imsrg_params['E3max'] = 28
  imsrg_params['hw'] = 10
  imsrg_params['BetaCM'] = 3
  imsrg_params['A'] = A
  imsrg_params['opnames'] = ['Rp2']
#   imsrg_params['opnames_decay'] = ['M0nu_GT_2.74_none']
  # imsrg_params['opfiles'] = opfiles = [['/work/submit/abelley/operators/M1_2BC_bare_hw10_emax12_e2max24.me2j.gz',"M1_2BC"]]
  imsrg_params['ref'] = Nucl
#   imsrg_params['valence_space'] = 'p-shell' # this is just a label when custom_valence_space is set
  imsrg_params['valence_space'] = 'PsdNsdfp-shell' # this is just a label when custom_valence_space is set
  imsrg_params['custom_valence_space'] = "O16,p0d5,p0d3,p1s1,n0d5,n0d3,n1s1,n0f7,n1p3"
  imsrg_params['label'] = 'SampleDelta'
  imsrg_params['denominator_delta'] = 10
  imsrg_params['denominator_delta_orbit'] = 'all'
  imsrg_params['run_cmd'] = """\
srun apptainer exec \\
  --bind /home/belleya/projects \
  --bind /home/belleya/scratch \
  /home/belleya/projects/rrg-holt/belleya/imsrg/imsrg.sif"""
  kshell_params = {}
  kshell_params['scratch_directory'] = f"/work/submit/{username}/work/test_mpi/"
  kshell_params['run_cmd'] = """\
mpirun -np $SLURM_NTASKS"""
  kshell_params = {}
  kshell_params['output_directory'] = f"/home/belleya/projects/rrg-holt/belleya/results/kshell/"
  kshell_params['scratch_directory'] = f"/home/belleya/scratch/kshell_works/"
  kshell_params['run_cmd'] = """srun"""
  kshell_params['header'] = (
      "#!/bin/bash\n"
      "module load openmpi\n"
  )
  for e, t, mem, num_samples in zip(emax, time, memory, samples_numbers):
    imsrg_params['emax'] = e
    index = np.array(df.index)
    rng = np.random.default_rng(seed=42)
    index = rng.choice(index, num_samples, replace=False, shuffle=False)

    for i in index:
      sample = df.iloc[i]
      SampleID = int(sample["SampleID"])
      weights = list(sample[LECs])

#   for i, sample in enumerate(samples):
#     sample = df.iloc[i]
#     SampleID = i
#     weights = sample

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

      kshell_params['header_density'] = f"""#!/bin/bash
#SBATCH --job-name=density_{Nucl}_{Nucl}_emax{imsrg_params['emax']}_Sample{SampleID}_%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
#SBATCH --mem=20G
#SBATCH --time=0:05:00
#SBATCH --account=rrg-holt
# ulimit -s unlimited
module load openmpi """


      imsrg_params['SampleID'] =  SampleID
      imsrg_params['LECs'] =  weights
      header_expvals = f"""#SBATCH --output={kshell_log_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_{imsrg_params['label']}_eval_%j.txt
#SBATCH --error={kshell_error_path}/{imsrg_params['ref']}_emax{imsrg_params['emax']}_{imsrg_params['label']}_eval_%j.txt
#SBATCH --account=rrg-holt"""

      imsrg_submit = Utils(Nucl, [state,state], imsrg_params, kshell_params)
      imsrg_submit.submit_all_combine_delta(f"{imsrg_submit.output_dir}/{imsrg_submit.filebase}_radii.csv", 
                                            header_expvals = header_expvals, verbose=True)

      # fn_ops = [f"{imsrg_submit.output_dir}{imsrg_submit.filebase}_{op}.snt" for op in imsrg_submit.opnames]
      # # tmp = [f"{imsrg_submit.output_dir}{imsrg_submit.filebase}_{op[1]}.snt" for op in imsrg_submit.opfiles]
      # # fn_ops.extend(tmp)
      # try:
      #   imsrg_submit.kshell.submit_all(f"{imsrg_submit.output_dir}/{imsrg_submit.filebase}_radii.csv", fn_ops,  ops_rankJ = [0],  verbose=True, header=header_expvals, gen_partition=True)
      # except:
      #   print(f"Error submitting {imsrg_submit.output_dir}/{imsrg_submit.filebase}_radii.csv")
      #   continue
      # count +=1
      # if count > 30:
      #   exit(0)
