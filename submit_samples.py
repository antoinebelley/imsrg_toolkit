import sys
import os

FILE2_THIS_FILE = os.path.abspath(__file__)
sys.path.append(os.path.dirname(FILE2_THIS_FILE))
from imsrg_toolkit.utils import Utils
from imsrg_toolkit.job_array import JobArrayChain
from imsrg_toolkit.settings import username, ROOT_DIR
import numpy as np
import pandas as pd
from pathlib import Path

#for O14
#emax=4 can request 2min with 256M and it is more than enough
#emax=6 can request 15min with 2G
#emax=8 can request 1:00 with 16G
#emax=10 can request 6:00 with 50G

# Use 'python3 submit_samples.py --dry-run' to generate everything without submitting
DRY_RUN = '--dry-run' in sys.argv
# Limit how many array tasks run simultaneously (None = no limit)
MAX_CONCURRENT = None

##########PARAMETERS TO CHANGE BEFORE RUN###################
# emax = [8]
# time = ["3:00:00"]
# memory = ["20G"]
emax = [10]
time = ["12:00:00"]
memory = ["100G"]
imsrg_log_path = f"/work/submit/{username}/results/imsrg_log/outputs/"
imsrg_error_path = f"/work/submit/{username}/results/imsrg_log/errors/"
kshell_log_path = f"/work/submit/{username}/results/kshell_log/outputs/"
kshell_error_path = f"/work/submit/{username}/results/kshell_log/errors/"
array_script_dir = f"/work/submit/{username}/work/job_arrays/"
mass =  [27]
Nucleus = "Al"

vs = 'sd-shell'
state = "2.5+1"
num_samples = 100
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
    # Resources come from the job-array headers below; this only sets
    # up the environment of each task.
    imsrg_params['header'] = (
        "#!/bin/bash\n"
        "export OMP_NUM_THREADS=24\n"
    )
    imsrg_params['run_cmd'] = """\
  srun apptainer exec \\
    --bind /home/submit \\
    --bind /work/submit \\
    --bind /scratch/submit \\
    --bind /ceph/submit \\
    /work/submit/abelley/pyimsrg.sif """

    kshell_params = {}
    kshell_params['scratch_directory'] = f"/work/submit/{username}/work/test3/"
    kshell_params['header'] = "#!/bin/bash\n"
    kshell_params['run_cmd'] = """\
  srun apptainer exec \\
    --bind /home/submit \\
    --bind /work/submit \\
    --bind /scratch/submit \\
    --bind /ceph/submit \\
    /work/submit/abelley/work/kshell/kshell.sif """

    # One job array per stage (one task per sample), chained with
    # aftercorr dependencies, instead of 4 individual jobs per sample.
    imsrg_array_header = (
        f"#SBATCH --job-name={Nucl}_e{e}_imsrg\n"
        f"#SBATCH --nodes=1\n"
        f"#SBATCH --ntasks=1\n"
        f"#SBATCH --cpus-per-task=24\n"
        f"#SBATCH --output={imsrg_log_path}/{Nucl}_emax{e}_imsrg_%A_%a.txt\n"
        f"#SBATCH --error={imsrg_error_path}/{Nucl}_emax{e}_imsrg_%A_%a.txt\n"
        f"#SBATCH --time={t}\n"
        f"#SBATCH --mem={m}\n"
    )
    diag_array_header = (
        f"#SBATCH --job-name=kshell_{Nucl}_e{e}_diag\n"
        f"#SBATCH --nodes=1\n"
        f"#SBATCH --ntasks=1\n"
        f"#SBATCH --cpus-per-task=10\n"
        f"#SBATCH --output={kshell_log_path}/{Nucl}_emax{e}_diag_%A_%a.txt\n"
        f"#SBATCH --error={kshell_error_path}/{Nucl}_emax{e}_diag_%A_%a.txt\n"
        f"#SBATCH --time=00:30:00\n"
    )
    density_array_header = (
        f"#SBATCH --job-name=kshell_{Nucl}_e{e}_density\n"
        f"#SBATCH --nodes=1\n"
        f"#SBATCH --ntasks=1\n"
        f"#SBATCH --cpus-per-task=10\n"
        f"#SBATCH --output={kshell_log_path}/{Nucl}_emax{e}_density_%A_%a.txt\n"
        f"#SBATCH --error={kshell_error_path}/{Nucl}_emax{e}_density_%A_%a.txt\n"
        f"#SBATCH --time=00:30:00\n"
    )
    expvals_array_header = (
        f"#SBATCH --job-name={Nucl}_e{e}_expvals\n"
        f"#SBATCH --output={kshell_log_path}/{Nucl}_emax{e}_expvals_%A_%a.txt\n"
        f"#SBATCH --error={kshell_error_path}/{Nucl}_emax{e}_expvals_%A_%a.txt\n"
    )

    chain = JobArrayChain(f"{Nucl}_e{e}_hw{imsrg_params['hw']}", array_script_dir)
    kshell_workdir = kshell_params['scratch_directory']
    stages = {
        'imsrg': chain.new_stage('imsrg', imsrg_array_header, workdir=kshell_workdir),
        'diag': chain.new_stage('diag', diag_array_header, workdir=kshell_workdir),
        'density': chain.new_stage('density', density_array_header, workdir=kshell_workdir),
        'expvals': chain.new_stage('expvals', expvals_array_header, workdir=kshell_workdir),
    }

    for i in index:
      sample = df.iloc[i]
      SampleID = int(sample["SampleID"])
      weights = list(sample[LECs])
      print("generated scripts for sample #", SampleID)

      imsrg_params['SampleID'] = SampleID
      imsrg_params['LECs'] = weights

      imsrg_submit = Utils(Nucl, [state, state], imsrg_params, kshell_params)
      scripts = imsrg_submit.gen_scripts_combine_delta(
          f"{imsrg_submit.output_dir}/{imsrg_submit.filebase}_R2p.csv",
          ops_rankJ=[0]
      )
      for stage_name, stage in stages.items():
        stage.add_task(SampleID, scripts[stage_name])

    chain.submit(max_concurrent=MAX_CONCURRENT, dry_run=DRY_RUN)
