import sys
# from imsrg_toolkit.imsrg import Imsrg
from imsrg_toolkit.kshell_utils import KshellWavefunctionScript, KshellDensityScript, KshellToolkit
from imsrg_toolkit.utils import Utils
import numpy as np
import pandas as pd



LECs = ['Ct1S0pp','Ct1S0np','Ct1S0nn','Ct3S1','C1S0','C3P0','C1P1','C3P1','C3S1','CE1','C3P2','c1','c2','c3','c4','cD','cE']
df = pd.read_csv("/work/submit/abelley/imsrg_toolkit/data/8000Samples.txt")

index = np.array(df.index)
rng = np.random.default_rng()
index = rng.choice(index, 1, replace=False, shuffle=False)
Nucl = "O14"

imsrg_params = {}
imsrg_params['emax'] = 4
imsrg_params['E3max'] = 12
imsrg_params['hw'] = 10
imsrg_params['A'] = 14
imsrg_params['opnames'] = ['Rp2']
imsrg_params['ref'] = Nucl
imsrg_params['valence_space'] = 'p-shell' # this is just a label when custom_valence_space is set
imsrg_params['label'] = 'SampleDelta'
imsrg_params['run_cmd'] = """\
  srun apptainer exec \\
    --bind /home/submit \\
    --bind /work/submit \\
    --bind /scratch/submit \\
    --bind /cvmfs \\
    --bind /ceph/submit \\
    /work/submit/abelley/pyimsrg.sif """

kshell_params = {}
kshell_params['run_cmd'] = """\
  srun apptainer exec \\
    --bind /home/submit \\
    --bind /work/submit \\
    --bind /scratch/submit \\
    --bind /cvmfs \\
    --bind /ceph/submit \\
    /work/submit/abelley/work/kshell/kshell.sif """
kshell_params['header'] = """#!/bin/bash
#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=/work/submit/abelley/results/kshell_log/outputs/test_out_%j.txt
#SBATCH --error=/work/submit/abelley/results/kshell_log/errors/test_err_%j.txt
#SBATCH --time=10:00 """

for i in index:
  sample = df.iloc[i]
  print(sample.keys)
  SampleID = int(sample["SampleID"])
  weights = list(sample[LECs])

  imsrg_params['header'] = f"""#!/bin/bash
#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=/work/submit/abelley/results/imsrg_log/outputs/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_%j.txt
#SBATCH --error=/work/submit/abelley/results/imsrg_log/errors/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_%j.txt
#SBATCH --time=00:02:00
#SBATCH --mem=256M

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=24
"""





  imsrg_params['SampleID'] =  SampleID
  imsrg_params['LECs'] =  weights
  header_expvals = f"""#SBATCH --output=/work/submit/abelley/results/kshell_log/outputs/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_eval_%j.txt
# #SBATCH --error=/work/submit/abelley/results/kshell_log/errors/{imsrg_params['ref']}_emax{imsrg_params['emax']}_Sample{SampleID}_eval_%j.txt"""

  imsrg_submit = Utils(Nucl, ["+1", "+1"], imsrg_params, kshell_params)
  imsrg_submit.submit_all_combine_delta(f"{imsrg_submit.output_dir}/{imsrg_submit.filebase}_R2p.csv", header_expvals = header_expvals)
