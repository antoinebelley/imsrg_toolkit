from imsrg_toolkit.kshell_utils import KshellToolkit
from pathlib import Path
import os
from subprocess import run, PIPE
from imsrg_toolkit.settings import *
from imsrg_toolkit.imsrg_params import *


# @dataclass(frozen=False)
# class ImsrgParams():
#   #Paths for inputs and outputs
#   scratch_directory: str = f'/work/submit/{username}/work/imsrg/'
#   output_directory_base: str = f'/home/submit/{username}/results/'
#   file2b_directory: str = INTERACTION_2B_PATH
#   file3b_directory: str = INTERACTION_3B_PATH

#   #Model space parameters
#   A: int = 6
#   emax: int = 2
#   E3max: int = 6
#   hw: int = 10
#   ref: str = 'He6'
#   valence_space: str = 'p-shell'
#   custom_valence_space: str = None 

#    #2B interaction parameters 
#   label: str = 'SampleDelta'
#   SampleID: str = None #This is for the interactions samples only
#   file2e1max: int = 14
#   file2e2max: int = 28
#   file2lmax: int = 14    

#   #3B interaction parameters
#   file3e1max: int = 16
#   file3e2max: int = 32
#   file3e3max: int = 28
#   file3_format: str = 'no2b'
#   file3_precision: str = 'half'

#   #Parameters for the BetaCM
#   BetaCM: float = 0.0
#   hwBetaCM: float = hw  # Negative value means use the frequency

#   #IMSRG solver parameters
#   basis : str = 'HF'
#   method: str = 'magnus'
#   denominator_partitioning: str = 'Epstein_Nesbet'
#   eta_criterion: float = 1e-6
#   smax: int = 500
#   dsmax: float = 0.5
#   ds0: float = 0.5
#   denominator_delta: float = 0
#   denominator_delta_orbit: str = None
#   domega: float = 0.2
#   omega_norm_max: float = 0.25
#   ode_tolerance: float = 1e-6
#   core_generator: str = 'atan'
#   valence_space_generator: str = 'shell-model-atan'

#   #Operators parameters
#   opfiles: list = []
#   opnames: list = []
#   write_HO_ops: bool = True
#   write_HF_ops: bool = True


#   def gen_filebase(self):
#     if not self.sampleID:
#       self.filebase = f"{self.valence_space}_{self.label}_{self.ref}_{self.method}_e{self.emax}_E{self.E3max}_hw{self.hw}"
#     else:
#       self.filebase = f"{self.valence_space}_{self.label}_{self.sampleID}_{self.ref}_{self.method}_e{self.emax}_E{self.E3max}_hw{self.hw}"

#   def __post_init__(self):
#     #generate the file name for the snt file
#     self.output_dir = f"{self.output_directory_base}/{self.ref}/{self.label}/"
#     self.gen_filebase(self.sampleID)
#     self.intfile = f"{self.output_dir}/{self.filebase}"


class Utils():
  def __init__(self, Nucl, state_list, imsrg_params, kshell_params, Nucl_daughter=None, submit_cmd='sbatch', HF=False):
    self.imsrg_params = imsrg_params
    self.HF = HF
    self.module_path = ROOT_DIR
    pars = ImsrgParams(**imsrg_params)
    fn_snt = pars.intfile
    self.opnames = pars.opnames
    self.opfiles = pars.opfiles
    self.output_dir = pars.output_dir
    fn_snt_path = Path(fn_snt)
    self.submit_cmd = submit_cmd
    if HF:
      fn_snt = fn_snt+"_HF"
    self.fn_snt = fn_snt+".snt"
    self.kshell = KshellToolkit(self.fn_snt, Nucl, state_list, Nucl_daughter=Nucl_daughter, submit_cmd=submit_cmd,  **kshell_params)
    self.filebase = fn_snt_path.name
    self.fn_py = self.filebase+".py"
    self.fn_sh = self.filebase+".sh"
    self.scratch_directory = f'/work/submit/{username}/work/imsrg/'


  def write_script_header(self):
    script = "#!/usr/bin/env python3\n"
    script += "import sys\n"
    script += f"sys.path.append('{self.module_path}')\n"
    script += "from imsrg_toolkit.imsrg import Imsrg\n"
    script += "from imsrg_toolkit.kshell_utils import KshellToolkit\n"
    script += "params = {\n"
    for key, value in self.imsrg_params.items():
      if key == 'header' or key == 'run_cmd' : continue
      if type(value) == str:
        script += f"\t '{key}': '{value}',\n"
      else:
        script += f"\t '{key}': {value},\n"
    script+='}\n'
    script += "params_kshell = {\n"
    for key, value in self.kshell.params.items():
      if key == 'header' or key == 'run_cmd': continue
      if type(value) == str:
        script += f"\t '{key}': '{value}',\n"
      else:
        script += f"\t '{key}': {value},\n"
    script+='}\n'

    script += "imsrg = Imsrg(**params)\n"

    return script


  def gen_imsrg_python_script(self, file2b, file3b):
    fn_pyimsrg = self.scratch_directory+self.fn_py
    script = self.write_script_header()
    script += f"imsrg.run('{file2b}', '{file3b}', HF = {self.HF})\n"
    script = self.add_kshell_partition(script)
    f = open(fn_pyimsrg, "w")
    f.write(script)
    f.close()
    os.chmod(fn_pyimsrg, 0o755)
    return fn_pyimsrg


  def gen_imsrg_python_script_combine_delta(self):
    fn_pyimsrg = self.scratch_directory+self.fn_py
    script = self.write_script_header()
    script += f"imsrg.run_combine_delta(HF = {self.HF})\n"
    script = self.add_kshell_partition(script)
    f = open(fn_pyimsrg, "w")
    f.write(script)
    f.close()
    os.chmod(fn_pyimsrg, 0o755)
    return fn_pyimsrg


  def add_kshell_partition (self, script):
    if self.HF:
      script += (f"kshell = KshellToolkit(imsrg.intfile+'_HF.snt', "
      f"Nucl = '{self.kshell.Nucl}', "
      f"state_list={self.kshell.state_list}, ")
      if self.kshell.Nucl_daughter != self.kshell.Nucl:
        script += f"Nucl_daughter='{self.kshell.Nucl_daughter}', "
      script += (f"submit_cmd='{self.kshell.submit_cmd}', "
      f"**{self.kshell.params})\n"
      )
    else:
      script += (f"kshell = KshellToolkit(imsrg.intfile+'.snt', "
      f"Nucl = '{self.kshell.Nucl}', "
      f"state_list={self.kshell.state_list}, ")
      if self.kshell.Nucl_daughter != self.kshell.Nucl:
        script += f"Nucl_daughter='{self.kshell.Nucl_daughter}', "
      script += (f"submit_cmd='{self.kshell.submit_cmd}', "
      f"**{self.kshell.params})\n"
      )
    script += "kshell.gen_partition()\n"
    if self.kshell.Nucl != self.kshell.Nucl_daughter:
      script += "kshell.gen_partition(ket=False)\n"
    return script


  def write_submission_script(self, fn_script, file):
    submit_script = f"{self.imsrg_params['header']}\n"
    submit_script += f"{self.imsrg_params['run_cmd']} python3 {file}\n"
    f = open(fn_script, "w")
    f.write(submit_script)
    f.close()
    os.chmod(fn_script, 0o755)


  def gen_imsrg_submit_script(self, file2b, file3b):
    fn_script = self.scratch_directory+self.fn_sh
    python_script = self.gen_imsrg_python_script(file2b, file3b)
    self.write_submission_script(fn_script, python_script)
    return fn_script


  def gen_imsrg_submit_script_combine_delta(self):
    fn_script = self.scratch_directory+self.fn_sh
    python_script = self.gen_imsrg_python_script_combine_delta()
    self.write_submission_script(fn_script, python_script)
    return fn_script


  def submit_job(self, file, verbose=False):
    jobid = run([self.submit_cmd, '--parsable', file], stdout=PIPE, text=True, check=True).stdout.rstrip()
    if verbose:
      print(f'Submitted imsrg with jobid {jobid}')
    return jobid


  def submit_imsrg(self, file2b, file3b, verbose=False):
    fn_sh = self.gen_imsrg_submit_script(file2b, file3b)
    return self.submit_job(fn_sh, verbose=verbose)


  def submit_imsrg_combine_delta(self, verbose = True):
    fn_sh = self.gen_imsrg_submit_script_combine_delta()
    return self.submit_job(fn_sh, verbose=verbose)


  def gen_oplist(self):
    if self.HF:
      fn_ops = [f"{self.output_dir}{self.filebase}_{op}_HF.snt" for op in self.opnames]
    else:
      fn_ops = [f"{self.output_dir}{self.filebase}_{op}.snt" for op in self.opnames]
    if len(self.opfiles) > 0 :
      if self.HF:
        tmp = [f"{self.output_dir}{self.filebase}_{op[1]}_HF.snt" for op in self.opfiles]
      else:
        tmp = [f"{self.output_dir}{self.filebase}_{op[1]}.snt" for op in self.opfiles]
      fn_ops.extend(tmp)
    return fn_ops


  def submit_all(self, file2b, file3b, fn_output, ops_rankJ=None, ops_rankP=None, ops_rankZ=None,  header_expvals=None, verbose=False):
    imsrg_id = self.submit_imsrg(file2b, file3b, verbose=verbose)
    fn_ops = self.gen_oplist()
    self.kshell.submit_all(fn_output, fn_ops, previous_jobid = imsrg_id, ops_rankJ = ops_rankJ, ops_rankP = ops_rankP, ops_rankZ = ops_rankZ, header = header_expvals, verbose=verbose)


  def submit_all_combine_delta(self, fn_output, ops_rankJ=None, ops_rankP=None, ops_rankZ=None,  header_expvals=None, verbose=False):
    imsrg_id = self.submit_imsrg_combine_delta(verbose=verbose)
    fn_ops = self.gen_oplist()
    self.kshell.submit_all(fn_output, fn_ops, previous_jobid = imsrg_id, ops_rankJ = ops_rankJ, ops_rankP = ops_rankP, ops_rankZ = ops_rankZ, header = header_expvals, verbose=verbose)