from imsrg_toolkit.kshell_utils import KshellToolkit
from pathlib import Path
import os
from subprocess import run, PIPE


class ImsrgParams():
  def __init__(self, sampleID = None, **kwargs):
    #### Here are the default parameters for the imsrg###
    ### TODO add all IMSRG parameters in the params
    #Paths to different directories that are used
    #TODO update those from a config file
    self.scratch_directory = '/work/submit/abelley/work/imsrg/'
    self.output_directory_base = '/home/submit/abelley/results/'
    self.file2b_directory = '/ceph/submit/data/group/ab-initio/me2j/'
    self.file3b_directory  = '/ceph/submit/data/group/ab-initio/me3j/'

    # Model space parameters
    self.Z = 2
    self.A = 6
    self.emax = 2
    self.E3max = 6
    self.hw = 10
    self.ref = 'He6'
    self.valence_space = 'p-shell'
    self.custom_valence_space = None

    #2B interaction parameters
    self.label = 'SampleDelta'
    self.file2e1max = 14
    self.file2e2max = 28
    self.file2lmax = 14

    #3B interaction parameters
    self.file3e1max = 16
    self.file3e2max = 32
    self.file3e3max = 28
    self.file3_format = 'no2b'
    self.file3_precision = 'half'

    #IMSRG solver parameters
    self.method = 'magnus'
    self.denominator_partitioning = 'Epstein_Nesbet'
    self.eta_criterion = 1e-6
    self.smax = 500
    self.dsmax = 0.5
    self.ds0 = 0.5
    self.denominator_delta = 0
    self.domega = 0.2
    self.omega_norm_max = 0.25
    self.ode_tolerance = 1e-6
    self.core_generator = 'atan'
    self.valence_space_generator = 'shell-model-atan'

    #Operators parameters
    self.opfiles = []
    self.opnames = []
    self.write_HO_ops = True
    self.write_HF_ops = True

    #If dictionay is given, update the attributes using the
    #dictionary keys and values.
    self.update_params(**kwargs)

    #generate the file name for the snt file
    self.output_dir = f"{self.output_directory_base}/{self.ref}/{self.label}/"
    self.gen_filebase(sampleID)
    self.intfile = f"{self.output_dir}/{self.filebase}"

  def update_params(self, **kwargs):
    for key, value in kwargs.items():
      setattr(self, key, value)


  def gen_filebase(self, sampleID = None):
    if not sampleID:
      self.filebase = f"{self.valence_space}_{self.label}_{self.ref}_{self.method}_e{self.emax}_E{self.E3max}_hw{self.hw}"
    else:
      self.filebase = f"{self.valence_space}_{self.label}_{sampleID}_{self.ref}_{self.method}_e{self.emax}_E{self.E3max}_hw{self.hw}"


class Utils():
  def __init__(self, Nucl, state_list, imsrg_params, kshell_params,  Nucl_daughter=None, submit_cmd='sbatch', module_path = "/work/submit/abelley/imsrg_toolkit/"):
    self.imsrg_params = imsrg_params
    self.module_path = module_path
    pars = ImsrgParams(**imsrg_params)
    fn_snt = pars.intfile
    self.opnames = pars.opnames
    self.output_dir = pars.output_dir
    fn_snt_path = Path(fn_snt)
    self.submit_cmd = submit_cmd
    self.fn_snt = fn_snt+".snt"
    self.kshell = KshellToolkit(self.fn_snt, Nucl, state_list, Nucl_daughter=Nucl_daughter, submit_cmd=submit_cmd,  **kshell_params)
    self.filebase = fn_snt_path.name
    self.fn_py = self.filebase+".py"
    self.fn_sh = self.filebase+".sh"
    self.scratch_directory = '/work/submit/abelley/work/imsrg/'
  
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
    script += f"imsrg.run('{file2b}', '{file3b}')\n"
    script = self.add_kshell_partition(script)
    f = open(fn_pyimsrg, "w")
    f.write(script)
    f.close()
    os.chmod(fn_pyimsrg, 0o755)
    return fn_pyimsrg


  def gen_imsrg_python_script_combine_delta(self, LECs, sampleID):
    fn_pyimsrg = self.scratch_directory+self.fn_py
    script = self.write_script_header()
    script += f"imsrg.run_combine_delta({LECs}, {sampleID})\n"
    script = self.add_kshell_partition(script)
    f = open(fn_pyimsrg, "w")
    f.write(script)
    f.close()
    os.chmod(fn_pyimsrg, 0o755)
    return fn_pyimsrg

  def add_kshell_partition (self, script):
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


  def gen_imsrg_submit_script_combine_delta(self, LECs, sampleID):
    fn_script = self.scratch_directory+self.fn_sh
    python_script = self.gen_imsrg_python_script_combine_delta(LECs, sampleID)
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
  

  def submit_imsrg_combine_delta(self, LECs, sampleID, verbose = True):
    fn_sh = self.gen_imsrg_submit_script_combine_delta(LECs, sampleID)
    jobid = run([self.submit_cmd, '--parsable', fn_sh], stdout=PIPE, text=True, check=True).stdout.rstrip()
    return self.submit_job(fn_sh, verbose=verbose)


  def submit_all(self, file2b, file3b, fn_output, ops_rankJ=None, ops_rankP=None, ops_rankZ=None,  header_expvals=None, verbose=False):
    imsrg_id = self.submit_imsrg(file2b, file3b, verbose=verbose)
    fn_ops = [f"{self.output_dir}{self.filebase}_{op}.snt" for op in self.opnames]
    self.kshell.submit_all(fn_output, fn_ops, previous_jobid = imsrg_id, ops_rankJ = ops_rankJ, ops_rankP = ops_rankP, ops_rankZ = ops_rankZ, header = header_expvals, verbose=verbose)

  def submit_all_combine_delta(self, LECs, sampleID, fn_output, ops_rankJ=None, ops_rankP=None, ops_rankZ=None,  header_expvals=None, verbose=False):
    imsrg_id = self.submit_imsrg_combine_delta(LECs, sampleID, verbose=verbose)
    fn_ops = [f"{self.output_dir}{self.filebase}_{op}.snt" for op in self.opnames]
    self.kshell.submit_all(fn_output, fn_ops, previous_jobid = imsrg_id, ops_rankJ = ops_rankJ, ops_rankP = ops_rankP, ops_rankZ = ops_rankZ, header = header_expvals, verbose=verbose)