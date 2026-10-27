from imsrg_toolkit.kshell_utils import KshellToolkit
from pathlib import Path
import os
from subprocess import run, PIPE
from imsrg_toolkit.settings import *
from imsrg_toolkit.imsrg_params import *


class Utils():
  def __init__(self, Nucl, state_list, imsrg_params, kshell_params, Nucl_daughter=None, submit_cmd='sbatch', HF=False):
    self.imsrg_params = imsrg_params
    self.HF = HF
    self.module_path = ROOT_DIR
    pars = ImsrgParams(**imsrg_params)
    fn_snt = pars.intfile
    self.opnames = pars.opnames
    self.opfiles = pars.opfiles
    self.opnames_decay = pars.opnames_decay
    self.output_dir = pars.output_dir
    fn_snt_path = Path(fn_snt)
    self.submit_cmd = submit_cmd
    if HF:
      fn_snt = fn_snt+"_HF"
    self.fn_snt = fn_snt+".snt"
    decay_mixed = False
    if len(pars.opnames_decay) > 0:
      decay_mixed = True
    self.kshell = KshellToolkit(self.fn_snt, Nucl, state_list, Nucl_daughter=Nucl_daughter, decay_mixed=decay_mixed, submit_cmd=submit_cmd,  **kshell_params)
    self.filebase = fn_snt_path.name
    self.fn_py = self.filebase+".py"
    self.fn_sh = self.filebase+".sh"
    self.imsrg_directory = pars.imsrg_scripts_directory


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
    skips  = ['header','header_daughter','header_density', 'run_cmd']
    script += "params_kshell = {\n"
    for key, value in self.kshell.params.items():
      if key in skips: continue
      if type(value) == str:
        script += f"\t '{key}': '{value}',\n"
      else:
        script += f"\t '{key}': {value},\n"
    script+='}\n'

    script += "imsrg = Imsrg(**params)\n"

    return script


  def write_script_header_anapole(self):
    script = "#!/usr/bin/env python3\n"
    script += "import sys\n"
    script += f"sys.path.append('{self.module_path}')\n"
    script += "from imsrg_toolkit.imsrg import PvImsrg\n"
    script += "from imsrg_toolkit.kshell_utils import KshellToolkit\n"
    script += "params = {\n"
    for key, value in self.imsrg_params.items():
      if key == 'header' or key == 'run_cmd' : continue
      if type(value) == str:
        script += f"\t '{key}': '{value}',\n"
      else:
        script += f"\t '{key}': {value},\n"
    script+='}\n'
    skips  = ['header','header_daughter','header_density', 'run_cmd']
    script += "params_kshell = {\n"
    for key, value in self.kshell.params.items():
      if key in skips: continue
      if type(value) == str:
        script += f"\t '{key}': '{value}',\n"
      else:
        script += f"\t '{key}': {value},\n"
    script+='}\n'

    script += "imsrg = PvImsrg(**params)\n"

    return script


  def gen_imsrg_python_script(self, file2b, file3b):
    fn_pyimsrg = self.imsrg_directory+self.fn_py
    script = self.write_script_header()
    script += f"imsrg.run('{file2b}', '{file3b}', HF = {self.HF})\n"
    script = self.add_kshell_partition(script)
    f = open(fn_pyimsrg, "w")
    f.write(script)
    f.close()
    os.chmod(fn_pyimsrg, 0o755)
    return fn_pyimsrg
  

  def gen_imsrg_python_script_anapole(self, file2b, file3b, scale=1000, staged=True):
    fn_pyimsrg = self.imsrg_directory+self.fn_py
    script = self.write_script_header_anapole()
    script += f"imsrg.run_anapole('{file2b}', '{file3b}', scale = {scale}, staged = {staged})\n"
    script = self.add_kshell_partition(script)
    f = open(fn_pyimsrg, "w")
    f.write(script)
    f.close()
    os.chmod(fn_pyimsrg, 0o755)
    return fn_pyimsrg
  

  def gen_imsrg_python_script_anapole_staged_split(self, file2b, file3b, VPNC_path=None):
    fn_pyimsrg = self.imsrg_directory+self.fn_py
    script = self.write_script_header_anapole()
    script += f"imsrg.run_anapole_staged_split('{file2b}', '{file3b}', VPNC_path='{VPNC_path}')\n"
    script = self.add_kshell_partition(script)
    f = open(fn_pyimsrg, "w")
    f.write(script)
    f.close()
    os.chmod(fn_pyimsrg, 0o755)
    return fn_pyimsrg


  def gen_imsrg_python_script_combine_delta(self):
    fn_pyimsrg = self.imsrg_directory+self.fn_py
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
    fn_script = self.imsrg_directory+self.fn_sh
    python_script = self.gen_imsrg_python_script(file2b, file3b)
    self.write_submission_script(fn_script, python_script)
    return fn_script

  def gen_imsrg_submit_script_anapole(self, file2b, file3b, staged=True, scale=1000):
    fn_script = self.imsrg_directory+self.fn_sh
    python_script = self.gen_imsrg_python_script_anapole(file2b, file3b, staged=staged, scale=scale)
    self.write_submission_script(fn_script, python_script)
    return fn_script
  
  def gen_imsrg_submit_script_anapole_staged_split(self, file2b, file3b, VPNC_path=None):
    fn_script = self.imsrg_directory+self.fn_sh
    python_script = self.gen_imsrg_python_script_anapole_staged_split(file2b, file3b, VPNC_path=VPNC_path)
    self.write_submission_script(fn_script, python_script)
    return fn_script


  def gen_imsrg_submit_script_combine_delta(self):
    fn_script = self.imsrg_directory+self.fn_sh
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


  def submit_anapole(self, file2b, file3b, verbose=False, staged=True, scale = 1000):
    fn_sh = self.gen_imsrg_submit_script_anapole(file2b, file3b, scale=scale, staged=staged)
    return self.submit_job(fn_sh, verbose=verbose)


  def submit_anapole_staged_split(self, file2b, file3b, verbose=False, VPNC_path=None):
    fn_sh = self.gen_imsrg_submit_script_anapole_staged_split(file2b, file3b, VPNC_path = VPNC_path)
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


  def gen_oplist_decay(self):
    if self.HF:
      fn_ops = [f"{self.output_dir}{self.filebase}_{op}_HF.snt" for op in self.opnames_decay]
    else:
      fn_ops = [f"{self.output_dir}{self.filebase}_{op}.snt" for op in self.opnames_decay]
    return fn_ops


  def submit_all(self, file2b, file3b, fn_output, ops_rankJ=None, ops_rankP=None, ops_rankZ=None,  header_expvals=None, verbose=False):
    imsrg_id = self.submit_imsrg(file2b, file3b, verbose=verbose)
    fn_ops = self.gen_oplist()
    self.kshell.submit_all(fn_output, fn_ops, previous_jobid = imsrg_id, ops_rankJ = ops_rankJ, ops_rankP = ops_rankP, ops_rankZ = ops_rankZ, header = header_expvals, verbose=verbose)


  def gen_scripts_combine_delta(self, fn_output, ops_rankJ=None, ops_rankP=None, ops_rankZ=None, header_expvals=None):
    """Generate all per-sample scripts (imsrg, diag, density, expvals) without
    submitting them, so they can be submitted as chained job arrays
    (see imsrg_toolkit.job_array)."""
    fn_sh = self.gen_imsrg_submit_script_combine_delta()
    fn_ops = self.gen_oplist()
    scripts = {'imsrg': [fn_sh]}
    scripts.update(self.kshell.gen_scripts(fn_output, fn_ops, ops_rankJ=ops_rankJ, ops_rankP=ops_rankP, ops_rankZ=ops_rankZ, header=header_expvals))
    return scripts


  def submit_all_combine_delta(self, fn_output, ops_rankJ=None, ops_rankP=None, ops_rankZ=None,
                               ops_rankJ_decay=None, ops_rankP_decay=None, ops_rankZ_decay=None, 
                               header_expvals=None, verbose=False):
    imsrg_id = self.submit_imsrg_combine_delta(verbose=verbose)
    fn_ops = self.gen_oplist()
    fn_ops_decay = self.gen_oplist_decay()
    self.kshell.submit_all(fn_output, fn_ops, fn_ops_decay = fn_ops_decay, previous_jobid = imsrg_id, 
                           ops_rankJ = ops_rankJ, ops_rankP = ops_rankP, ops_rankZ = ops_rankZ,
                           ops_rankJ_decay = ops_rankJ_decay, ops_rankP_decay = ops_rankP_decay, ops_rankZ_decay = ops_rankZ_decay, 
                           header = header_expvals, verbose=verbose)
  

  def submit_all_anapole(self, file2b, file3b, fn_output, ops_rankJ=None, ops_rankP=None, ops_rankZ=None,  header_expvals=None, verbose=False, staged=True, scale = 1000):
    imsrg_id = self.submit_anapole(file2b, file3b, verbose=verbose, staged = staged, scale=scale)
    fn_ops = fn_ops = [f"{self.output_dir}{self.filebase}_Anapolepp.snt"]
    self.kshell.submit_anapole(fn_output, fn_ops, previous_jobid = imsrg_id, ops_rankJ = ops_rankJ, ops_rankP = ops_rankP, ops_rankZ = ops_rankZ, header = header_expvals, verbose=verbose, scale=scale)


  def submit_all_anapole_staged_split(self, file2b, file3b, fn_output,   
                        header_expvals=None, verbose=False, VPNC_path=None):
    imsrg_id = self.submit_anapole_staged_split(file2b, file3b, 
                                                verbose=verbose, VPNC_path=VPNC_path)
    names = ['best', 'hpi1', 'hrho0', 'hrho1', 'hrho2', 'hrho1p', 'homega0', 'homega1']                                                
    fn_ops = [f"{self.output_dir}{self.filebase}_Anapolepp_{name}.snt" for name in names]
    ops_rankJ = [1 for i in range(len(names))]
    ops_rankZ = [0 for i in range(len(names))]
    ops_rankP = [1 for i in range(len(names))]
    self.kshell.submit_anapole(fn_output, fn_ops, 
                                previous_jobid = imsrg_id, ops_rankJ = ops_rankJ, 
                                ops_rankP = ops_rankP, ops_rankZ = ops_rankZ, 
                                header = header_expvals, verbose=verbose)