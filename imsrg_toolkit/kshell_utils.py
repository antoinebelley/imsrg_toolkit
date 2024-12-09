import sys, os
from pathlib import Path
import numpy as np
import re
from imsrg_toolkit.periodictable import periodic_table
from textwrap import dedent
from subprocess import run, PIPE
from imsrg_toolkit.TransitionDensity import TransitionDensity
from imsrg_toolkit.Operator import Operator
import itertools
import pandas as pd




def _ZNA_from_str(Nucl):
    """
    ex.) Nucl="O16" -> Z=8, N=8, A=16
    """
    isdigit = re.search(r'\d+', Nucl)
    A = int( isdigit.group() )
    asc = Nucl[:isdigit.start()] + Nucl[isdigit.end():]
    asc = asc.lower()
    asc = asc[0].upper() + asc[1:]
    Z = periodic_table.index(asc)
    N = A-Z
    return Z, N, A




def state_string(state, A):
  """
  example:
  0+1 -> j0p
  0+2 -> j0p
  2+2 -> j4p
  0.5-2 -> j1n
  1.5-2 -> j3n
  +1 -> m0p or m1p
  -1 -> m0n or m1n
  """
  if(state.find("+")!=-1): _str = state.split("+")
  if(state.find("-")!=-1): _str = state.split("-")
  if(len(_str)!=2): raise ValueError("Input value is not correct: "+state)
  try:
    n = int(_str[1])
  except:
    raise ValueError("Input value is not correct: "+state)
  if(_str[0] == ""):
    if( state.find("+")!=-1): state_str = "m0p"
    if( state.find("-")!=-1): state_str = "m0n"
    if(A%2==1):
        if( state.find("+")!=-1): state_str = "m1p"
        if( state.find("-")!=-1): state_str = "m1n"
  elif(_str != ""):
    try:
        j_double = int(2*float(_str[0]))
    except:
        raise ValueError("Input value is not correct: "+state)
    if( state.find("+")!=-1): state_str = "j{:d}p".format(j_double)
    if( state.find("-")!=-1): state_str = "j{:d}n".format(j_double)
  return state_str




def _str_J_to_Jfloat(string):
    """
    '0' -> 0
    '1' -> 1
    '1/2' -> 0.5
    '3/2' -> 1.5
    """
    if(string.find("/")!=-1): return float(string[:-2])*0.5
    return float(string)




class KshellScript():
  def __init__(self, fn_snt):
    self.Nucl = "He6"
    self.Z, self.N, self.A = _ZNA_from_str(self.Nucl)
    self.path_to_kshell =  " "
    self.header = " "
    self.output_directory =  "/home/submit/abelley/results/kshell/"
    self.scratch_directory = "/work/submit/abelley/work/"
    fn_snt_path = Path(fn_snt)
    self.filebase = fn_snt_path.name[:-4]
    self.fn_base = self.Nucl+ "_" + self.filebase
    self.run_cmd = 'srun'
    self.fn_snt = fn_snt
    

  def update_params(self, **kwargs):
    for key, value in kwargs.items():
      setattr(self, key, value)
    self.Z, self.N, self.A = _ZNA_from_str(self.Nucl)




class KshellWavefunctionScript(KshellScript):
  def __init__(self, fn_snt, **kwargs):
    super().__init__(fn_snt)
    self.states = "+1" 
    self.hw_truncation = None
    self.ph_truncation = None
    self.snt_prm = {}
    self.read_snt()
    self.nf = self.element2nf()
    self.update_params(**kwargs)
  

  def read_comment_skip(self, fp):
    while 1:
        arr = fp.readline().split()
        if not arr: return None
        if arr[0] == '!namelist':
            if arr[2] != '=': raise 'ERROR namlist line in snt'
            var_dict[arr[1]]  = ' '.join(arr[3:])
        for i in range(len(arr)): 
            if arr[i][0]=="!" or arr[i][0]=="#": 
                arr = arr[:i]
                break
        if not arr: continue
        try:
            return [int(i) for i in arr]
        except ValueError:
            try:
                return [int(i) for i in arr[:-1]]+[float(arr[-1])]
            except ValueError:
                return arr


  def read_snt(self):
    fp = open(self.fn_snt, 'r')
    np, nn, ncp, ncn  = self.read_comment_skip(fp)
    norb, lorb, jorb, torb = [], [], [], []
    npn = [np, nn]
    nfmax = [0, 0]
    for i in range(np+nn):
        arr = self.read_comment_skip(fp)
        if i+1 != int(arr[0]): 
            print( 'snt index error', i, arr[0] )
            raise 
        norb.append( int(arr[1]) )
        lorb.append( int(arr[2]) )
        jorb.append( int(arr[3]) )
        torb.append( int(arr[4]) )
        nfmax[(int(arr[4])+1)//2] += int(arr[3]) + 1
    fp.close()
    self.snt_prm['ncore'] = (ncp, ncn)
    self.snt_prm['n_jorb'] = (np, nn)
    self.snt_prm['norb'] = norb
    self.snt_prm['lorb'] = lorb
    self.snt_prm['jorb'] = jorb
    self.snt_prm['torb'] = torb
    self.snt_prm['nfmax'] = nfmax


  def element2nf(self):
    digits = []
    letters = []
    for char in self.Nucl:
      if char.isdigit():
        digits.append(char)
      else:
        letters.append(char)
    ele = ''.join(digits + letters)
    isdigit = re.search(r'\d+', ele)
    if not isdigit:
        print( '\n *** Invalid: unknown element ***', ele )
        return False
    mass = int( isdigit.group() )
    asc = ele[:isdigit.start()] + ele[isdigit.end():]
    asc = asc.lower()
    asc = asc[0].upper() + asc[1:]
    if not asc in periodic_table: 
        print( '*** Invalid: unknown element ***', asc )
        return False
    z = periodic_table.index(asc)
    corep, coren = self.snt_prm['ncore']
    
    if corep >= 0: nf1 =  z - corep
    else:          nf1 = -z - corep
    if coren >= 0: nf2 =   mass - z  - coren
    else:          nf2 = -(mass - z) - coren
        
    # print( '\n number of active particles ', nf1, nf2 )
    
    if nf1 < 0 or nf2 < 0 or \
       nf1 > self.snt_prm['nfmax'][0] or \
       nf2 > self.snt_prm['nfmax'][1]:
        print( '*** ERROR: nuclide out of model space ***' )
        return False
    return (nf1, nf2)


  def get_occupation(self, hw_ex=False):
    H = Operator()
    H.read_operator_file(self.fn_snt,A=self.A)
    logs = []
    states = self.states.split(",")
    for state in states:
      str_state = state_string(self.states, self.A)
      log = f"{self.scratch_directory}/log_{self.fn_base}_{str_state}.txt"
      logs.append(log)
    e_data = {}
    Njpi = {}
    for log in logs:
        if(not os.path.exists(log)):
            print(f"{log} is not found")
            continue
        f = open(log,"r")
        while True:
            line = f.readline()
            if(not line): break
            dat = line.split()
            if(len(dat) < 2): continue
            if(dat[1] == "<H>:"):
                dat = line.split()
                n_eig= int(dat[0])
                ene  = float(dat[2]) + H.get_0bme()
                mtot = int(dat[6][:-2])
                J = dat[6]
                if(self.A%2==0): J = str(int(dat[6][:-2])//2)
                prty = int(dat[8])
                if(not (J,prty) in Njpi): Njpi[(J,prty)]=1
                else: Njpi[(J,prty)]+=1
                hws = None
                while ene in e_data: ene += 0.000001
                line = f.readline()
                dat = line.split()
                if(dat[0]=="<Hcm>:"): tt = int(dat[5][:-2])
                if(dat[0]=="<TT>:"): tt = int(dat[3][:-2])
                line = f.readline()
                data = line.split()
                if(line[0:7] ==" <p Nj>"):
                    plist = []
                    for i in range(len(data)-2):
                        plist.append(float(data[i+2]))
                line = f.readline()
                data = line.split()
                if(line[0:7] ==" <n Nj>"):
                    nlist = []
                    for i in range(len(data)-2):
                        nlist.append(float(data[i+2]))
                if(hw_ex):
                    while len(line)!=0:
                        line = f.readline()
                        data = line.split()
                        if(line[0:4] ==" hw:"):
                            hws = {}
                            for i in range(len(data)-1):
                                hw, prob = data[i+1].split(":")
                                hws[int(hw)] = float(prob)
                            break
                if(hws!=None): e_data[ (J,prty,Njpi[(J,prty)]) ] = (ene, log, tt, plist, nlist, hws)
                if(hws==None): e_data[ (J,prty,Njpi[(J,prty)]) ] = (ene, log, tt, plist, nlist)
        f.close()
    return e_data


  def get_wf_index(self):
    jpn_to_idx = {}
    logs = set()
    idxs = {}
    data = self.get_occupation()
    for key in data.keys():
        fn_log = data[key][1]
        if( fn_log in logs ):
            idxs[ fn_log ] += 1
        else:
            idxs[ fn_log ] = 1
            logs.add( fn_log )
        jpn_to_idx[key] = (fn_log, idxs[ fn_log ])
    return jpn_to_idx


  def gen_partition(self, parity):
    #parity : "1" or "-1" 
    from imsrg_toolkit import gen_partition
    fn_ptn = self.scratch_directory + self.fn_base
    if parity == 1:
      fn_ptn += "_p"
    elif parity == -1 :
      fn_ptn += "_n"
    fn_ptn += ".ptn"
    self.fn_ptn = fn_ptn
    if self.hw_truncation == None and self.ph_truncation == None:
      tmod = 0
      truncation_params = None
    elif self.hw_truncation == None and self.ph_truncation != None:
      #TODO this actually need to be implemented in gen_partition.py
      tmod = 1
      truncation_params = self.ph_truncation
    if self.hw_truncation != None and seld.ph_truncation == None:
      tmod = 2
      truncation_params = self.hw_truncation
    #TODO add other options for the truncations of the model space
    gen_partition.main(self.fn_snt, fn_ptn, self.nf, parity, tmod, truncation_params)


  def gen_script(self, gen_partition = False):
    if self.states[0] == "+" or self.states[0]=="-":
      m = 0
      if self.A%2 ==1:
        m = 1
    else:
      J = float(findall(r"[-+]?\d*\.*\d+", self.states)[0])
      m = int(2*J)
    str_state = state_string(self.states, self.A)
    if gen_partition and str_state[-1] == 'p':
      self.gen_partition(1)
    elif gen_partition and str_state[-1] == 'n':
      self.gen_partition(-1)
    jdouble = 'true'
    if str_state[0] ==  "m":
      jdouble = 'false'
    s=f"{self.header}\n"
    s += dedent(f"""
      # ---------- {self.fn_base} --------------
      cat > {self.fn_base}_{str_state}.input <<EOF
      &input
      beta_cm = 0
      eff_charge = 1.5, 0.5, 
      fn_int = "{self.fn_snt}"
      fn_ptn = "{self.fn_ptn}"
      fn_save_wave = "{self.fn_base}_{str_state}.wav"
      gl = 1.0, 0.0, 
      gs = 3.91, -2.678, 
      hw_type = 1
      is_double_j = .{jdouble}.
      max_lanc_vec = 200
      maxiter = 300
      mode_lv_hdd = 0
      mtot = {m}
      n_eigen = 1
      n_restart_vec = 10
      &end
      EOF
    """)
    s+= f"{self.run_cmd} ./kshell.exe {self.fn_base}_{str_state}.input > log_{self.fn_base}_{str_state}.txt 2>&1\n" 
    s+= dedent(f"""
      rm -f tmp_snapshot_{Path(self.fn_ptn).name}_0_* tmp_lv_{Path(self.fn_ptn).name}_0_* {self.fn_base}_{str_state}.input 

      ./collect_logs.py log_*{self.fn_base}* > summary_{self.fn_base}.txt
      cp summary_{self.fn_base}.txt {self.output_directory}
    """)
    s = dedent(s)
    fn_script = f"{self.scratch_directory}/{self.fn_base}.sh"
    f = open(fn_script, "w")
    f.write(s)
    f.close()
    os.chmod(fn_script, 0o755)
    return fn_script


  def summary_to_dictionary(self, comment_snt="!"):
    fn_summary = f'{self.output_directory}/summary_{self.fn_base}.txt'
    H = Operator()
    H.read_operator_file(self.fn_snt,comment=comment_snt,A=self.A)
    if(not os.path.exists(fn_summary)): return {}
    f = open(fn_summary,'r')
    lines = f.readlines()
    f.close()
    edict={}
    for line in lines:
        data = line.split()
        try:
            N = int(data[0])
            J = data[1]
            P = data[2]
            i = int(data[3])
            e = float(data[5])
            eex = float(data[6])
            edict[(J,P,i)] = e + H.get_0bme()
        except:
            continue
    return edict




class KshellDensityScript(KshellScript):
  def __init__(self, fn_snt, Nucl_daughter=None, **kwargs):
    super().__init__(fn_snt)
    self.state_list = ["+1", "+1"]
    if Nucl_daughter == None:
      self.Nucl_daughter = self.Nucl
    else:
      self.Nucl_daughter = Nucl_daughter
    self.Z_daughter, self.N_daughter, self.A_daughter = _ZNA_from_str(self.Nucl_daughter)
    self.fn_script = f"{self.scratch_directory}/density_{self.filebase}_{self.Nucl_daughter}{state_string(self.state_list[1], self.A_daughter)}_{self.Nucl}{state_string(self.state_list[0], self.A)}.sh"
    self.fn_density = f"{self.scratch_directory}/density_{self.filebase}_{self.Nucl_daughter}{state_string(self.state_list[1], self.A_daughter)}_{self.Nucl}{state_string(self.state_list[0], self.A)}.txt"
    self.update_params(**kwargs)
    

  def gen_script(self, fn_ptn, fn_ptn_daughter=None):
    if not fn_ptn_daughter:
      fn_ptn_daughter = fn_ptn
    s = f"{self.header}\n"
    s += dedent(f"""
      cat >density_{self.filebase}_{self.Nucl_daughter}{state_string(self.state_list[1], self.A_daughter)}_{self.Nucl}{state_string(self.state_list[0], self.A)}.input <<EOF
      &input
        fn_int   = "{self.fn_snt}"
        fn_ptn_l = "{fn_ptn_daughter}"
        fn_ptn_r = "{fn_ptn}"
        fn_load_wave_l = "{self.Nucl_daughter}_{self.filebase}_{state_string(self.state_list[1], self.A_daughter)}.wav"
        fn_load_wave_r = "{self.Nucl}_{self.filebase}_{state_string(self.state_list[0], self.A)}.wav"
        hw_type = 2
        eff_charge = 1.5, 0.5
        gl = 1.0, 0.0
        gs = 3.91, -2.678
        is_tbtd = .true.
      &end
      EOF
      """)
    s += dedent(f"{self.run_cmd} ./transit.exe density_{self.filebase}_{self.Nucl_daughter}{state_string(self.state_list[1], self.A_daughter)}_{self.Nucl}{state_string(self.state_list[0], self.A)}.input > density_{self.filebase}_{self.Nucl_daughter}{state_string(self.state_list[1], self.A_daughter)}_{self.Nucl}{state_string(self.state_list[0], self.A)}.txt 2>&1\n")
    s += f"rm density_{self.filebase}_{self.Nucl_daughter}{state_string(self.state_list[1], self.A_daughter)}_{self.Nucl}{state_string(self.state_list[0], self.A)}.input"
    f = open(self.fn_script, "w")
    f.write(s)
    f.close()
    os.chmod(self.fn_script, 0o755)
    return self.fn_script




class KshellToolkit():
  def __init__(self, fn_snt, Nucl, state_list, Nucl_daughter=None, submit_cmd='sbatch', **kwargs):
    self.Nucl = Nucl
    self.Z, self.N, self.A = _ZNA_from_str(self.Nucl)
    self.fn_snt = fn_snt
    self.submit_cmd = submit_cmd
    self.state_list = state_list
    self.params = kwargs
    self.kshell_ket = KshellWavefunctionScript(fn_snt, Nucl = Nucl, states=state_list[-1], **kwargs)
    if Nucl_daughter != None:
      self.Nucl_daughter = Nucl_daughter
      self.Z_daughter, self.N_daughter, self.A_daughter = _ZNA_from_str(self.Nucl_daughter)
      self.kshell_bra = KshellWavefunctionScript(fn_snt, Nucl = Nucl_daughter, states=state_list[0], **kwargs)
      self.density_script = KshellDensityScript(fn_snt, Nucl_daughter = Nucl_daughter, Nucl = Nucl, state_list = state_list, **kwargs)
    else:
      self.Nucl_daughter = Nucl
      self.kshell_bra = self.kshell_ket
      self.density_script = KshellDensityScript(fn_snt, Nucl = Nucl, state_list=state_list, **kwargs)
    self.outputs = []
    

  def gen_partition(self, parity, ket=True):
    if ket:
      self.kshell_ket.gen_partition(parity)
    else:
      self.kshell_bra.gen_parititon(parity)


  def submit_diag(self, gen_partition=True, previous_jobid = -1, verbose = False):
    if gen_partition:
      if state_string(self.state_list[-1], self.A)[-1] == 'p': 
        parity = 1
      else:
        parity = -1
      self.gen_partition(parity)
    ket_sh = self.kshell_ket.gen_script()
    jobid_ket = run([self.submit_cmd, '--parsable', f'--dependency=afterok:{previous_jobid}', ket_sh], stdout=PIPE, text=True, check=True).stdout.rstrip()
    if verbose:
      print(f'Submitted ket diagonalization with jobid {jobid_ket}')
    if self.Nucl != self.Nucl_daughter:
      if gen_partition:
        if state_string(self.state_list[0], self.A_daughter)[-1] == 'p': 
          parity = 1
        else:
          parity = -1
        self.gen_partition(parity, ket=False)
      bra_sh = self.kshell_bra.gen_script()
      jobid_bra = run([self.submit_cmd, '--parsable', f'--dependency=afterok:{previous_jobid}', bra_sh], stdout=PIPE, text=True, check=True).stdout.rstrip()
      if verbose: 
        print(f'Submitted bra diagonalization with jobid {jobid_bra}')
      return [jobid_bra, jobid_ket]
    else:
      return [jobid_ket]
    

  def submit_density(self, previous_jobids=-1, verbose = False):
    if self.Nucl_daughter != self.Nucl:
      fn_sh = self.density_script.gen_script(self.kshell_ket.fn_ptn)
    else:
      fn_sh = self.density_script.gen_script(self.kshell_ket.fn_ptn, self.kshell_bra.fn_ptn)
    if previous_jobids != -1:
      previous_jobids = ':'.join(map(str, previous_jobids))
    jobid_density = run([self.submit_cmd, '--parsable', f'--dependency=afterok:{previous_jobids}', fn_sh], stdout=PIPE, text=True).stdout.rstrip()
    # else:
    #   jobid_density = run([self.submit_cmd, '--parsable', fn_sh], stdout=PIPE, text=True).stdout
    if verbose:
      print(f'Submitted density with jobid {jobid_density}')
    return jobid_density
  

  def calc_opexpvals(self, fn_op, op_rankJ = 0, op_rankP = 1, op_rankZ = 0):
    op = Operator(filename=fn_op, rankJ=op_rankJ, rankP=op_rankP, rankZ=op_rankZ)
    if self.kshell_bra.A < self.kshell_ket.A or self.kshell_bra.Z < self.kshell_ket.Z:
      kshell_bra, kshell_ket = self.kshell_ket, self.kshell_bra
    else:
      kshell_bra = self.kshell_bra
      kshell_ket = self.kshell_ket
    wf_index_bra = kshell_bra.get_wf_index()
    wf_index_ket = kshell_ket.get_wf_index()
    energies_bra = kshell_bra.summary_to_dictionary()
    energies_ket = kshell_ket.summary_to_dictionary()
    exp_vals = pd.DataFrame()
    for state_bra, state_ket in itertools.product(list(wf_index_bra.keys()), list(wf_index_ket.keys())):
      Jbra, Pbra, nn_bra = state_bra
      Jket, Pket, nn_ket = state_ket
      Jfbra = _str_J_to_Jfloat(Jbra)
      Jfket = _str_J_to_Jfloat(Jket)
      if(Pbra * Pket * op_rankP == -1): continue
      if( not int(abs(Jfbra-Jfket)) <= op_rankJ <= int(Jfbra+Jfket) ): continue
      if Pbra == 1: Pbra = "+"
      elif Pbra == -1 : Pbra = "-"
      if Pket == 1: Pket = "+"
      elif Pket == -1 : Pket = "-"
      en_bra = energies_bra[(Jbra,Pbra,nn_bra)]
      en_ket = energies_ket[(Jket,Pket,nn_ket)]
      Density = TransitionDensity(filename=self.density_script.fn_density, Jbra=Jfket, wflabel_bra=wf_index_bra[state_bra][-1], \
                  Jket=Jfbra, wflabel_ket=wf_index_ket[state_ket][-1])
      output = Density.eval(op)
      output = [fn_op,self.Nucl_daughter,Jbra,Pbra,nn_bra,en_bra,self.Nucl,Jket,Pket,nn_ket,en_ket,*output]
      self.outputs.append(output)
    

  def gen_df_from_outputs(self):
    self.df = pd.DataFrame(self.outputs)
    self.df.columns = ["fn_op", "Nucl bra","J bra","P bra","n bra","Energy bra","Nucl ket","J ket","P ket","n ket","Energy ket","Zero","One","Two"]


  def gen_expvals_script(self, fn_ops, ops_rankJ=None, ops_rankP=None, ops_rankZ=None, header=None):
    if ops_rankJ == None:
      ops_rankJ = [0 for _ in  fn_ops]
    if ops_rankP == None:
      ops_rankP = [1 for _ in  fn_ops]
    if ops_rankZ == None:
      ops_rankZ = [0 for _ in  fn_ops]

    fn_eval = self.kshell_bra.scratch_directory+self.kshell_bra.filebase+"_eval.py"
    eval_script = "#!/usr/bin/env python3\n"
    if header != None:
      evals_script += f"{header}\n"
    eval_script += "from imsrg_toolkit.kshell_utils import kshell_toolkit\n"
    eval_script += "params = {\n"
    for key, value in self.params.items():
      if key == 'header' or key == 'run_cmd' : continue
      if type(value) == str:
        eval_script += f"\t '{key}': '{value}',\n"
      else:
        eval_script += f"\t '{key}': {value},\n"
    eval_script+='}\n'
    eval_script += f"vals = kshell_toolkit('{self.fn_snt}', '{self.Nucl}', {self.state_list}, **params)\n"
    for op, op_rankJ, op_rankP, op_rankZ  in zip(fn_ops, ops_rankJ, ops_rankP, ops_rankZ):
      eval_script += f"vals.calc_opexpvals('{op}', op_rankJ = {op_rankJ}, op_rankP = {op_rankP}, op_rankZ = {op_rankZ})\n"
    eval_script += "vals.gen_df_from_outputs()"
    f = open(fn_eval, "w")
    f.write(eval_script)
    f.close()
    os.chmod(fn_eval, 0o755)
    return fn_eval


  def submit_expvals(self, fn_ops, ops_rankJ=None, ops_rankP=None, ops_rankZ=None,  previous_jobid = -1, verbose = False, header=None):
    fn_sh = self.gen_expvals_script(fn_ops,  ops_rankJ=ops_rankJ, ops_rankP=ops_rankP, ops_rankZ=ops_rankZ, header = header)
    jobid = run([self.submit_cmd, '--parsable', f'--dependency=afterok:{previous_jobid}', fn_sh], stdout=PIPE, text=True, check=True).stdout.rstrip()
    if verbose:
      print(f'Submitted expvals with jobid {jobid}')
    return jobid


  def submit_all(self, fn_ops = [], ops_rankJ=None, ops_rankP=None, ops_rankZ=None,  previous_jobid = -1, verbose = False, header=None):
    os.chdir(self.kshell_ket.scratch_directory)
    #Submit the diagonalization
    diag_ids = self.submit_diag(previous_jobid=previous_jobid, verbose = verbose)
    #Submit the density
    density_id = self.submit_density(previous_jobids = diag_ids,  verbose = verbose)
    if len(fn_ops) > 0:
      #Submit the exp vals calculations
      self.submit_expvals(fn_ops, ops_rankJ = ops_rankJ, ops_rankP = ops_rankP, ops_rankZ = ops_rankZ, previous_jobid = density_id, verbose=verbose, header=header)
    