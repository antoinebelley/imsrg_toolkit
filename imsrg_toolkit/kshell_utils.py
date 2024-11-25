import sys, os
from pathlib import Path
import numpy as np
import re
from imsrg_toolkit.periodictable import periodic_table
from textwrap import dedent
from subprocess import run, PIPE




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




class kshell_script():

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




class kshell_wavefunction_script(kshell_script):

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




class kshell_density_script(kshell_script):
  
  def __init__(self, fn_snt, Nucl_daughter=None, **kwargs):
    super().__init__(fn_snt)
    self.state_list = ["+1", "+1"]
    if not Nucl_daughter:
      self.Nucl_daughter = self.Nucl
    else:
      self.Nucl_daughter = Nucl_daughter
    self.Z_daughter, self.N_daughter, self.A_daughter = _ZNA_from_str(self.Nucl_daughter)
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
    fn_script = f"{self.scratch_directory}/density_{self.filebase}_{self.Nucl_daughter}{state_string(self.state_list[1], self.A_daughter)}_{self.Nucl}{state_string(self.state_list[0], self.A)}.sh"
    f = open(fn_script, "w")
    f.write(s)
    f.close()
    os.chmod(fn_script, 0o755)
    return fn_script




class kshell_tookit():

   def __init__(self, fn_snt, Nucl, state_list, Nucl_daughter=None, **kwargs):
    self.Nucl = Nucl
    self.submit_cmd = 'sbatch'
    self.state_list = state_list
    self.kshell_ket = kshell_wavefunction_script(fn_snt, Nucl = Nucl, **kwargs)
    if Nucl_daughter != None:
      self.Nucl_daughter == Nucl
      self.kshell_bra = kshell_wavefunction_script(fn_snt, Nucl = Nucl_daughter, **kwargs)
    else:
      self.kshell_bra == self.kshell_ket
    

  def gen_partion(self, parity, ket=True):
    if ket:
      self.kshell_ket.gen_partition(parity)
    else:
      self.kshell_bra.gen_parititon(parity)


  def calc_diag(self, gen_partition=True, submit = True):
    if gen_parition:
      if str_state(state_list[-1], self.Nucl)[-1] == 'p': 
        parity = 1
      else:
        parity = -1
      self.gen_partition(parity)
    self.ket_diag = self.kshell_ket.gen_script()
    if submit:
        jobid_ket = run([self.submission_cmd, '--parsable',self.ket_diag], stdout=PIPE, text=True).stdout
    if self.Nucl != self.Nucl_daughter:
      if gen_parition:
        if str_state(state_list[-1], self.Nucl_daughter)[-1] == 'p': 
          parity = 1
        else:
          parity = -1
        self.gen_partition(parity, ket=False)
      self.bra_diag = self.kshell_bra.gen_script()
      if submit:
        jobid_bra = run([self.submission_cmd, '--parsable',self.bra_diag], stdout=PIPE, text=True).stdout
      return jobid_bra, jobid_ket
    else:
      return jobid_ket
    

  def calc_density():
    pass


  def calc_opexpvals():
    if self.kshell_bra.A < self.kshell_ket.A or self.kshell_bra.Z < self.kshell_ket.Z:
      kshell_bra, kshell_ket = self.kshell_ket, self.kshell_bra