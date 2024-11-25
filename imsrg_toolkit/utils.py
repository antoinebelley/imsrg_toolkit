from pyIMSRG import *
from sys import stdout
import sys, os
from pathlib import Path
import numpy as np
import re
from imsrg_toolkit.periodictable import periodic_table


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




class imsrg_toolkit():

  def __init__(self, **kwargs):
    #### Here are the default parameters for the imsrg###
    ### TODO add all IMSRG parameters in the params
    #Paths to different directories that are used
    #TODO update those from a config file
    self.scratch_directory = '/work/submit/abelley/work/'
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


  def update_params(self, **kwargs):
    for key, value in kwargs.items():
      setattr(self, key, value)


  def set_imsrgsolver(self):
    """
      Intialize the imsrgsolver based on the parameters saved in params.
    """
    self.imsrgsolver.SetMethod(self.method)
    self.imsrgsolver.SetDenominatorPartitioning(self.denominator_partitioning)
    self.imsrgsolver.SetEtaCriterion(self.eta_criterion)
    self.imsrgsolver.SetSmax(self.smax)
    self.imsrgsolver.SetDsmax(self.dsmax)
    self.imsrgsolver.SetDs(self.ds0)
    self.imsrgsolver.SetDenominatorDelta(self.denominator_delta)
    self.imsrgsolver.SetdOmega(self.domega)
    self.imsrgsolver.SetOmegaNormMax(self.omega_norm_max)
    self.imsrgsolver.SetODETolerance(self.ode_tolerance)


  def init_modelspace(self):
    vs = self.valence_space
    if self.custom_valence_space:
      vs = self.custom_valence_space
    self.ms = ModelSpace(self.emax, self.ref, vs)
    self.ms.SetE3max(self.E3max)
    lmax = self.emax
    self.ms.SetLmax(lmax)
    self.ms.SetHbarOmega(self.hw)
    self.ms.SetTargetMass(self.A)


  def read_interaction(self, file2b, file3b):
    Hbare = Operator(self.ms,0,0,0,3)
    Hbare.SetHermitian()
    if self.file3_format == 'no2b':
      Hbare.ThreeBody.SetMode('no2b')
    if self.file3_precision == 'half':
      Hbare.ThreeBody.SetMode("no2bhalf")
    self.rw.ReadBareTBME_Darmstadt(self.file2b_directory+file2b, Hbare, 
                                self.file2e1max, 
                                self.file2e2max, 
                                self.file2lmax)
    Hbare.ThreeBody.ReadFile([self.file3b_directory+file3b], 
                                  [self.file3e1max, 
                                  self.file3e2max, 
                                  self.file3e3max
                                  ]
                                  )
    Hbare += Trel_Op(self.ms)
    return Hbare
  

  def read_interaction_combine_delta(self, LECs):
    #Array containing the 2b file name
    files2b_delta = files_2b = ["TwBME-HO_NN-only_DN2LO_ALL_0_bare_hw10_emax14_e2max28.me2j.gz"]
    files_2b.append("TwBME-HO_NN-only_DN2LO_Ct_1S0pp_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_Ct_1S0np_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_Ct_1S0nn_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_Ct_3S1_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_C_1S0_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_C_3P0_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_C_1P1_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_C_3P1_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_C_3S1_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_C_3S1_3D1_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_C_3P2_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_c1_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_c2_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_c3_bare_hw10_emax14_e2max28.me2j.gz")
    files_2b.append("TwBME-HO_NN-only_DN2LO_c4_bare_hw10_emax14_e2max28.me2j.gz")
    # LECs for the 2b part, the constant part need to be reweighted by the other LECs
    # as it is included in all the files by defaults in NuHamil (i.e. it is overcounted).
    LECs_2b = [1-np.sum(LECs[:-2])]
    for i in range(len(LECs[:-2])):
      LECs_2b.append(LECs[i])
    #Array containing the 3B files. 
    files_3b = ["NO2B_half_ThBME_3NFJmax15_c1_1_c3_0_c4_0_cD_0_cE_0_NonLocal4_394_IS_hw10_ms16_32_28.stream.bin"]
    files_3b.append("NO2B_half_ThBME_3NFJmax15_c1_0_c3_1_c4_0_cD_0_cE_0_NonLocal4_394_IS_hw10_ms16_32_28.stream.bin")
    files_3b.append("NO2B_half_ThBME_3NFJmax15_c1_0_c3_0_c4_1_cD_0_cE_0_NonLocal4_394_IS_hw10_ms16_32_28.stream.bin")
    files_3b.append("NO2B_half_ThBME_3NFJmax15_c1_0_c3_0_c4_0_cD_1_cE_0_NonLocal4_394_IS_hw10_ms16_32_28.stream.bin")
    files_3b.append("NO2B_half_ThBME_3NFJmax15_c1_0_c3_0_c4_0_cD_0_cE_1_NonLocal4_394_IS_hw10_ms16_32_28.stream.bin") 
    #LECs for the 3B part. Need to remove value from c3 due to convention
    LECs_3b = [LECs[11],LECs[13]-2.972246,LECs[14]+1.486123,LECs[15],LECs[16]]
    #Initialized the Hamiltonian operator  
    Hbare = Operator(self.ms,0,0,0,3)
    Hbare.SetHermitian()
    #Initialized a temporary Hamiltonian operator
    #It will be used to add the different components
    #to the Hamiltonian.
    Hbare_temp = Operator(self.ms,0,0,0,3)
    Hbare_temp.SetHermitian()
    #Set parameters for 3B part of the interaction
    if self.file3_format == 'no2b':
      Hbare.ThreeBody.SetMode('no2b')
      Hbare_temp.ThreeBody.SetMode('no2b')
    if self.file3_precision == 'half':
      Hbare.ThreeBody.SetMode("no2bhalf")
      Hbare_temp.ThreeBody.SetMode("no2bhalf")
    #Read constant part of the Hamiltonian and multiply it by the reweighted LEC
    #i.e. the one where we remove the double counting coming from the other files.
    self.rw.ReadBareTBME_Darmstadt(self.file2b_directory+files_2b[0], Hbare, 
                                self.file2e1max, 
                                self.file2e2max, 
                                self.file2lmax)
    Hbare *= LECs_2b[0]
    # Loops over all other LECs and read the associated file, multiply it by the LEC
    # Then erase the operator to save memory.
    for i,lec in enumerate(LECs_2b[1:]):
      self.rw.ReadBareTBME_Darmstadt(self.file2b_directory+files_2b[i+1], Hbare_temp,
                                self.file2e1max, 
                                self.file2e2max, 
                                self.file2lmax)
      Hbare_temp.ScaleTwoBody(lec)
      Hbare += Hbare_temp
      Hbare_temp.Erase()
    # Same Loop but for the 3B term. 
    for i,lec in enumerate(LECs_3b):
      Hbare_temp.ThreeBody.ReadFile([self.file3b_directory+files_3b[i]], 
                                    [self.file3e1max, 
                                    self.file3e2max, 
                                    self.file3e3max
                                    ]
                                    )
      Hbare_temp *= lec
      Hbare += Hbare_temp
      Hbare_temp.Erase()
    Hbare += Trel_Op(self.ms)
    return Hbare


  def print_estimatePT(self, HNO):
    #Give estimate with perturbation theory to make sure everything is ok
    print("Perturbative estimates of gs energy:")
    EMP2 = HNO.GetMP2_Energy()
    EMP2_3B = HNO.GetMP2_3BEnergy()
    print(f"EMP2 = {EMP2}")
    print(f"EMP2_3B =  {EMP2_3B}")
    Emp_3 = HNO.GetMP3_Energy()
    EMP3 = Emp_3[0]+Emp_3[1]+Emp_3[2]
    print(f"E3_pp = {Emp_3[0]}  E3_hh = {Emp_3[1]}  E3_ph = {Emp_3[2]} EMP3 =  {EMP3}")
    print(f"To 3rd order, E = {HNO.ZeroBody + EMP2 + EMP3 + EMP2_3B}")
    stdout.flush()  #So that print statement are in the right order


  def write_op_to_file(self, op, opname, extra=None):
    filename = f"{self.intfile}_{opname}"
    if extra:
      filename += f"_{extra}"
    filename += ".snt"
    if (op.GetJRank() == 0 and op.GetTRank() == 0 and op.GetParity() == 0):
      self.rw.WriteTokyo(op, filename, "op")
    else:
      self.rw.WriteTensorTokyo(filename,op)


  def evolve_operators(self):
    if len(self.opnames) != 0:
      for i,opname in enumerate(self.opnames):
          print(f"Starting to evolve {opname}:")
          op = OperatorFromString(self.ms, opname)
          if self.write_HO_ops:
            print(f"Writing HO operators to {self.output_dir}")
            self.write_op_to_file(op, opname, extra = "HO")
          op = self.hf.TransformToHFBasis(op)
          if self.write_HF_ops:
            op = op.DoNormalOrderingCore()
            print(f"Writing HF operators to {self.output_dir}")
            self.write_op_to_file(op, opname, extra = "HF")
            op = op.UndoNormalOrdering()
          op = op.DoNormalOrdering()
          op = self.imsrgsolver.Transform(op)
          op = op.UndoNormalOrdering()
          op = op.DoNormalOrderingCore()
          print( opname , 'zero body = ', op.ZeroBody)
          self.write_op_to_file(op, opname)
    if len(self.opfiles) != 0:
      #TODO add function to read operator from file
      print("Opeator from file not yet implemented.")
      print("Continuing...")


  def evolve_Hamiltonian(self, HNO):
    #Initialize the IMSRGSolver instance and set parameters
    self.imsrgsolver = IMSRGSolver(HNO)
    self.imsrgsolver.SetHin(HNO)
    self.imsrgsolver.SetReadWrite(self.rw)
    #Set parameters for the IMSRG solver
    self.set_imsrgsolver()

    #Decouple the core
    self.imsrgsolver.SetGenerator(self.core_generator)
    self.imsrgsolver.Solve()

    #Update IMSRG params for the decoupling of the valence-space
    self.imsrgsolver.SetSmax( 2*self.smax)
    self.imsrgsolver.SetGenerator(self.valence_space_generator)
    self.imsrgsolver.Solve()

    #### Get evolved Hamiltonian NO wrt the core
    HNO = self.imsrgsolver.GetH_s()
    HNO = HNO.UndoNormalOrdering()
    HNO = HNO.DoNormalOrderingCore()
    return HNO


  def run(self, file2b, file3b):
    #Initiate the ReadWrite class to access files
    self.rw = ReadWrite()

    #Create the model space for the nuclei
    self.init_modelspace()

    #Create the input hamiltonian from the 2b and 3b file
    Hbare = self.read_interaction(file2b, file3b)

    #Solve HFMBPT to obtain reference state
    self.hf = HFMBPT( Hbare )
    self.hf.Solve()
    HNO = self.hf.GetNormalOrderedH(2)

    #Give estimate with perturbation theory to make sure everything is ok
    self.print_estimatePT(HNO)
    
    #Do the IMSRG evolution of the Hamiltonian
    HNO = self.evolve_Hamiltonian(HNO)

    # Write things to disk
    self.output_dir = f"{self.output_directory_base}/{self.ref}/{self.label}/"
    Path(self.output_dir).mkdir(parents=True, exist_ok=True)
    print(f'Writing output file to {self.output_dir}')
    self.gen_filebase()
    self.intfile = f"{self.output_dir}/{self.filebase}"
    self.rw.WriteTokyo(HNO,self.intfile+".snt", "")
    stdout.flush()
    
    # Evolve operators 
    self.evolve_operators()
       
    # # #Write_kshell_jobs to be submitted after the imsrg has ran
    # for states in params["state_lists"]:
    #   kshl_r, f_diag_r = write_kshell_diag(params['path_to_kshell'], intfile+".snt", params['Nucl'], params['hw_truncation'], params['ph_truncation'], params['header'], gen_partition=True, states=states[0])
    #   if params['Nucl_daughter']:
    #     kshl_l, f_diag_l = write_kshell_diag(params['path_to_kshell'], intfile+".snt", params['Nucl_daughter'], params['hw_truncation'], params['ph_truncation'], params['header'], gen_partition=True, states=states[1])
  

  def run_combine_delta(self, LECs, sampleID):
    #Initiate the ReadWrite class to access files
    self.rw = ReadWrite()

    #Create the model space for the nuclei
    self.init_modelspace()

    #Create the input hamiltonian from the 2b and 3b file
    Hbare = self.read_interaction_combine_delta(LECs)

    #Solve HFMBPT to obtain reference state
    self.hf = HFMBPT( Hbare )
    self.hf.Solve()
    HNO = self.hf.GetNormalOrderedH(2)

    #Give estimate with perturbation theory to make sure everything is ok
    self.print_estimatePT(HNO)

    #Do the IMSRG evolution of the Hamiltonian
    HNO = self.evolve_Hamiltonian(HNO)

    # Write things to disk
    self.output_dir = f"{self.output_directory_base}/{self.ref}/{self.label}/"
    Path(self.output_dir).mkdir(parents=True, exist_ok=True)
    print(f'Writing output file to {self.output_dir}')
    self.gen_filebase(sampleID)
    self.intfile = f"{self.output_dir}/{self.filebase}"
    self.rw.WriteTokyo(HNO,self.intfile+".snt", "")
    stdout.flush()
    
    # Evolve operators 
    self.evolve_operators()


  def gen_filebase(self, sampleID = None):
    if not sampleID:
      self.filebase = f"{self.valence_space}_{self.label}_{self.ref}_{self.method}_e{self.emax}_E{self.E3max}_hw{self.hw}"
    else:
      self.filebase = f"{self.valence_space}_{self.label}_{sampleID}_{self.ref}_{self.method}_e{self.emax}_E{self.E3max}_hw{self.hw}"


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
    s=f"""{self.header}

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
{self.run_cmd} ./kshell.exe {self.fn_base}_{str_state}.input > log_{self.fn_base}_{str_state}.txt 2>&1 

rm -f tmp_snapshot_{Path(self.fn_ptn).name}_0_* tmp_lv_{Path(self.fn_ptn).name}_0_* {self.fn_base}_{str_state}.input 


./collect_logs.py log_*{self.fn_base}* > summary_{self.fn_base}.txt
cp summary_{self.fn_base}.txt {self.output_directory}"""
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
    s = f"""{self.header}


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
{self.run_cmd} ./transit.exe density_{self.filebase}_{self.Nucl_daughter}{state_string(self.state_list[1], self.A_daughter)}_{self.Nucl}{state_string(self.state_list[0], self.A)}.input > density_{self.filebase}_{self.Nucl_daughter}{state_string(self.state_list[1], self.A_daughter)}_{self.Nucl}{state_string(self.state_list[0], self.A)}.txt 2>&1
rm density_{self.filebase}_{self.Nucl_daughter}{state_string(self.state_list[1], self.A_daughter)}_{self.Nucl}{state_string(self.state_list[0], self.A)}.input"""
    fn_script = f"{self.scratch_directory}/density_{self.filebase}_{self.Nucl_daughter}{state_string(self.state_list[1], self.A_daughter)}_{self.Nucl}{state_string(self.state_list[0], self.A)}.sh"
    f = open(fn_script, "w")
    f.write(s)
    f.close()
    os.chmod(fn_script, 0o755)
    return fn_script




class kshell_tookit():

  def __init__(self, fn_snt, Nucl, state_list, Nucl_daughter=None, **kwargs):
    self.Nucl = Nucl
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
      
  def calc_diag(self, gen_partition=True, run = 'True'):
    if gen_parition:
      if str_state(state_list[-1], self.Nucl)[-1] == 'p': 
        parity = 1
      else:
        parity = -1
      self.gen_partition(parity)
    self.ket_diag = self.kshell_ket.gen_script()
    if self.Nucl != self.Nucl_daughter:
      if gen_parition:
        if str_state(state_list[-1], self.Nucl_daughter)[-1] == 'p': 
          parity = 1
        else:
          parity = -1
        self.gen_partition(parity, ket=False)
      self.bra_diag = self.kshell_bra.gen_script()
    


  def calc_density():
    pass

  def calc_opexpvals():
    if self.kshell_bra.A < self.kshell_ket.A or self.kshell_bra.Z < self.kshell_ket.Z:
      kshell_bra, kshell_ket = self.kshell_ket, self.kshell_bra