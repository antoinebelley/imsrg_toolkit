from pyIMSRG import *
from sys import stdout
import sys, os
from pathlib import Path
import numpy as np
from imsrg_toolkit.settings import *
from imsrg_toolkit.imsrg_params import *





class Imsrg(ImsrgParams):
  def __init__(self, **kwargs):
    # #TODO Clean this up to use the ImsrgParams class in utils.py
    # #### Here are the default parameters for the imsrg###
    # ### TODO add all IMSRG parameters in the params
    # #Paths to different directories that are used
    # #TODO update those from a config file
    # self.scratch_directory = f'/work/submit/{username}/work/imsrg/'
    # self.output_directory_base = f'/home/submit/{username}/results/'
    # self.file2b_directory = '/ceph/submit/data/group/ab-initio/me2j/'
    # self.file3b_directory  = '/ceph/submit/data/group/ab-initio/me3j/'
    

    # # Model space parameters
    # self.A = 6
    # self.emax = 2
    # self.E3max = 6
    # self.hw = 10
    # self.ref = 'He6'
    # self.valence_space = 'p-shell'
    # self.custom_valence_space = None

    # #2B interaction parameters
    # self.label = 'SampleDelta'
    # self.file2e1max = 14
    # self.file2e2max = 28
    # self.file2lmax = 14

    # #3B interaction parameters
    # self.file3e1max = 16
    # self.file3e2max = 32
    # self.file3e3max = 28
    # self.file3_format = 'no2b'
    # self.file3_precision = 'half'

    # #IMSRG solver parameters
    # self.method = 'magnus'
    # self.denominator_partitioning = 'Epstein_Nesbet'
    # self.eta_criterion = 1e-6
    # self.smax = 500
    # self.dsmax = 0.5
    # self.ds0 = 0.5
    # self.denominator_delta = 0
    # self.domega = 0.2
    # self.omega_norm_max = 0.25
    # self.ode_tolerance = 1e-6
    # self.core_generator = 'atan'
    # self.valence_space_generator = 'shell-model-atan'

    # #Operators parameters
    # self.opfiles = []
    # self.opnames = []
    # self.write_HO_ops = True
    # self.write_HF_ops = True

    # #If dictionay is given, update the attributes using the
    # #dictionary keys and values.
    # self.update_params(**kwargs)

    # self.output_dir = f"{self.output_directory_base}/{self.ref}/{self.label}/"
    super().__init__(**kwargs)
    Path(self.scratch_directory).mkdir(parents=True, exist_ok=True)
    Path(self.output_dir).mkdir(parents=True, exist_ok=True)


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


  def read_interaction_combine_delta(self):
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
    LECs_2b = [1-np.sum(self.LECs[:-2])]
    for i in range(len(self.LECs[:-2])):
      LECs_2b.append(self.LECs[i])
    #Array containing the 3B files.
    files_3b = ["NO2B_half_ThBME_3NFJmax15_c1_1_c3_0_c4_0_cD_0_cE_0_NonLocal4_394_IS_hw10_ms16_32_28.stream.bin"]
    files_3b.append("NO2B_half_ThBME_3NFJmax15_c1_0_c3_1_c4_0_cD_0_cE_0_NonLocal4_394_IS_hw10_ms16_32_28.stream.bin")
    files_3b.append("NO2B_half_ThBME_3NFJmax15_c1_0_c3_0_c4_1_cD_0_cE_0_NonLocal4_394_IS_hw10_ms16_32_28.stream.bin")
    files_3b.append("NO2B_half_ThBME_3NFJmax15_c1_0_c3_0_c4_0_cD_1_cE_0_NonLocal4_394_IS_hw10_ms16_32_28.stream.bin")
    files_3b.append("NO2B_half_ThBME_3NFJmax15_c1_0_c3_0_c4_0_cD_0_cE_1_NonLocal4_394_IS_hw10_ms16_32_28.stream.bin")
    #LECs for the 3B part. Need to remove value from c3 due to convention
    LECs_3b = [self.LECs[11],self.LECs[13]-2.972246,self.LECs[14]+1.486123,self.LECs[15],self.LECs[16]]
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
      if op.GetJRank() == 0:
        op.MakeReduced()
      self.rw.WriteTensorTokyo(filename,op)


  def evolve_operators(self, HF=False):
    self.ops = []
    self.op_strings = [] #Name at the end of the snt file
    if len(self.opnames) != 0:
      for opname in self.opnames:
        print(f"Starting to evolve {opname}:")
        op = OperatorFromString(self.ms, opname)
        self.ops.append(op)
        self.op_strings.append(opname)
    if len(self.opfiles) != 0:
      for opfile in self.opfiles:
        op = self.rw.ReadOperator2b_Miyagi(opfile[0], self.ms)
        self.ops.append(op)
        self.op_strings.append(opfile[1])
    for op, name in zip(self.ops, self.op_strings):
      #Write HO operators
      if self.write_HO_ops:
        op = op.DoNormalOrderingCore()
        print(f"Writing HO operators to {self.output_dir}")
        self.write_op_to_file(op, name, extra = "HO")
        op = op.UndoNormalOrderingCore()
      #Transform operators to HF basis
      op = self.hf.TransformToHFBasis(op)
      #Write HF operator to file.
      if self.write_HF_ops or HF:
        op = op.DoNormalOrderingCore()
        print(f"Writing HF operators to {self.output_dir}")
        self.write_op_to_file(op, name, extra = "HF")
        if HF:
          return
        op = op.UndoNormalOrderingCore()
      op = op.DoNormalOrdering()
      op = self.imsrgsolver.Transform(op)
      op = op.UndoNormalOrdering()
      op = op.DoNormalOrderingCore()
      print( opname , 'zero body = ', op.ZeroBody)
      self.write_op_to_file(op, name)


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

  def run(self, file2b, file3b, HF=False):
    #Initiate the ReadWrite class to access files
    self.rw = ReadWrite()

    #Create the model space for the nuclei
    self.init_modelspace()

    #Create the input hamiltonian from the 2b and 3b file
    Hbare = self.read_interaction(file2b, file3b)
    if self.BetaCM!=0:
      Hbare += self.BetaCM*OperatorFromString(self.ms, f'HCM_{self.hwBetaCM}')

    #Solve HFMBPT to obtain reference state
    self.hf = HFMBPT( Hbare )
    self.hf.Solve()
    HNO = self.hf.GetNormalOrderedH(2)
    if self.BetaCM != 0:
      HNO -= self.BetaCM * 1.5*self.hwBetaCM



    #If we only want the HF results, stop the calculation of the interacion here
    if HF:
      HNO = HNO.DoNormalOrderingCore()
      print(f'Writing output file to {self.output_dir}')
      self.rw.WriteTokyo(HNO,self.intfile+"_HF.snt", "")
      stdout.flush()
    else:
      #Give estimate with perturbation theory to make sure everything is ok
      self.print_estimatePT(HNO)

      #Do the IMSRG evolution of the Hamiltonian
      HNO = self.evolve_Hamiltonian(HNO)

      # Write things to disk
    
      print(f'Writing output file to {self.output_dir}')
      self.intfile = f"{self.output_dir}/{self.filebase}"
      self.rw.WriteTokyo(HNO,self.intfile+".snt", "")
      stdout.flush()

    # Evolve operators
    self.evolve_operators(HF=HF)


  def run_combine_delta(self, HF=False):
    #Initiate the ReadWrite class to access files
    self.rw = ReadWrite()

    #Create the model space for the nuclei
    self.init_modelspace()

    #Create the input hamiltonian from the 2b and 3b file
    Hbare = self.read_interaction_combine_delta()
    if self.BetaCM!=0:
      Hbare += self.BetaCM*OperatorFromString(self.ms, f'HCM_{self.hwBetaCM}')

    #Solve HFMBPT to obtain reference state
    self.hf = HFMBPT( Hbare )
    self.hf.Solve()
    HNO = self.hf.GetNormalOrderedH(2)
    if self.BetaCM != 0:
      HNO -= self.BetaCM * 1.5*self.hwBetaCM
    if HF:
      HNO = HNO.DoNormalOrderingCore()
      Path(self.output_dir).mkdir(parents=True, exist_ok=True)
      print(f'Writing output file to {self.output_dir}')
      self.intfile = f"{self.output_dir}/{self.filebase}"
      self.rw.WriteTokyo(HNO,self.intfile+"_HF.snt", "")
      stdout.flush()
    
    else:
      #Give estimate with perturbation theory to make sure everything is ok
      self.print_estimatePT(HNO)

      #Do the IMSRG evolution of the Hamiltonian
      HNO = self.evolve_Hamiltonian(HNO)

      # Write things to disk
      self.output_dir = f"{self.output_directory_base}/{self.ref}/{self.label}/"
      Path(self.output_dir).mkdir(parents=True, exist_ok=True)
      print(f'Writing output file to {self.output_dir}')
      self.intfile = f"{self.output_dir}/{self.filebase}"
      self.rw.WriteTokyo(HNO,self.intfile+".snt", "")
      stdout.flush()

    # Evolve operators
    self.evolve_operators(HF=HF)



# class PvImsrg(Imsrg):
#   def __init__(staged=False, **kwargs):
#     super.__init__(**kwargs)
#     staged = staged

#   def run_anapole_moment(self, file2b, file3b HF=False, scale=1000):
#     #Initiate the ReadWrite class to access files
#     self.rw = ReadWrite()

#     #Create the model space for the nuclei
#     self.init_modelspace()

#     #Create the input hamiltonian from the 2b and 3b file
#     Hbare = self.read_interaction(file2b, file3b, HF = HF)
#     if self.BetaCM!=0:
#       Hbare += self.BetaCM*OperatorFromString(self.ms, f'HCM_{self.hwBetaCM}')

#     #Solve HFMBPT to obtain reference state
#     self.hf = HFMBPT( Hbare )
#     self.hf.Solve()
#     HNO = self.hf.GetNormalOrderedH(2)
#     if self.BetaCM != 0:
#       HNO -= self.BetaCM * 1.5*self.hwBetaCM
    
#     VPNC_path = f"{INTERACTION_2B_PATH}/PVTCint_DDH_best-LO-NonLocal2-500_TwBME-HO_NN-only_N3LO_EM500_bare_hw15_emax12_e2max24.me2j.gz"
#     VPNC = self.rw.ReadOperator2b_Miyagi(VPNC_path, self.ms)
#     VPNC.SetAntiHermitian()
#     VPNC *= scale



#     #If we only want the HF results, stop the calculation of the interaction here
#     if HF:
#       HNO = HNO.DoNormalOrderingCore()
#       print(f'Writing output file to {self.output_dir}')
#       self.rw.WriteTokyo(HNO,self.intfile+"_HF.snt", "")
#       stdout.flush()
#     else:
#       #Give estimate with perturbation theory to make sure everything is ok
#       self.print_estimatePT(HNO)

#       #Do the IMSRG evolution of the Hamiltonian
#       HNO = self.evolve_Hamiltonian(HNO)

#       # Write things to disk
    
#       print(f'Writing output file to {self.output_dir}')
#       self.intfile = f"{self.output_dir}/{self.filebase}"
#       self.rw.WriteTokyo(HNO,self.intfile+".snt", "")
#       stdout.flush()

#     # Evolve operators
#     self.evolve_operators(HF=HF)