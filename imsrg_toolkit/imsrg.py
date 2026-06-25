from pyIMSRG import *
from sys import stdout
import sys, os
from pathlib import Path
import numpy as np
from imsrg_toolkit.settings import *
from imsrg_toolkit.imsrg_params import *





class Imsrg(ImsrgParams):
  def __init__(self, **kwargs):
    #Initiate the class based on the abstract data class
    #which contains all of the parameters.
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
    self.imsrgsolver.SetDenominatorDeltaOrbit(self.denominator_delta_orbit)
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
    if len(self.opnames_decay) != 0:
      for opname in self.opnames_decay:
        print(f"Starting to evolve {opname}:")
        op = OperatorFromString(self.ms, opname)
        self.ops.append(op)
        self.op_strings.append(opname)
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
      #If 3f2, Factorized double commutators should arleady be set
      #to True.
      if self.approx == 'imsrg3f2' and op.GetJRank() == 0:
        Commutator.FactorizedDoubleCommutator.SetUse_1b_Intermediates(True)
        Commutator.FactorizedDoubleCommutator.SetUse_2b_Intermediates(True)
      else:
        Commutator.FactorizedDoubleCommutator.SetUse_1b_Intermediates(False)
        Commutator.FactorizedDoubleCommutator.SetUse_2b_Intermediates(False)
      op = self.imsrgsolver.Transform(op)
      op = op.UndoNormalOrdering()
      op = op.DoNormalOrderingCore()
      print( opname , 'zero body = ', op.ZeroBody)
      self.write_op_to_file(op, name)


  def evolve_Hamiltonian(self, HNO, renormal_order=True):
    #Initialize the IMSRGSolver instance and set parameters
    self.imsrgsolver = IMSRGSolver(HNO)
    self.imsrgsolver.SetHin(HNO)
    self.imsrgsolver.SetReadWrite(self.rw)
    #Set parameters for the IMSRG solver
    self.set_imsrgsolver()

    if self.approx == 'imsrg3f2':
      # We only include the factorized commutators with 1b intermediates
      # during the flow, because they are cheper and the 2b intermediates
      # have a much less significant influence on Hod, so they can safely
      # be added in at the end of the flow.
      # We also use the Hunter-Gatherer mode, which works best if we
      # take OmegaNormMax to be smaller than 1
      self.imsrgsolver.SetHunterGatherer(True)
      self.imsrgsolver.SetOmegaNormMax(0.1)
      BCH.SetUseFactorizedCorrection(True)
      Commutator.FactorizedDoubleCommutator.SetUse_1b_Intermediates(True)
      Commutator.FactorizedDoubleCommutator.SetUse_2b_Intermediates(False)

    #Decouple the core
    self.imsrgsolver.SetGenerator(self.core_generator)
    self.imsrgsolver.Solve()

    #Update IMSRG params for the decoupling of the valence-space
    if self.valence_space == self.ref:
      Hs = self.imsrgsolver.GetH_s()
      triples = 0
      if self.approx == 'imsrg3f2':
          Commutator.FactorizedDoubleCommutator.SetUse_1b_Intermediates(True)
          Commutator.FactorizedDoubleCommutator.SetUse_2b_Intermediates(True)
          Hs = self.imsrgsolver.Transform(HNO)
          triples = self.imsrgsolver.CalculatePerturbativeTriples()
      Eimsrg = Hs.ZeroBody
      print('IMSRGEnergy = {:.6f}  +  {:.6f}  = {:.6f}'.format(Eimsrg,triples,Eimsrg+triples))

    else:
      self.imsrgsolver.SetSmax( 2*self.smax)
      self.imsrgsolver.SetGenerator(self.valence_space_generator)
      self.imsrgsolver.Solve()
      #### Get evolved Hamiltonian NO wrt the core
      Hs = self.imsrgsolver.GetH_s()
      if self.approx == 'imsrg3f2':
        Commutator.FactorizedDoubleCommutator.SetUse_1b_Intermediates(True)
        Commutator.FactorizedDoubleCommutator.SetUse_2b_Intermediates(True)
        Hs = self.imsrgsolver.Transform(HNO)
        triples = self.imsrgsolver.CalculatePerturbativeTriples()
        print('Adding triples correction  = {:.6f}'.format(triples),flush=True)
        Hs.ZeroBody += triples

    if Hs.GetParticleRank()>2:
      Hs.ThreeBody.Erase()

    if renormal_order:
      Hs = Hs.UndoNormalOrdering()
      Hs = Hs.DoNormalOrderingCore()
    return Hs


  def run(self, file2b, file3b, HF=False, verbose=True):
    if verbose == True:
      print({k: v for k, v in self.__dict__.items() if not k.startswith("__")}, flush=True)
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
    self.hf.PrintSPEandWF()
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



  
class PvImsrg(Imsrg):
  def __init__(staged=False, **kwargs):
    super().__init__(**kwargs)


  def op_HF_NO(self, op):
    # Transform operators to HF basis
    if op.GetJRank() == 0 and (op.GetTRank() != 0 or op.GetParity() != 0):
      op.MakeNotReduced()
    op = self.hf.TransformToHFBasis(op)
    # Write HF operator to file.
    op = op.DoNormalOrdering()
    if op.GetJRank() == 0 and (op.GetTRank() !=0 or op.GetParity()!=0):
      op.MakeReduced()
    return op
  

  def add_operator(self, op, opname):
    self.ops.append(op)
    self.op_strings.append(opname)


  def __write_op_to_file__(self, fn_base, op, opname=None, path=None, extra=None):
    if path:
      fn = f'{path}/{fn_base}'
    else:
      fn = f'{fn_base}'
    if opname:
      fn += f"_{opname}"
    if extra:
      fn += f"_{extra}"
    fn += ".snt"
    op_copy = Operator(op)
    if op_copy.GetJRank() == 0 and (op_copy.GetTRank() != 0 or op_copy.GetParity() != 0):
      op_copy.MakeNotReduced()
    op_copy = op_copy.UndoNormalOrdering()
    if op_copy.GetJRank() == 0 and (op_copy.GetTRank() != 0 or op_copy.GetParity() != 0):
      op_copy.MakeReduced()
    if op_copy.GetJRank() == 0 and (op_copy.GetTRank() != 0 or op_copy.GetParity() != 0):
      op_copy.MakeNotReduced()
    op_copy = op_copy.DoNormalOrderingCore()
    if op_copy.GetJRank() == 0 and (op_copy.GetTRank() != 0 or op_copy.GetParity() != 0):
      op_copy.MakeReduced()
    if op_copy.GetJRank() == 0 and op_copy.GetTRank() == 0 and op_copy.GetParity() == 0:
      self.rw.WriteTokyo(op_copy, fn, "")
    else:
      self.rw.WriteTensorTokyo(fn, op_copy)
    return fn


  def write_to_file(self, HNO, fn_base = None, VPV=None, path=None, extra=None):
    files = []
    if not fn_base:
      vs = self.vs
      fn_base = f"{vs}_{self.ref}_hw{self.hw}_emax{self.emax}_E{self.E3max}_s{self.solver_params['smax']}_ds{self.solver_params['dsmax']}"
      if self.BetaCM != 0:
        fn_base += f'_BCM{self.BetaCM}'
    print(fn_base)
    files.append(self.__write_op_to_file__(fn_base, HNO, path=path, extra=extra))
    if VPV:
      files.append(self.__write_op_to_file__(fn_base, VPV, opname = "VPV", path=path, extra=extra))
    for op, opname in zip(self.ops, self.op_strings):
      files.append(self.__write_op_to_file__(
          fn_base, op, opname=opname, path=path, extra=extra))
    return files


  def get_HNO(self,  H0):
    self.hf = HartreeFock(H0)
    self.hf.Solve()
    self.hf.PrintSPEandWF()
    HNO = self.hf.GetNormalOrderedH(2)
    return HNO
  

  def operators_to_HF(self):
    self.ops = [self.op_HF_NO(op) for op in self.ops]


  def evolve_operators_PV(self):
    self.ops = [self.imsrgsolver.Transform(
        op, op_PV) for op, op_PV in zip(*[iter(self.ops)]*2)]
    self.ops = [item for t in self.ops for item in t]


  def evolve_operators(self):
    self.ops = [self.imsrgsolver.Transform(op)for op in self.ops]


  def evolve_Hamiltonian_PV(self, HNO, VPV, core_generator = None, valence_space_generator=None):
    if not core_generator:
      core_generator = self.core_generator
    if not valence_space_generator:
      valence_space_generator = self.valence_space_generator
    print(core_generator)
    print(valence_space_generator)
    stdout.flush()
    # Initialize the IMSRGSolver instance and set parameters
    self.imsrgsolver = IMSRGSolverPV(HNO, VPV)
    self.imsrgsolver.SetReadWrite(self.rw)
    #Set parameters for the IMSRG solver
    self.set_imsrgsolver()

    #Decouple the core
    self.imsrgsolver.SetGeneratorPV(core_generator)
    self.imsrgsolver.Solve_magnus_euler()

    #Update IMSRG params for the decoupling of the valence-space
    self.imsrgsolver.SetSmax( 2*self.smax)
    self.imsrgsolver.SetGeneratorPV(valence_space_generator)
    self.imsrgsolver.Solve_magnus_euler()

    #### Get evolved Hamiltonian
    HNO = self.imsrgsolver.GetH_s()
    VPV = self.imsrgsolver.GetVPT_s()

    return HNO, VPV


  def run_anapole(self, file2b, file3b, scale=1000, staged=True):
    #Initiate the ReadWrite class to access files
    self.rw = ReadWrite()

    #Create the model space for the nuclei
    self.init_modelspace()

    #Create the input hamiltonian from the 2b and 3b file
    Hbare = self.read_interaction(file2b, file3b)
    if self.BetaCM!=0:
      Hbare += self.BetaCM*OperatorFromString(self.ms, f'HCM_{self.hwBetaCM}')

    #Create anapole operators
    Anapole = OperatorFromString(self.ms, 'Anapole')
    Anapolepp = Operator(self.ms, 1, 0, 0, 2)
    self.add_operator(Anapolepp, "Anapolepp")
    self.add_operator(Anapole, "Anapole")

    #Read the VPV interaction
    VPNC_path = f"{INTERACTION_2B_PATH}/PVTCint_DDH_best-LO-NonLocal2-500_TwBME-HO_NN-only_N3LO_EM500_bare_hw{self.hw}_emax12_e2max24.me2j.gz"
    VPNC = self.rw.ReadOperator2b_Miyagi(VPNC_path, self.ms)
    VPNC.SetAntiHermitian()
    VPNC *= scale

    #Solve HFMBPT to obtain reference state
    HNO = self.get_HNO(Hbare)
    if self.BetaCM != 0:
      HNO -= self.BetaCM * 1.5*self.hwBetaCM
    #Give estimate with perturbation theory to make sure everything is ok
    self.print_estimatePT(HNO)
    
    VPNC = self.op_HF_NO(VPNC)
    self.operators_to_HF()
    
    if staged:
      #We first do the decoupling using the standard IMSRG
      HNO = self.evolve_Hamiltonian(HNO, renormal_order=False)
      VPNC = self.imsrgsolver.Transform(VPNC)
      self.evolve_operators()
    
    HNO, VPNC = self.evolve_Hamiltonian_PV(HNO, VPNC, core_generator='white', valence_space_generator='shell-model')
    self.evolve_operators_PV()
    
    self.output_dir = f"{self.output_directory_base}/{self.ref}/{self.label}/"
    Path(self.output_dir).mkdir(parents=True, exist_ok=True)
    print(f'Writing output file to {self.output_dir}')
    stdout.flush()
    fn_snt, fn_VPNC, fn_Anapolepp, fn_Anapole = self.write_to_file(
      HNO, VPV=VPNC, path=os.path.abspath(self.output_dir), fn_base=self.filebase)
    return fn_snt, fn_Anapolepp