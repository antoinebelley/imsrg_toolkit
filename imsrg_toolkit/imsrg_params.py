from imsrg_toolkit.settings import *
from dataclasses import dataclass, field

@dataclass(frozen=False)
class ImsrgParams():
  #Paths for inputs and outputs
  scratch_directory: str = f'/work/submit/{username}/work/imsrg/'
  output_directory_base: str = f'/home/submit/{username}/results/'
  file2b_directory: str = INTERACTION_2B_PATH
  file3b_directory: str = INTERACTION_3B_PATH
  #Model space parameters
  A: int = 6
  emax: int = 2
  E3max: int = 6
  hw: int = 10
  ref: str = 'He6'
  valence_space: str = 'p-shell'
  custom_valence_space: str = None 
   #2B interaction parameters 
  label: str = 'SampleDelta'
  SampleID: str = None #This is for the interactions samples only
  LECs: list[float] = field(default_factory=list)
  file2e1max: int = 14
  file2e2max: int = 28
  file2lmax: int = 14    
  #3B interaction parameters
  file3e1max: int = 16
  file3e2max: int = 32
  file3e3max: int = 28
  file3_format: str = 'no2b'
  file3_precision: str = 'half'
  #Parameters for the BetaCM
  BetaCM: float = 0.0
  hwBetaCM: float = hw  # Negative value means use the frequency
  #IMSRG solver parameters
  basis : str = 'HF'
  method: str = 'magnus'
  denominator_partitioning: str = 'Epstein_Nesbet'
  eta_criterion: float = 1e-6
  smax: int = 500
  dsmax: float = 0.5
  ds0: float = 0.5
  denominator_delta: float = 0
  denominator_delta_orbit: str = None
  domega: float = 0.2
  omega_norm_max: float = 0.25
  ode_tolerance: float = 1e-6
  core_generator: str = 'atan'
  valence_space_generator: str = 'shell-model-atan'
  #Operators parameters
  opfiles: list[str] = field(default_factory=list)
  opnames: list[str] = field(default_factory=list)
  write_HO_ops: bool = True
  write_HF_ops: bool = True
  #Strings for submission
  run_cmd: str = None
  header: str = None


  def gen_filebase(self):
    if not self.SampleID:
      self.filebase = f"{self.valence_space}_{self.label}_{self.ref}_{self.method}_e{self.emax}_E{self.E3max}_hw{self.hw}"
    else:
      self.filebase = f"{self.valence_space}_{self.label}_{self.SampleID}_{self.ref}_{self.method}_e{self.emax}_E{self.E3max}_hw{self.hw}"


  def __post_init__(self):
    #generate the file name for the snt file
    self.output_dir = f"{self.output_directory_base}/{self.ref}/{self.label}/"
    self.gen_filebase()
    self.intfile = f"{self.output_dir}/{self.filebase}"