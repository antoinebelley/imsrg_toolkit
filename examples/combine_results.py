import pandas as pd
import glob
import re
import pandas as pd
import os, sys
import numpy as np
HOME = os.path.expanduser("~")
sys.path.append(HOME)

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)





A = 24
Z = 13
N = A-Z
Nucl = f"Al{A}"
directory = f"/home/submit/abelley/results/{Nucl}/SampleDelta/"
output_file = f"/work/submit/abelley/results/{Nucl}_radii.csv"



LEC_labels = ['Ct1S0pp',
                 'Ct1S0np',
                 'Ct1S0nn',
                 'Ct3S1',
                 'C1S0',
                 'C3P0',
                 'C1P1',
                 'C3P1',
                 'C3S1',
                 'CE1',
                 'C3P2',
                 'c1',
                 'c2',
                 'c3',
                 'c4',
                 'cD',
                 'cE'
              ]


col_names = ['Nucl bra',
            'J bra',  
            'P bra',
            'n bra',
            'Energy bra', 
            'Nucl ket',
            'J ket',  
            'P ket',
            'n ket',
            'Energy ket',
            'Zero',
            'One',
            'Two',
            'Rch'
            ]

col_labels = LEC_labels + col_names 

class IMSRGResultsDF():
  """
  Data frame object containing all the reuslts up to date for the decay.
  """
  def __init__(self, file = None):
    self.intialize_df(file)

  def create_df(self):
    #Creates an empty data frame so that we can add data to it as we go
    names = ['Sample', 'emax',] 
    empty_index = pd.MultiIndex.from_tuples([], names=names)
    columns = col_labels
    df = pd.DataFrame(index=empty_index, columns=columns)
    self.df = df

  def intialize_df(self, file = None):
    if file == None:
      #Creates an empty data frame so that we can add data to it as we go
      self.create_df()
    else:
      try:
        self.df = pd.read_csv(file, dtype={'Sample': int, 'emax':int})
        names = ['Sample', 'emax',] 
        self.df.set_index(names, inplace=True)
        self.df.sort_index(inplace=True)
        print("Read df from file.")
      except:
        print("Output csv file not found. Creating new dataframe.")
        self.create_df()



  def add_to_dataframe(self, Nucl, Sample, emax, data):
    index = (Sample, emax)
    for key in data.keys():
      if isinstance(data[key], str):
        if key != "fn_op":
          data[key] = data[key].strip()
      try:
        if self.df.loc[index, key] == data[key]:
          continue
      except KeyError:
        pass
      self.df.loc[index, key] = data[key]
    self.df = self.df.fillna('')
    self.df.drop_duplicates(inplace=True)
    self.df.sort_index(inplace=True)
    

  def to_csv(self, file):
    # self.df = self.df.sort_values(['Sample'])
    self.df.sort_index(inplace=True)
    self.df.drop_duplicates(inplace=True)
    self.df.to_csv(file)

  def __str__(self):
    return self.df.__str__()



def Rp2_to_Rch2( Rp2, Z, N, CODATA=True ):
    """
    inputs:
        Rp2: mean squared point proton radius
        Z: proton number
        N: neutron number
    output:
        mean squared charge radii
    """
    if(CODATA):
        rcp2 = 0.8783**2 # CODATA
        rcn2 = -0.115    # CODATA
    else:
        rcp2 = 0.709     # Nature 466, 213 (2010).
        rcn2 = -0.106    # Phys. Rev. Lett. 124, 082501

    DF = 0.033
    return Rp2 + rcp2 + N/Z * rcn2 + DF


#Function to read the 0b term from the IMSRG
def read_0b_energy(fn):
  fn = fn.replace(".csv", ".snt")
  f = open(fn, 'r')
  lines = f.readlines()
  zerob = float(re.findall(r'[+-]?[0-9]+.+[0-9]', lines[4])[0])
  return zerob




LECs = pd.read_csv("/work/submit/abelley/imsrg_toolkit/data/8000Samples.txt", 
                     usecols = list(range(18)))
LECs = LECs.set_index("SampleID")

old_LECs_duplicate = [40383, 51445,  66697,  97836, 144298, 163571, 172911, 215846, 237612, 252207,
  280489, 28255, 388492, 444839, 537066, 546455, 556707, 557713, 589500, 647468,
  690146, 704538, 708334, 709194, 732716, 777891, 829596]


failed_samples = []

results = IMSRGResultsDF(file=output_file)

count = 0

for name in glob.glob(f"{directory}*.csv"):
  # if count>=1: break
  # try:
    fn = name.replace(directory,'')
    params = fn.split("_")
    for param in params:
      if param.isnumeric():
        SampleId = int(param)
      elif param.startswith('e'):
        emax = param.replace('e','')
        emax = int(emax)
    if SampleId in old_LECs_duplicate:
      continue
    LEC_cols = LECs.loc[SampleId].to_dict()
    dict_cols = pd.read_csv(name, dtype={
      'Nucl bra': str,
      'J bra' : str,  
      'P bra': str,
      'n bra': str,
      'Energy bra': float, 
      'Nucl ket': str,
      'J ket' : str,  
      'P ket': str,
      'n ket': str,
      'Energy ket': float,
      'Zero': float,
      'One': float,
      'Two': float
    }, usecols=[i for i in range(15)][2:] )
    dict_cols = dict_cols.rename(columns=lambda x: x.strip())
    dict_cols = dict_cols.iloc[0].to_dict()
    dict_cols.update(LEC_cols)
    Rp2 = dict_cols['Zero'] + dict_cols['One'] + dict_cols['Two']

    dict_cols['Rch'] = np.sqrt(Rp2_to_Rch2(Rp2,Z,N))
    results.add_to_dataframe(Nucl, SampleId, emax, dict_cols)
  # except:
  #   failed_samples.append(SampleId)
   
print(failed_samples) 
  
results.to_csv(output_file)