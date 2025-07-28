import os
#Get the unsername of the the user on the cluster
username = os.environ['USER']
#Finds the absolute path where imsrg_toolkit is installed
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#Paths to the interactions files
INTERACTION_2B_PATH = '/ceph/submit/data/group/ab-initio/me2j/' 
INTERACTION_3B_PATH = '/ceph/submit/data/group/ab-initio/me3j/'
