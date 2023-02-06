#! /usr/bin/python3

################################
####### Slurm Variables ########
################################

# Number of Segments (usually in nanoseconds):
Segments="50"
# Number of CPUs to be used in parallel (sum of all replicas)
ncpus="16"

## Choose the list of machines (partition)
## More than one partition can be supplied, separated by commas and without spaces
## GPU:   All machines with GPU
## MD:    All machines available for MD without GPU
## MD32:  List of machines with 32 cores/threads
## MD24:  List of older machines with 24 cores/threads
Partition="MD32c,MD32f"
#
## Requeue: the job is restarted in another node in case the assigned one reboots
## Options are 1 or 0
requeue=1
#
## Email Address to receive notification of completion:
## Multiple email addresses can be provided separated by comma.
## If left blank, no notification will be sent.
Email=""
#
# Files
# File Name: "pHREprefix_sn.extension" ex: pHRE_GLU.mdp
sn="penta"
REUSprefix="REUS"
USprefix="US"


#########################
###### REUS Inputs ######
#########################
#Reps: must be list of the number of reps to be done 
reps        = ['1','2','3']

# Directory Name: "prefix + pHs"
# pHs values :: MUST BE LIST of pH values
# pH values :: MUST BE STRINGS equal do the name of pH dirs
pHprefix     = 'pH'
pHs        = ['05.00']

# US Directories Name: "USprefix + USs"
# pHs values :: MUST BE LIST of pH values
# pH values :: MUST BE STRINGS equal do the name of pH dirs
USs = {
    '-10': 1000,
    '-08': 1000,
    '-06': 1000,
    '-04': 1000,
    '-02': 1000,
    '+00': 1000,
    '+02': 1000,
    '+04': 1000,
    '+06': 1000,
    '+08': 500,
    '+10': 500,
    '+12': 500,
    '+14': 500,
    '+16': 500,
    '+18': 500,
    '+20': 500,
}
              

# Files
# File Name: "pHREprefix_sn.extension" ex: pHRE_GLU.mdp
sn         = 'penta'
REUSprefix = 'REUS'

# Required Files Location
# Attention: pHmdp file should be prepared to run ONE cphmd segment
sts        = '/home/noliveira/CpHMD/St-54a7_Fit_DelPHi'  ### path to the ST files
gmx        = '/gromacs/gromacs-5.1.5_pH_I/bin/gmx' ### path to the desired gromacs  
delphi     = '/home/noliveira/CpHMD/DelphiTools_v2.0'  ### path to the Delphitools 
cphmd      = '/home/noliveira/CpHMD/CpH-MD_v2.1_US_verlet/scripts/CpHMD.sh' ### path to the CpHMD script
cphdir     = '/home/noliveira/CpHMD/CpH-MD_v2.1_US_verlet'  ### path to the CpHMD directory
pHmdp      = 'USTEMPLATE/' + sn + '_000.pHmdp' ### path to the pHmpd file
mdp        = 'USTEMPLATE/' + sn + '.mdp'  ### path to the mdp file
pbp        = 'USTEMPLATE/' + sn + '.pbp'  ### path to the pbp file
sites      = 'USTEMPLATE/' + sn + '.sites'  ### path to the sites file
top        = 'USTEMPLATE/' + sn + '.top'  ### path to the top file
pdb        = 'USTEMPLATE/' + sn + '.pdb'  ### path to the pdb file
ndx        = 'USTEMPLATE/' + sn + '.ndx'  ### path to the ndx file
### If you have gros in a separate folder
gro        = 'USTEMPLATE/GRO/' + sn  ### path to the gro file
# Location of the Gro Whole file needed for for pdb2gmx 
# not to break the special bonds (FULL PATH):
GroW       = 'USTEMPLATE/GRO/' + sn   ### path to the GroW file
# Location of the rules file (FULL PATH):
Rules      = '/home/noliveira/CpHMD/CpH-MD_v2.1_US_verlet/lipids_54A7.rules'  ### path to the rules files

# RE Input Parameters
effsteps   = 10000
nstxtcout  = 5000
dt         = 0.002
ncycles    = 50     #number of cycles to make 1 ns: 1000/effsteps*dt (50)
REUSfreq   = 1      
T          = 310.0  #should be written as float 

# US Input Parameters
pullgroup2  = 'memb'  #Reference group
pullgroup1  = 'asp'  #Group where the force will be applied
pullvec   = '0.0 0.0 1.0' #axys where the force will be applied

