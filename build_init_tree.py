#! /usr/bin/python3.8

#####################
# Imports
###
import os
import subprocess
from REUS_input import *

#####################
# Functions
###
def create_dirs(REUSprefix, pHprefix, USprefix, USs, sn, mdp, pbp, sites, cphmd, gro, pH, rep):
  """
  Creates the folders for the different pH values and then creates each US 
  folder inside, filling those with the files needed to run the input:
    REUSprefix is str used to create filenames (prefix_sysname.extension)
    pHprefix is str of pH directory names prefix (prefix_pH)
    USprefix is str of US directory names prefix (prefix_US)
    USs is a dictionary of US and K to build 
    sn is str of sysname used to create filenames (prefix_sysname.extension)
    mdp is str with .mdp file path
    pbp is str with .pbp file path
    sites is str with .sites file path
    cphmd is str with CpHMD.sh file path
    gro is the name of the folder with all the 
    pH is str with pH value in directory name (prefix_pH)
    rep is the str with the number of the rep to be added to the folder 
  """
  if os.path.isdir(pHprefix + pH + "_r" + rep):
    raise Exception('Couldn\'t create {0}{1}. Directory already exists.'
                    .format(pHprefix, pH))
  os.system('mkdir -p {0}{1}_r{2}'.format(pHprefix,pH,rep))
  os.system('ln -rs REUS.py {0}{1}_r{2}'.format(pHprefix,pH,rep))
  os.system('ln -rs REUS_input.py {0}{1}_r{2}'.format(pHprefix,pH,rep))
  os.system('ln -rs srun_REUS.sh {0}{1}_r{2}'.format(pHprefix,pH,rep))
  #os.system('cp srun_REUS.dat {0}{1}_r{2}'.format(pHprefix,pH,rep))
  for US in USs.keys():
    os.system('mkdir {0}{1}_r{2}/{3}{4}'.format(pHprefix, pH, rep, USprefix, US))
##### REMOVER PARA SISTEMAS NORMAIS ######
    os.system('cp -d {}_{}.mdp {}{}_r{}/{}{}/{}_{}.mdp'
              .format(gro, US, pHprefix, pH, rep, USprefix, US, REUSprefix, sn))
    #os.system('cp -d {} {}{}/{}{}/{}_{}.mdp'
    #          .format(mdp, pHprefix, pH, USprefix, US, REUSprefix, sn))
    os.system('cp -d {} {}{}_r{}/{}{}/{}_{}.pbp'
              .format(pbp, pHprefix, pH, rep, USprefix, US, REUSprefix, sn))
    os.system('cp -d {} {}{}_r{}/{}{}/{}_{}.sites'
              .format(sites, pHprefix, pH, rep, USprefix, US, REUSprefix, sn))
    os.system('cp -d {} {}{}_r{}/{}{}/CpHMD.sh'
              .format(cphmd, pHprefix, pH, rep, USprefix, US))
    os.system('cp -d {} {}{}_r{}/{}{}/{}_{}.top'
              .format(top, pHprefix, pH, rep, USprefix, US, REUSprefix, sn))
    os.system('cp -d {} {}{}_r{}/{}{}/{}_{}_PDBin.pdb'
              .format(pdb, pHprefix, pH, rep, USprefix, US, REUSprefix, sn))
    os.system('cp -d {} {}{}_r{}/{}{}/{}_{}.ndx'
              .format(ndx, pHprefix, pH, rep, USprefix, US, REUSprefix, sn))
    
    
    ### Name specific files ####
    os.system('cp -d {}_{}.gro {}{}_r{}/{}{}/{}_{}_GROW.gro'
              .format(gro, US, pHprefix, pH, rep, USprefix, US, REUSprefix, sn))
    os.system('cp -d {}_{}.gro {}{}_r{}/{}{}/{}_000.gro'
              .format(gro, US, pHprefix, pH, rep, USprefix, US, sn))

def change_pHmdp(pHmdp, REUSprefix, sn, pHprefix, pH, rep, USprefix, US, effsteps, cphdir, sts, gmx, delphi, T):
  """
  Rewrites the pHmdp file in the pH folder to be consistent with the inputs
  provided in pHRE_input
  input:
    pHmdp is str of pHmdp file path
    REUSprefix is str used to create filenames (prefix_sysname.extension)
    sn is str of sysname used to create filenames (prefix_sysname.extension)
    prefix is str of directory names prefix (prefix_pH)
    pH is str with pH value in directory name (prefix_pH)
    effsteps is int of Effective MD Steps in CpHMD
    sts is str of St files path
    gmx is str of Gromacs path
    delphi is the location of the delphi tools folder
    T is the temperature of the system 
  """
  os.system('rm -rdf TMP')
  with open(pHmdp) as f, open('TMP', 'a') as f_new:
    for line in f:
      #Set the input gro name
      if 'export GROin=' in line:
        local_gro = '{0}_{1}_000.gro'.format(REUSprefix, sn)
        line      = 'export GROin={0}\n'.format(local_gro)
      #set the pH of the run
      elif 'export pH=' in line:
        line = 'export pH={0}\n'.format(pH)
      #set the number of effective steps
      elif 'export EffectiveSteps=' in line:
        line = 'export EffectiveSteps={0}\n'.format(effsteps)
      #set the CpH to use
      elif 'export CpHDIR=' in line:
        line = 'export CpHDIR={0}\n'.format(cphdir)
      #set the ST directory
      elif 'export StDIR=' in line:
        line = 'export StDIR={0}\n'.format(sts)
      #set the Delphi
      elif 'export DelphiDir=' in line:
        line = 'export DelphiDir={0}\n'.format(delphi)
      #set the groDir
      elif 'export GroDIR=' in line:
        line = 'export GroDIR={0}\n'.format(gmx)
      #set the temp
      elif 'export temp=' in line:
        line = 'export temp={0}\n'.format(T)

      #set the topology
      elif 'export TOPin=' in line:
        line = 'export TOPin={0}_{1}.top\n'.format(REUSprefix, sn)
      #set the pdb
      elif 'export PDBin=' in line:
        line = 'export PDBin={0}_{1}_PDBin.pdb\n'.format(REUSprefix, sn)
      #set the index
      elif 'export NDXin=' in line:
        line = 'export NDXin={0}_{1}.ndx\n'.format(REUSprefix, sn)
      #set the GRO whole
      elif 'export GROwhole=' in line:
        line = 'export GROwhole={0}_{1}_GROW.gro\n'.format(REUSprefix, sn)
      #set the topology
      #elif 'export nCPU=' in line:
      #  line = 'export nCPU={0}\n'.format(CPU)
        
      f_new.write(line)
    print("test")
  for US in USs.keys():
    print(US)
    os.system('cp TMP {0}{1}_r{2}/{3}{4}/{5}_{6}_001.pHmdp'
              .format(pHprefix, pH, rep, USprefix, US, REUSprefix, sn))


def change_mdp(mdp, REUSprefix, sn, pHprefix, pH, rep, USprefix, US, nstxtcout, dt, force, pullgroup1, pullgroup2, pullvec, init):
  """
  Rewrites the mdp file to be consistent with the inputs
  provided in pHRE_input
  input:
    mdp is str with .mdp file path
    REUSprefix is str used to create filenames (prefix_sysname.extension)
    sn is str of sysname used to create filenames (prefix_sysname.extension)
    prefix is str of directory names prefix (prefix_pH)
    pH is str with pH value in directory name (prefix_pH)
    nstxtcout is int with the value of nstxtcout to be changed in the .mdp file
    dt is int with the value of dt to be changed in the .mdp file
  """
  coord=int()

  os.system('rm -rdf TMP-mdp')
  with open(mdp) as f, open('TMP-mdp', 'a') as f_new:
    for line in f:
      if 'nstxout-compressed' in line:
        line = 'nstxout-compressed = {0}\n'.format(nstxtcout)
      elif 'dt' in line:
        line = 'dt = {0} ; ps!\n'.format(dt)
      #elif 'pull-group1-name' in line:
      #  line = 'pull-group1-name = {0}\n'.format(pullgroup1)
      #elif 'pull-group2-name' in line:
      #  line = 'pull-group2-name = {0}\n'.format(pullgroup2)
      #elif 'pull-coord1-vec' in line:
      #  line = 'pull-coord1-vec = {0}\n'.format(pullvec)
      elif 'pull-coord1-k' in line:
        line = 'pull-coord1-k = {0}\n'.format(force)
      elif 'pull-coord1-init' in line:
        line = 'pull-coord1-init = {0}\n'.format(init)
      f_new.write(line)
    
  os.rename('TMP-mdp', '{0}{1}_r{2}/{3}{4}/{5}_{6}.mdp'
            .format(pHprefix, pH, rep, USprefix, US, REUSprefix, sn))

#####################
# Main code
###

if __name__ == "__main__":
  for pH in pHs:
    for rep in reps:
      create_dirs(REUSprefix, pHprefix, USprefix, USs, sn, mdp, pbp, sites, cphmd, gro, pH, rep)
      change_pHmdp(pHmdp, REUSprefix, sn, pHprefix, pH, rep, USprefix, USs, effsteps, cphdir, sts, gmx, delphi, T)
      for US in USs.keys():
        force = USs[US]
        init = float(US)/10
        mdf = '{}{}_r{}/{}{}/{}_{}.mdp'.format(pHprefix, pH, rep, USprefix, US, REUSprefix, sn)
        change_mdp(mdf, REUSprefix, sn, pHprefix, pH, rep, USprefix, US, nstxtcout, dt, force, pullgroup1, pullgroup2, pullvec, init)
 
