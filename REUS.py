#! /usr/bin/python3

#####################
# Imports
###
import os
import sys
import math
import numpy as np
import random
import datetime
import subprocess
from multiprocessing import Process
from REUS_input import *

# Note: Using GMX5 the trjconv flags app and quiet do not exist

#####################
# Functions
###
def check_input(USprefix, REUSprefix, sts, gmx, pHmdp, mdp, pbp, sites, cphmd,
                effsteps, nstxout, dt, ncycles, REUSfreq, T, USs, argv, sn):
  """
  Checks if the user defined input arguments are as expected

  USprefix is str of US folders prefix name
  REUSprefix is str of CpH and RE template files
  sts is str of .st files directory path
  gmx is str of gromacs bin directory with trjconv and eneconv
  pHmdp is str of template pHmdp file (is only used by build_init_tree.py)
  mdp is str of template mdp file     (is only used by build_init_tree.py)
  pbp is str of template pbp file     (is only used by build_init_tree.py)
  sites is str of template sites file
  cphmd is str of CpHMD.sh template file
  effsteps is int of Effective MD Steps in CpHMD
  nstxout is int with the value of nstxout-compressed the .mdp file
  dt is float with the value of dt of the .mdp file
  ncycles is int of the number of cycles of REUS
  REUSfreq is int of the frequency of the RE step
  T is the temperature of the run (used in the exchange calculation function)
  Kb is the boltzmann constant (used in the exchange calculation function)
  USs is a dict of str of the simulated US values and force constants
  argv is a list with two elements
    in which the second element is the current block number
  sn is str of the system name
  """
  if len(argv) != 2:
    raise NameError('Please provide the number of the block to '
                    'perform REUS\n'
                    'Example of usage: python REUS.py 001')
  else:
    curr = str(argv[1])
    try:
      if int(curr) > 999 or int(curr) < 1:
        raise ValueError('curr should be a value between 0 and 999')
    except:
      raise ValueError('curr should be a value between 0 and 999')

  str_args = {'USprefix': USprefix, 'REUSprefix': REUSprefix, 'sts': sts,
              'gmx': gmx, 'pHmdp': pHmdp, 'mdp': mdp, 'pbp': pbp,
              'sites': sites, 'cphmd': cphmd}
  for arg in list(str_args.keys()):
    if type(str_args[arg]) is not str:
      raise TypeError('{0} should be defined in the '
                      'REUS_input.py file as a string'.format(arg))
    elif len(str_args[arg]) == 0:
      raise TypeError('{0} should be defined in the '
                      'REUS_input.py file as a string'.format(arg))

  int_args = {'effsteps': effsteps, 'nstxout': nstxout,
              'ncycles': ncycles, 'REUSfreq': REUSfreq}
  for arg in list(int_args.keys()):
    if type(int_args[arg]) is not int:
      raise TypeError('{0} should be defined in the '
                      'REUS_input.py file as an int'.format(arg))
  float_args = {'dt': dt, 'T': T}
  for arg in list(float_args.keys()):
    if type(float_args[arg]) is not float:
      raise TypeError('{0} should be defined in the '
                      'REUS_input.py file as an float'.format(arg))

  if type(USs) is not dict:
    raise TypeError('USs should be defined in the '
                    'REUS_input.py file as a dictionary'
                    'of ref positions and US force values')

  if ncycles < 1:
    raise ValueError('dt has to be bigger than 0')

  if REUSfreq < 1:
    raise ValueError('REUSfreq has to be bigger than 0')

  for US_value in USs.keys():
      US= '{0}{1}'.format(USprefix, US_value)
      if US not in os.listdir(os.getcwd()):
        raise ValueError('{0} is not a directory'.format(US))
      file_prefix = '{0}/{1}_{2}'.format(US, REUSprefix, sn)
      file_extensions = ('_001.pHmdp', '.mdp', '.pbp', '.sites')
      for file_extension in file_extensions:
        file_name = '{0}{1}'.format(file_prefix, file_extension)
        if not os.path.isfile(file_name):
          raise AttributeError('File {0} should exist\nYour REUSprefix is {1}'
                              .format(file_name, REUSprefix))
      file_name = '{0}/CpHMD.sh'.format(US)
      if not os.path.isfile(file_name):
        raise AttributeError('File {0} should exist\nYour REUSprefix is {1}'
                            .format(file_name, REUSprefix))

      with open('{0}/{1}_{2}_001.pHmdp'
                .format(US, REUSprefix, sn)) as f:
        for line in f:
          line = line.strip()
          if 'export GROin=' in line:
            local_gro = '{0}_{1}_000.gro'.format(REUSprefix, sn).strip()
            if local_gro != line.split('=')[-1]:
              raise ValueError('Inconsistent gro name in .pHmdp file in {0}'
                              .format(US))
          elif 'export EffectiveSteps=' in line:
            local_effsteps = line.split('=')[-1].replace(' ', '').strip()
            if local_effsteps != str(effsteps):
              raise ValueError('Inconsistent effsteps values in '
                              '.pHmdp and REUS_input files in {0}'.format(US))
          elif 'export StDIR=' in line:
            local_sts = line.split('=')[-1].replace(' ', '').strip()
            if local_sts != sts:
              raise ValueError('Inconsistent sts values in '
                              '.pHmdp and REUS_input files in {0}'.format(US))
          elif 'export GroDIR=' in line:
            local_gmx = line.split('=')[-1].replace(' ', '').strip()
            if local_gmx != gmx:
              raise ValueError('Inconsistent gmx values in '
                              '.pHmdp and REUS_input files in {0}'.format(US))
          elif 'export InitCycle=' in line:
            initcycle = line.split('=')[1].replace(' ', '')
            if initcycle != '1':
              print('Warning: InitCycle in .pHmdp not 1')

          elif 'export EndCycle=' in line:
            endcycle = line.split('=')[1].replace(' ', '')
            if endcycle != '1':
              print('Warning: EndCycle in .pHmdp not 1')

      ### Add check for correct umbrella
      USv=int(US_value)/10

      with open('{0}/{1}_{2}.mdp'
                .format(US, REUSprefix, sn)) as f:
        for line in f:
          if ';' in line:
            line = line.split(';')[0]
          if 'nstxout-compressed' in line:
            local_nstxout = line.split(';')[0].split('=')[-1].replace(' ', '').strip()
            if local_nstxout != str(nstxout):
              raise ValueError('Inconsistent nstxout-compressed values in '
                              '.mdp and REUS_input files in {0}'.format(US))
          elif 'dt' in line:
            local_dt = line.split(';')[0].split('=')[-1].replace(' ', '').strip()
            if local_dt != str(dt):
              raise ValueError('Inconsistent dt values in '
                              '.mdp and REUS_input files in {0}'.format(US))
          elif 'pull-coord1-init' in line:
            local_init = line.split(';')[0].split('=')[-1].replace(' ', '').strip()
            if local_init != str(USv):
              raise ValueError('Inconsistent reference value ({0}) in '
                              '.mdp and REUS_input files in {1}'.format(local_init,US))                   

  dir_list = (sts, gmx)
  for dir_i in dir_list:
    if dir_i[-1] != '/':
      dir_i = dir_i + "/"
    if not os.path.isdir(sts):
      raise ValueError('sts should be the path of a non-empty directory')


  return curr


def remove_oldfiles(curr, sn, USs, USprefix, REUSprefix):
  """
  Removes the old GLOB files by appending .old to their names:
    GLOBLOG_XXX is log file for the pHRE
    GLOBOUT_XXX is the standard output file
    GLOBERR is the error file
  Removes the old sn_curr.* files in each pH folder
  input:
    curr is str of current cycle
    REUSprefix is str of CpH and RE template files
    sn is str of the system name
    USs is a dictionary of str of the simulated US values
    and their K values.
    USprefix is str of US folders USprefix name
  output:
    {'LOG': GLOBLOG, 'OUT': GLOBOUT, 'ERR': 'GLOBERR'}
    GLOBfiles_dic is dict of the GLOB files names
  """
  GLOBLOG = 'GLOBLOG_' + curr
  GLOBOUT = 'GLOBOUT_' + curr
  GLOBDEB = 'GLOBDEB_' + curr
  GLOBfiles_dict = {'LOG': GLOBLOG, 'OUT': GLOBOUT, 'ERR': 'GLOBERR', 'DEB': GLOBDEB}
  for GLOBfile in list(GLOBfiles_dict.values()):
    if os.path.isfile(GLOBfile):

        subprocess.run("cat {0} >> {0}.old".format(GLOBfile))
        subprocess.run("/bin/rm {0}".format(GLOBfile))

  extensions_list = ('edr', 'gro', 'mocc', 'occ', 'xtc', 'tpr')
  for US in USs.keys():
    for extension in extensions_list:
      f = '{0}{1}/{2}_{3}_{4}.{5}'.format(USprefix, US, REUSprefix, sn, curr, extension)
      if os.path.isfile(f):
        subprocess.run("mv {0} {0}.old".format(f))

  return GLOBfiles_dict

def read_sites(sites):
  """
  Reads number of residues in the sites file
  returns the type of simulation being MD if sites
  is empty (Unnused still)
  """
  f=os.stat(sites).st_size

  if f==0:
    simtype="MD"
  else:
    simtype="CpHMD"  
   
  return simtype


def gen_vars_n_copy_gros(dt, nstxout, effsteps, cphmd, ncycles,
                          USs, USprefix, sn, curr, REUSprefix):
  """
  Generate useful variables
  Update the .gro file to be used in this REUS segment
  The gro file used in the REUS segment is called REUSprefix_sn_000.gro
  """
  prev = '{0:03d}'.format(int(curr) - 1)
  for US in USs.keys():
    cp_gro = "cp {0}{1}/{2}_{3}.gro {0}{1}/{4}_{2}_000.gro".format(USprefix,
              US, sn, prev, REUSprefix)
    subprocess.run(cp_gro, shell=True)
  # These variables will be used by trjconv
  # to split trajectories during exchanges
  b1 = dt * nstxout
  #e1 = dt * effsteps / 2
  #b2 = dt * effsteps / 2 + nstxout * dt
  e2 = dt * effsteps

  # Usefull list to log exchanges and template of the lists
  USs_tmp = [i for i in USs.keys()]
  
  USs_exc = [i for i in USs.keys()]

  # Initial time of current block (ps)
  #initt      = int(prev) * effsteps * dt * ncycles + nstxout * dt #This one may only be for the pHRE
  initt      = int(prev) * effsteps * dt * ncycles
  # There are two types of replica exchanges
  # (TYPE 1: 1-2,3-4 OR TYPE 2: 2-3,4-5)
  pHexcTYPE = 1

  # CpHMD.sh file name
  cphmdsh    = cphmd.split('/')[-1].replace(' ', '')

  # Current Path
  pwd = os.getcwd()

  return b1, e2, USs_exc, USs_tmp, pwd, cphmdsh, initt, pHexcTYPE


def CpHMD_housekeeping(path, runname):
  """
  Add the content of the .log file into a general LOGFILE
  Add the content of the .info file into a general INFOFILE
  Add the content of the ERROR file into a general ERRORFILE
  Remove initial gro in order to avoid repeated (and wrong) CpHMD segments
  in the next steps
  """
  # file_list -> (old_file, logfile)
  # file_list[0] -> old_file file_list[1] -> logfile
  file_list = (('{}_001.log'.format(runname), 'LOGFILE'),
               ('{}_001.info'.format(runname), 'INFOFILE'),
               ('ERROR', 'ERRORFILE'))

  for file_i in file_list:
    fill = 'cd {0} && cat {1} >> {2}'.format(path, file_i[0], file_i[1])
    subprocess.run(fill, shell=True)
    rm = 'cd {0} && /bin/rm -f {1}'.format(path, file_i[0])
    subprocess.run(rm, shell=True)

  rm_gro = 'cd {} && /bin/rm -f {}_000.gro >> INFOFILE'.format(path, runname)
  subprocess.run(rm_gro, shell=True)


def CpHMD_seg(path, runname, cphmdsh):
  """
  Runs a single CpHMD segment (PB->MC->MMsol->MDwhole)
  And does some housekeeping
  input: 
    path is str of the directory where CpHMD is being executed
    runname is str with the prefix of the CpHMD files
    cphmdsh is str of the name of the CpHMD.sh script
  output:
    LOGFILE, ERRORFILE and INFOFILE are updated
    runname_001.{edr,gro,mocc,occ,xtc,tpr} are created
  """
  # Run CpHMD segment
  run_CpHMD = 'cd {0} && nice -n 19 ./{1} {2}_001.pHmdp >> LOGFILE 2>> ERRORFILE'\
              .format(path, cphmdsh, runname)
  subprocess.run(run_CpHMD, shell=True)

  # Remove some unimportant files
  useless_files_list = ['TMP_delphi.crg', 'TMP_delphi.out',
                        'TMP_*.gro', 'TMP_protein.gro', 'aux.pbp',
                        runname + '_mod.sites', runname + '_stmod.sites',
                        runname + '_000.gro']
  for f in useless_files_list:
    rm_f = 'cd {0} && /bin/rm -f {1}'.format(path, f)
    subprocess.run(rm_f, shell=True)

  # Rename pbp file back to its original name
  rename_pbp = 'cd {0} && mv DELPHI.pbp {1}'.format(path, runname + '.pbp')
  subprocess.run(rename_pbp, shell=True)

  # Report to LOGFILE, INFOFILE and ERRORFILE
  CpHMD_housekeeping(path, runname)

  return None


def run_CpHMD_seg(pwd, REUSprefix, sn, cphmdsh, USs, USprefix):
  """
  Starts one thread of CpHMD per pH values in pHs
  Waits for all the threads to end
  """
  threads = []
  
  for US in USs.keys():
    path = pwd + '/' + USprefix + US
    threads.append(Process(target=CpHMD_seg,
                          args=(path, REUSprefix + '_' + sn, cphmdsh)))
  for t in threads:
    t.start()
  for t in threads:
    t.join()


def log_in_GLOB(GLOBLOG, initcurr, effsteps, dt, ncycles, REUSfreq, flag):
  """
  Reports to GLOBLOG
  input:
    flag is str which can be 'REUS' when a REUS is going to start
                             'CpHMD' when a CpHMD step is going to start
                             'first' when the first block is going to start
  """
  with open(GLOBLOG, 'a') as fo:
    if flag == 'REUS':
      log = 'REUS step at t = {0} ps :: Starting at {1}\n'\
            .format(initcurr + effsteps * dt,
                    datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S"))
    elif flag == 'CpHMD':
      log = 'CpHMD step (t = {0} ps) :: Starting at {1}\n'\
            .format(initcurr,
                    datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S"))
    elif flag == 'first':
      log = 'REUS cycle lasts {0} ps\n'\
            'RE steps are made every {1} ps\n'\
            .format(effsteps * ncycles * dt, effsteps * REUSfreq * dt)

    else:
      raise NameError('Variable flag not defined. '
                      'The programmer is be to blamed')
    
    fo.write(log)


def get_positions(USs, USprefix, REUSprefix, sn):
  """
  Reads the pullx file (the file of the CpHMD-RE) and returns
  the number of position of each simulated US window

  input:
  #   pHprefix is str of the pH directory names prefix
  #   pHs is a list of str of the simulated pH values
  #   USprefix is str of the US directory names prefix
  #   USs is a list of str of the simulated US values
  #   REUSprefix is str of the pHRE files prefix
  #   sn is str of the system name
  #   crg_ndx is list of ints with the number of tautomers of a site
  #   natCA is list of str of the nature of the titrable species
  #         the order of the list is given by the appearence in the .sites file
  ouput:
    pos is a list of the last position of the reference molcule
    the order is the same as in the USs dictionary
  """

  pos=[]
  for US in USs.keys():
    US_pos = read_pullx('{0}{1}/{2}_001_pullx.xvg'.format(USprefix, US,REUSprefix + '_' + sn))
    pos.append(US_pos)

  return pos


def read_pullx(pullx_file):
  """
  Reads the current pullx file to extract the last line
  Returns the position of the US.
  """
  with open(pullx_file) as fi:
    lastline=fi.readlines()[-1]
    a=lastline.split()
    lastpos=a[1]

  return lastpos


def run_REUS(GLOBLOG, GLOBOUT, pHexcTYPE, USs, USs_exc, USs_tmp, USprefix, REUSprefix,
             sn, curr, gmx, b1, e2, initt, initcurr):
  """
  Wrapper function of all the REUS related functions and reports to GLOBLOG
  """

  pos = get_positions(USs, USprefix, REUSprefix, sn)
  with open(GLOBLOG, 'a') as fo:
    fo.write('REUS step of TYPE {0}.\n'.format(pHexcTYPE))

    # Tests all possible exchanges for a given type
    # exc_list = (Accepted pHs, Tested_pHs, Probability of Accepted pHs)
    exc_list = REUS_step(USs, pos, pHexcTYPE)
    for i in range(len(exc_list[1])):
      fo.write('P({0}) = {1}.\n'.format(exc_list[1][i], exc_list[2][i]))


      #fo.write('Positions {0}.\n'.format(exc_list[3][i]))

    # Update of the xtc, edr, occ and mocc of both pH values to be exchanged
    for n, exc in enumerate(exc_list[0]):
      fo.write('Positions: {0} .\n'.format(exc_list[3][n]))
      fo.write('{0} exchange was accepted.\n'.format(exc))
      
      perf_exchange(exc, USprefix, REUSprefix + '_' + sn,
                    sn + '_' + curr, gmx, b1, e2, initt, initcurr,GLOBOUT)

  # Report exchanges
  with open(GLOBOUT, 'a') as fo:
    fo.write('{0} '.format(initcurr + effsteps * dt))
    for n, exc in enumerate(exc_list[0]):
      US1, US2 = exc.split('_')
      US1i, US2i = USs_tmp.index(US1), USs_tmp.index(US2)
      USs_exc[US2i], USs_exc[US1i] = USs_exc[US1i], USs_exc[US2i]
    
    for US in USs_exc:
      fo.write('{0} '.format(US))
    
    fo.write('\n')
  # Change the type of exchanges for the next step
  if (pHexcTYPE == 1):
    pHexcTYPE = 2
  else:
    pHexcTYPE = 1

  return pHexcTYPE


def REUS_step(USs, Pos, typ):
  """
  Performs a full USRE step
  i.e. tries all possible exchanges and returns a list with the accepted ones
  It can try two types of exchanges (TYPE 1: 1-2,3-4 OR TYPE 2: 2-3,4-5)
  input:
    USs is a dictionary of str of the simulated US window values
    Pos is a list of the positions of the reference molecule
            the order is the same as in USs
    typ is int representing the type of the exchange. Belongs to {1, 2}
  output:
    (excs, excs_all, probs_all)
    excs is list of str of accepted US exchanges value
    excs_all is list of str of tested US exchanges
    probs_all is list of float of probabilities of exchange excs_all[i]
    umbposi is a list of the current positions of the US windows
  """
  excs      = []
  excs_all  = []
  probs_all = []
  umbposi   = []
  if (len(Pos) != len(USs.keys())):
    raise ValueError('Number of US windows is different to the number of '
                     'entries in positions list.')
  for t in range(typ - 1, len(USs.keys()) - 1, 2):
    Km     = list(USs.values())[t]
    refm   = list(USs.keys())[t]
    posi   = Pos[t]

    Kn     = list(USs.values())[t + 1]
    refn   = list(USs.keys())[t + 1]
    posj   = Pos[t + 1]
  

    exc = REUS_exc(Km, refm, posi, Kn, refn, posj, T)
    excs_all.append(refm + '_' + refn)
    probs_all.append(exc[1])
    
    if exc[0]:
      excs.append(list(USs.keys())[t] + '_' + list(USs.keys())[t + 1])
      umbposi.append(posi + '_' + posj)

  return excs, excs_all, probs_all, umbposi


def REUS_exc(Km, refm, posi, Kn, refn, posj, T):
  """
  Decides if the US exchange is accepted
  Returns a boolean
  function = exp( -ln(10) * (pH1 - pH2) * (nprots2 - nprots1))
  input:
    pH1 is str of pH value
    nprots1 is int of the number of protons in pH1
    pH2 is str of pH value
    nprots2 is int of the number of protons in pH2
  output:
    (boolean, prob)
    boolean -> True if exchange is accepted
            -> False if exchange is rejected
    prob is the probability of exchange between pH1 and pH2
  """

  ### Kb usado 0.001985875 Kcal /(mol K)
  ##### Valores necessários Kb, T, Km (constante força umb m), Kn (constante força umb n), 
  ##### refm (posição referencia da umb m), refn (posição referencia da umb n), 
  ##### posj (posição final da umbrella m obtida do pullx), posi (posição final da umbrella n obtida do pullx)

  Umi= 0.5 * Km * ((float(posi) - float(refm)/10)**2)
  Umj= 0.5 * Km * ((float(posj) - float(refm)/10)**2)

  Unj= 0.5 * Kn * ((float(posj) - float(refn)/10)**2)
  Uni= 0.5 * Kn * ((float(posi) - float(refn)/10)**2)


####### (d Um)-(d Un) = (Umj - Umi) - (Unj - Uni)
  function = np.exp(-(1/(0.001985875*T)) * (Umj - Umi - Unj + Uni) )

  """ function = np.exp(-(1/(0.001985875*T)) * (
                     0.5 * Km * (float(posj) - float(refm)/10)*(float(posj) - float(refm)/10) - #Umj
                     0.5 * Km * (float(posi) - float(refm)/10)*(float(posi) - float(refm)/10) - #Umi
                     0.5 * Kn * (float(posj) - float(refn)/10)*(float(posj) - float(refn)/10) + #Unj
                     0.5 * Kn * (float(posi) - float(refn)/10)*(float(posi) - float(refn)/10)   #Uni
                     )) """

  #### Debug variables ####
  subprocess.run("echo {0} {1} {2} {3} {4} {5} {6} >> GLOBDEB".format(refm, Km, posj, refn, Kn, posi, T), 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)
  subprocess.run("echo Umbs {0} {1} {2} {3} >> GLOBDEB".format(Umj, Umi, Unj, Uni), 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)

  prob = min(1, function)
  if (prob == 1):
    return True, prob
  else:
    ref = random.random()
    if (prob > ref):
      return True, prob
    else:
      return False, prob


def perf_exchange(exc, USprefix, REUSrunname, runname,
                  gmx, b1, e2, initt, initcurr, GLOBOUT):
  """
  Performs an exchange between two pH values
  REUSrunname_001.gro files are moved to the other directory
  1st half of edr and xtc files are added to the runname files
  mocc and occ files are added to both pH values
  input:
    exc is str of two pH values to exchange separated by a underscore
    prefix is str of pH folders prefix name
    REUSrunname is str of pHREprefix_sn ex: REUS_penta
    runname is str of systemname and the current frame ex:penta_001
    gmx is str of gromacs bin directory with trjconv and eneconv
    b1 is int of the time of first frame of the CpH step
    e2 is int of the last frame of the CpH step
    initt is int of the initial time of current block (ps)
    initcurr is int of initial time of current CpHRE cycle (ps)
  """
  US1, US2 = exc.split('_')

  US_list = ((US1,US2),(US2,US1))
  for i in US_list:
    US_i  = i[0]
    US_ii = i[1]

    if US_i == US_ii:
      raise NameError('The exchange US are the same')

    # exchange gro files from US1 to US2 and vice-versa
    subprocess.run("mv {0}{1}/{2}_001.gro {0}{3}/{2}_000.gro"
                .format(USprefix, US_i, REUSrunname, US_ii), 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)

  ################################################################################
  ################################################################################
  ######### Deal with mdp references of the membrane being different #############
  ################################################################################
  ################################################################################

  #Get the reference 
  with open('{0}{1}/{2}.mdp'.format(USprefix, US1, REUSrunname)) as f:
      for line in f:
        if 'pull-group2-name' in line:
          local_refi = line.split(';')[0].split('=')[-1].replace(' ', '').strip()
    
  with open('{0}{1}/{2}.mdp'.format(USprefix, US2, REUSrunname)) as f:
    for line in f:
        if 'pull-group2-name' in line:
          local_refii = line.split(';')[0].split('=')[-1].replace(' ', '').strip()


  if local_refi != local_refii:
    #### make mdp tmp on US2 based on US1 ###
    with open('{0}{1}/{2}.mdp'.format(USprefix, US2, REUSrunname)) as f, open('{0}{1}/{2}_TMP.mdp'.format(USprefix, US2, REUSrunname), 'a') as f_new:
      for line in f:
        if 'pull-group2-name' in line:
          line = 'pull-group2-name = {0}\n'.format(local_refi)
        f_new.write(line)
    #### make mdp tmp on US1 based on US2 ###
    with open('{0}{1}/{2}.mdp'.format(USprefix, US1, REUSrunname)) as f, open('{0}{1}/{2}_TMP.mdp'.format(USprefix, US1, REUSrunname), 'a') as f_new:
      for line in f:
        if 'pull-group2-name' in line:
          line = 'pull-group2-name = {0}\n'.format(local_refii)
        f_new.write(line)
    ### mv tmp mdp to the base mdp in each file ###
    subprocess.run('mv {0}{1}/{2}.mdp {0}{1}/{2}.mdp.bk'.format(USprefix, US1, REUSrunname), 
                shell=True)
    subprocess.run('mv {0}{1}/{2}_TMP.mdp {0}{1}/{2}.mdp'.format(USprefix, US1, REUSrunname), 
                shell=True)

    subprocess.run('mv {0}{1}/{2}.mdp {0}{1}/{2}.mdp.bk'.format(USprefix, US2, REUSrunname), 
                shell=True)
    subprocess.run('mv {0}{1}/{2}_TMP.mdp {0}{1}/{2}.mdp'.format(USprefix, US2, REUSrunname), 
                shell=True)
    ##Report changed on glob
    with open(GLOBOUT, 'a') as fo:
      fo.write('mdp was exchanged between {0} and {1}\n'.format(US1, US2))
      
###################################################################################
###################################################################################
###################################################################################


  for i in (US1, US2):
    # if runname.tpr does not exist, create it from the REUS.tpr
    # else remove the REUS.tpr
    if not os.path.isfile("{0}{1}/{2}.tpr".format(USprefix, i, runname)):
      subprocess.run('mv {0}{1}/{2}_001.tpr {0}{1}/{3}.tpr'.format(USprefix, i,
                                                              REUSrunname,
                                                              runname), shell=True)
    else:
      subprocess.run('/bin/rm {0}{1}/{2}_001.tpr'.format(USprefix, i, REUSrunname), shell=True)

    # Cat with the mocc and occ files on both US being exchanged
    for extension in ('mocc', 'occ'):
      subprocess.run("cat {0}{1}/{2}_001.{4} >> {0}{1}/{3}.{4}"
                  .format(USprefix, i, REUSrunname, runname, extension), shell=True)
      subprocess.run("/bin/rm -f {0}{1}/{2}_001.{3}"
                .format(USprefix, i, REUSrunname, extension), shell=True)

    # update the runname.xtc file and delete the pHRE xtc file
    if os.path.isfile('{0}{1}/{2}.xtc'.format(USprefix, i, runname)):
      b = b1
      subprocess.run("{0} trjconv -f {1}{2}/{3}_001.xtc -o {1}{2}/{4}_tmp.xtc \
             -t0 {5} -b {6}>> GLOBERR 2>> GLOBERR".format(gmx, USprefix, i, REUSrunname, runname, initcurr,b), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)
 
      subprocess.run("{0} trjcat -f {1}{2}/{3}_tmp.xtc {1}{2}/{3}.xtc -o {1}{2}/{3}.xtc \
             >> GLOBERR 2>> GLOBERR".format(gmx, USprefix, i, runname), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)
     
    else:
      b = b1
      subprocess.run("{0} trjconv -f {1}{2}/{3}_001.xtc -o {1}{2}/{4}.xtc \
             -t0 {5} -b {6}>> GLOBERR 2>> GLOBERR".format(gmx, USprefix, i, REUSrunname, runname, initcurr,b), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)
    
    subprocess.run("/bin/rm -f {0}{1}/{2}_001.xtc".format(USprefix, i, REUSrunname), shell=True)
    subprocess.run("/bin/rm -f {0}{1}/{2}_tmp.xtc".format(USprefix, i, runname), shell=True)

    # update the runname.edr file and delete the REUS edr file
    subprocess.run('{0} eneconv -f {1}{2}/{3}_001.edr -o {1}{2}_TMP1.edr \
               >> GLOBERR 2>> GLOBERR'
              .format(gmx, USprefix, i, REUSrunname), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)

    if os.path.isfile('{0}{1}/{2}.edr'.format(USprefix, i, runname)):
      subprocess.run('echo "{0}\n{1}" | {2} eneconv -f {3}{4}/{5}.edr \
                {3}{4}_TMP1.edr -o {3}{4}_TMP3.edr -settime \
                >> GLOBERR 2>> GLOBERR'.format(initt, initcurr, gmx,
                                               USprefix, i, runname), shell=True)
                                      
    else:
      subprocess.run("mv {0}{1}_TMP1.edr {0}{1}_TMP3.edr".format(USprefix, i), shell=True)
  
    subprocess.run("mv {0}{1}_TMP3.edr {0}{1}/{2}.edr".format(USprefix, i, runname), shell=True)
    subprocess.run("/bin/rm -f {0}{1}_TMP?.edr".format(USprefix, i), shell=True)
    subprocess.run("/bin/rm -f {0}{1}/{2}_001.edr".format(USprefix, i, REUSrunname), shell=True)
  
  
    ##### Insert a way to concatenate correctly the pullx and pullf

    for extens in ('pullx', 'pullf'):
      subprocess.run("awk -v curr={0} '!/^#/ && !/^@/ {{print $1+curr ,$2}}' {1}{2}/{3}_001_{5}.xvg >>  {1}{2}/{4}_{5}.xvg".format(initcurr, USprefix, i, REUSrunname, runname, extens), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)
    
      subprocess.run("/bin/rm -f {0}{1}/{2}_001_{3}.xvg".format(USprefix, i, REUSrunname, extens), shell=True)


def run_noexchange(USs, USprefix, REUSrunname, runname,
                   gmx, initt, initcurr, GLOBLOG, b1):
  """
  input:
    pHs is a list of str of the simulated pH values
    prefix is str of pH folders prefix name
    pHRErunname is str of pHREprefix_sn ex: pHRE_cyclic_PEI
    runname is str of systemname ex: cyclic_PEI
    gmx is str of gromacs bin directory with trjconv and eneconv
    initt is int of the initial time of current block (ps)
    initcurr is int of initial time of current CpHRE cycle (ps)
    GLOBLOG is a str of the log file name for the pHRE
  """
  for US in USs.keys():
    if os.path.isfile('{0}{1}/{2}_000.gro'.format(USprefix, US, REUSrunname)):
      with open(GLOBLOG, 'a') as fo:
        fo.write('Files from US {0} were already treated in an exchange.\n'
                 .format(US))
    else:
      perf_noexchange(US, USprefix, REUSrunname, runname, gmx, initt,
                      initcurr, GLOBLOG, b1)
  subprocess.run("echo >> {0}".format(GLOBLOG), shell=True)


def perf_noexchange(US, USprefix, REUSrunname, runname, gmx, initt,
                    initcurr, GLOBLOG, b1):
  """
  Concatenates edr, xtc, occ, and mocc files
  Note: xtc is actually concatenated, in edr file the times of the previous
        frames are being changed as well
  REUSrunname_001.gro is moved to REUSrunname_000.gro
  Input:
    pHs is a list of str of the simulated pH values
    prefix is str of pH folders prefix name
    pHRErunname is str of pHREprefix_sn ex: pHRE_cyclic_PEI
    runname is str of systemname ex: cyclic_PEI
    gmx is str of gromacs bin directory with trjconv and eneconv
    initt is int of the initial time of current block (ps)
    initcurr is int of initial time of current CpHRE cycle (ps)
    GLOBLOG is a str of the log file name for the pHRE
  """
  # update the gro file to be used in the REUS cycle
  subprocess.run("mv {0}{1}/{2}_001.gro {0}{1}/{2}_000.gro".format(USprefix, US, REUSrunname), shell=True)
  # create tpr file if it does not exist
  if not os.path.isfile('{0}{1}/{2}.tpr'.format(USprefix, US, runname)):
    subprocess.run("mv {0}{1}/{2}_001.tpr {0}{1}/{3}.tpr".format(USprefix, US,
                                                            REUSrunname,
                                                            runname), shell=True)
  else:
    subprocess.run("/bin/rm {0}{1}/{2}_001.tpr".format(USprefix, US, REUSrunname), shell=True)

# update the runname.occ file and delete the pHRE occ file
  for extension in ('mocc', 'occ'):
    subprocess.run("cat {0}{1}/{2}_001.{4} >> {0}{1}/{3}.{4}"
                .format(USprefix, US, REUSrunname, runname, extension), shell=True)
    subprocess.run("/bin/rm -f {0}{1}/{2}_001.{3}"
              .format(USprefix, US, REUSrunname, extension), shell=True)

  # update the runname.xtc file and delete the pHRE xtc file
  if os.path.isfile('{0}{1}/{2}.xtc'.format(USprefix, US, runname)):
    b = b1

    subprocess.run("{0} trjconv -f {1}{2}/{3}_001.xtc -o {1}{2}/{4}_tmp.xtc \
           -t0 {5} -b {6}>> GLOBERR 2>> GLOBERR".format(gmx, USprefix, US, REUSrunname, runname, initcurr,b), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)
 
    subprocess.run("{0} trjcat -f {1}{2}/{3}_tmp.xtc {1}{2}/{3}.xtc -o {1}{2}/{3}.xtc \
           >> GLOBERR 2>> GLOBERR".format(gmx, USprefix, US, runname), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)
     
  else:
    b = b1
    
    subprocess.run("{0} trjconv -f {1}{2}/{3}_001.xtc -o {1}{2}/{4}.xtc \
           -t0 {5} -b {6}>> GLOBERR 2>> GLOBERR".format(gmx, USprefix, US, REUSrunname, runname, initcurr,b), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)

  #subprocess.run("/bin/rm -f {0}{1}/{2}_001.xtc".format(USprefix, US, REUSrunname), shell=True)
  subprocess.run("mv {0}{1}/{2}_001.xtc {0}{1}/{2}_{3}_OLD.xtc".format(USprefix, US, REUSrunname, initcurr), shell=True)
  subprocess.run("/bin/rm -f {0}{1}/{2}_tmp.xtc".format(USprefix, US, runname), shell=True)

  # update the runname.edr file and delete the pHRE edr file
  subprocess.run("{0} eneconv -f {1}{2}/{3}_001.edr \
          -o {1}{2}_TMP1.edr >> GLOBERR 2>> GLOBERR".format(gmx, USprefix, US, REUSrunname), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)

  if os.path.isfile("{0}{1}/{2}.edr".format(USprefix, US, runname)):
    subprocess.run('echo "{0}\n{1}" | {2} eneconv -f {3}{4}/{5}.edr \
              {3}{4}_TMP1.edr -o {3}{4}_TMP3.edr -settime \
              >> GLOBERR 2>> GLOBERR'.format(initt, initcurr, gmx,
                                             USprefix, US, runname), shell=True)
                                      
  else:
    subprocess.run("mv {0}{1}_TMP1.edr {0}{1}_TMP3.edr".format(USprefix, US), shell=True)
  
  subprocess.run("mv {1}{0}_TMP3.edr {1}{0}/{2}.edr".format(US, USprefix, runname), shell=True)
  subprocess.run("/bin/rm -f {0}{1}_TMP?.edr".format(USprefix, US), shell=True)
  subprocess.run("/bin/rm -f {0}{1}/{2}_001.edr".format(USprefix, US, REUSrunname), shell=True)

  ##### Insert a way to concatenate correctly the pullx and pullf
  for extens in ('pullx', 'pullf'):  
    subprocess.run("awk -v curr={0} '!/^#/ && !/^@/ {{print $1+curr ,$2}}' {1}{2}/{3}_001_{5}.xvg >>  {1}{2}/{4}_{5}.xvg".format(initcurr, USprefix, US, REUSrunname, runname, extens), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)

def rename_gros(USprefix, USs, REUSrunname, sn, curr):
  """
  Rename .gro files from REUS nomemclature to their current block
           REUSrunname_000.gro -> sn_curr.gro
  Example: 'REUStest_000.gro'  -> 'test_069.gro'
  """
  for US in USs.keys():
    subprocess.run('cp {0}{1}/{2}_000.gro {0}{1}/{3}_{4}.gro'
              .format(USprefix, US, REUSrunname, sn, curr), shell=True)

if __name__ == "__main__":
  # (re)define sites variable since the template might not be available
  sites = '{0}{1}/{2}_{3}.sites'.format(USprefix,list(USs.keys())[0],REUSprefix,sn)
  # Check input files
  curr = check_input(USprefix, REUSprefix, sts, gmx, pHmdp, mdp, pbp,
                     sites, cphmd, effsteps, nstxtcout, dt,
                     ncycles, REUSfreq, T, USs, sys.argv, sn)
  ### Run some preparation steps ###
  # Remove old GLOB files
  GLOB = remove_oldfiles(curr, sn, USs, USprefix, REUSprefix)

  # Generate some useful variables and give correct name to gro file
  b1, e2, USs_exc, USs_tmp, pwd,\
  cphmdsh, initt, pHexcTYPE  = gen_vars_n_copy_gros(dt, nstxtcout, effsteps,
                                                    cphmd, ncycles, USs, USprefix,
                                                    sn, curr, REUSprefix)
  log_in_GLOB(GLOB['LOG'], '', effsteps, dt, ncycles, REUSfreq, 'first')


  ###   Main cycle   ###
  for step in range(ncycles):
    # Initial time of current CpHRE cycle (ps)
    initcurr = initt + step * effsteps * dt
    #if initcurr == 20:
    #  indicador = subprocess.run("sid | egrep 'test_pHRE' | awk '{print $1}'")
    #  subprocess.run("scancel {0}".format(indicador), shell=True)
    ### Run a CpHMD segment at each pH value ###
    log_in_GLOB(GLOB['LOG'], initcurr, effsteps, dt,
                ncycles, REUSfreq, 'CpHMD')
    run_CpHMD_seg(pwd, REUSprefix, sn, cphmdsh, USs, USprefix)

    ### Run a pHRE step ###
    if ((step + 1) % REUSfreq == 0):
      log_in_GLOB(GLOB['LOG'], initcurr, effsteps, dt,
                  ncycles, REUSfreq, 'REUS')

      pHexcTYPE = run_REUS(GLOB['LOG'], GLOB['OUT'], pHexcTYPE, USs, USs_exc, USs_tmp, USprefix, REUSprefix, sn, curr, gmx, b1, e2,
                           initt, initcurr)

    ### Run a non exchange step in each pH dir ###
    run_noexchange(USs, USprefix, REUSprefix + '_' + sn, sn + '_' + curr,
                   gmx, initt, initcurr, GLOB['LOG'], b1)

  ### Rename last gro files ###
  rename_gros(USprefix, USs, REUSprefix + '_' + sn, sn, curr)
