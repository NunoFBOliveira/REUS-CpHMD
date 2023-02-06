#!/bin/bash -e
#
#
# Define runtime variables
awk '1;/REUS Inputs/{exit}' REUS_input.py > tmp.dat
source tmp.dat
prefix=$USprefix
sysname=$sn
runname="${REUSprefix}_${sysname}"
rundir=`pwd -P`
runfile=`basename $0`
#
# cleaning old queue files
rm -f $runname.slurm
#

USs=`sed -n '/{/,/}/p' REUS_input.py | awk '$1!="USs" && $1!="}" {print $1}' | awk -F "'" '{printf("%4s", $2)}'`

###Check if requested cpus is not higher than the used ###
## count number of entries in the US window dictionary

nw=`sed -n '/{/,/}/p' REUS_input.py | grep -o : | wc -l`

ucpus=`echo $ncpus $nw | awk '{print int($1/$2)}'`

### Compute the number of cpus to use in the DelphiT ###
# first check if ucpus is not already an available number on the delphi list #
if [[ $ucpus == @(1|2|4|6|8|10|12|16|24|32|64) ]]
then
    dcpus=$ucpus
else

    dcpus=`echo "1|2|4|6|8|10|12|16|24|32|64" | awk -v val=$ucpus -v RS="|" '{D=(val - $1)*(val - $1);if((VAL<=$1)&&(((!SET) || (D < DIFF)))){DIFF=D;X=$1;SET=1;}}END{if(!SET) print "No value found"; else print X }' `
  
fi

eff=`echo $nw $ncpus $ucpus | awk '{printf("%.0f\n", (($1*$3)/$2)*100)}'`

if [ "$eff" -lt 80 ]
then
    echo "REUS efficiency:" $eff
    echo "--------- Warning ---------"
    echo "Job efficiency is below 80%, perhaps there is a better number of requested CPUs to use in each run!"
    echo ""
else
    echo "REUS efficiency:" $eff
fi
       
#
# Write executable script:
cat <<EOF > $runname.slurm
#!/bin/bash -e

# Finds current block
i=1
j=001

while [ -f GLOB/GLOBLOG_\${j} ] ; do
        i=\$((i+1))
        j=\`printf "%03d\n" \${i}\`
done
#

k=\$((i-1))
l=\`printf "%03d\n" \${k}\`

########################
# Info about the Job:
echo "Job executed in Machine \$HOSTNAME" > $rundir/${runname}_\${j}.blockinfo
echo "Job executed with $ncpus processors" >> $rundir/${runname}_\${j}.blockinfo
echo "Job efficiency of \$eff using $ucpus (MD) and $dcpus (PB) processors per run" >> $rundir/${runname}_\${j}.blockinfo
echo "Job executed in DIR: /tmp/${USER}_REUS\$$ " >> $rundir/${runname}_\${j}.blockinfo

echo "" >> $rundir/${runname}_\${j}.blockinfo
echo -e "Job started on: " >> $rundir/${runname}_\${j}.blockinfo
date >> $rundir/${runname}_\${j}.blockinfo

# line needed to mount /programs using autofs
ls /programs/CpH-MD/ >/dev/null 2>&1; sleep 2

# Copy important files for local directory; first gro is called 
mkdir -p /tmp/${USER}_REUS\$$
cp REUS.py REUS_input.py /tmp/${USER}_REUS\$$
for US in $USs
do
mkdir /tmp/${USER}_REUS\$$/${prefix}\${US}
##### Insert the new CPU number in the pHmdp #####
## First the CPU to use on the MD part
sed -i "s/export nCPU=.*/export nCPU=$ucpus/g" ${prefix}\${US}/${REUSprefix}_${sysname}_001.pHmdp
## Then update the portion of the CPU to use on the DelphiT
sed -i "s/export dCPU=.*/export dCPU=$dcpus/g" ${prefix}\${US}/${REUSprefix}_${sysname}_001.pHmdp


##################
cp ${prefix}\${US}/${sysname}_\${l}.gro /tmp/${USER}_REUS\$$/${prefix}\${US}/${sysname}_\${l}.gro
cp ${prefix}\${US}/CpHMD.sh             /tmp/${USER}_REUS\$$/${prefix}\${US}
cp ${prefix}\${US}/${REUSprefix}_${sysname}{_001.pHmdp,.mdp,.pbp,.top,_GROW.gro,.ndx,_PDBin.pdb,.sites} /tmp/${USER}_REUS\$$/${prefix}\${US}


done

# Enter local directory
cd /tmp/${USER}_REUS\$$
sleep 5

# Run REUS-MD segment:
nice -n 19 /tmp/${USER}_REUS\$$/REUS.py \${j} > ${runname}_\${j}.err 2>&1 

# Copy files and return to shared directory
for US in $USs
do
cp -df ${prefix}\${US}/${sysname}_\${j}* ${rundir}/${prefix}\${US}
if (for f in ${prefix}\${US}/${sysname}_\${j}*; do diff \$f ${rundir}/\$f; done); then
echo "Files from US \${US} were correctly copied." >> $rundir/${runname}_\${j}.blockinfo
else
echo "Error in file copy... please check local files" >> $rundir/${runname}_\${j}.blockinfo
exit 1
fi

### Portion to ensure the compatibility of pull-references ###
if [ \`awk '/pull-group2-name/ {print \$3}' ${prefix}\${US}/${REUSprefix}_${sysname}.mdp\` != \`awk '/pull-group2-name/ {print \$3}' ${rundir}/${prefix}\${US}/${REUSprefix}_${sysname}.mdp\` ];then 
mv ${rundir}/${prefix}\${US}/${REUSprefix}_${sysname}.mdp ${rundir}/${prefix}\${US}/${REUSprefix}_${sysname}_\${j}.mdp.old
cp -f ${prefix}\${US}/${REUSprefix}_${sysname}.mdp ${rundir}/${prefix}\${US}/${REUSprefix}_${sysname}.mdp
echo "mdp file from US \${US} was exchanged." >> $rundir/${runname}_\${j}.blockinfo
fi
#####################################################

done
#####################
mkdir -p ${rundir}/GLOB
cp -df {GLOBLOG,GLOBOUT}_\${j} ${rundir}/GLOB
if (for f in {GLOBLOG,GLOBOUT}_\${j}; do diff \$f ${rundir}/GLOB/\$f; done); then
cd ${rundir}
rm -rf /tmp/${USER}_REUS\$$
else
echo "Error in file copy... please check local files" >> $rundir/${runname}_\${j}.blockinfo
exit 1
fi
#####################

#
echo "" >> $rundir/${runname}_\${j}.blockinfo
echo -e "Job finished on: " >> $rundir/${runname}_\${j}.blockinfo
date >> $rundir/${runname}_\${j}.blockinfo
#
# Launch next job before exiting:
if [ \${i} -lt $Segments ] # usually 1 segment == 1 ns
then
    cd $rundir; ./$runfile
fi

EOF
    chmod +x $runname.slurm

#
    if [ $requeue == 1 ]
    then 
	sbatch --requeue -p $Partition -N 1 -n $ncpus -o $runname.sout -e $runname.serr $runname.slurm
    else
	sbatch -p $Partition -N 1 -n $ncpus -o $runname.sout -e $runname.serr $runname.slurm
    fi
    echo ""
    echo "Job submitted to Partition(s): $Partition with $ncpus Processors"

#
## End of Script
#
exit
