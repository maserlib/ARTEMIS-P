#!/bin/bash
#SBATCH --job-name=ARTEMIS-P
#SBATCH --nodes=1 --ntasks-per-node=16 
#SBATCH --time=120 
#SBATCH --partition=padc_medium
#SBATCH --mail-user=baptiste.cecconi@obspm.fr
#SBATCH --mail-type=ALL
#SBATCH --mem=1gb 
#SBATCH --tmp=10gb 

printenv

# Definition du repertoire local scratch
SCRATCH=/scratch/$USER/run.${SLURM_JOBID}
HOME_DIR=$HOME

# Creation du repertoire $SCRATCH sur le noeud reserve
mkdir -p $SCRATCH
cd $SCRATCH
# Recopie des donnees et de l'executable dans le repertoire $SCRATCH
cp ${HOME_DIR}/ARTEMIS-P/RaytracingCode/raytracing $SCRATCH
cp ${HOME_DIR}/ARTEMIS-P/RaytracingCode/init_environ.txt $SCRATCH
cp ${HOME_DIR}/ARTEMIS-P/RaytracingCode/init_raytracing.txt $SCRATCH

./raytracing

mv ${SCRATCH} /poubelle/$USER
rm -rf ${SCRATCH}

exit 0
 
