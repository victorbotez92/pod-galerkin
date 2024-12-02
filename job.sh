#!/bin/bash
#SBATCH --job-name=run_pod_glk # suite_to_binaries nom du job
#SBATCH --partition=cpu_med # cpu_long  #fusion_shm_big     # Nom d'une partition pour une exécution cpu
#SBATCH --ntasks=200 #2040          # Nombre total de processus MPI
#SBATCH --hint=nomultithread       # 1 processus MPI par coeur physique (pas d'hyperthreading)
#SBATCH --time=04:00:00            # Temps d’exécution maximum demande (HH:MM:SS)
#SBATCH --output=JobLogs/run.log  # Nom du fichier de sortie
#SBATCH --error=JobLogs/run.err   # Nom du fichier d'erreur (ici commun avec la sortie)
#SBATCH --exclusive  

# /!\ Attention, la ligne suivante est trompeuse mais dans le vocabulaire
# de Slurm "multithread" fait bien référence à l'hyperthreading.

# on se place dans le répertoire de soumission
cd ${SLURM_SUBMIT_DIR}

# nettoyage des modules charges en interactif et herites par defaut
module purge
#module load intel-all/19.0.4
module load parmetis/4.0.3/intel-19.0.3.199-intel-mpi-int32-real64
module load fftw/3.3.10/intel-20.0.4.304
module load petsc/3.18.1/intel-20.0.4.304-intel-mpi
module load intel-parallel-studio/cluster.2020.4/intel-20.0.4.304



set -x
# exécution du code
cp $HOME/POD_GALERKIN/EXECUTABLE/a.exe a_pod_galerking.exe

date
srun ./a_pod_galerking.exe
date
#ls -lrt
