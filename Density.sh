#!/bin/bash                               
#$ -l rmem=48G 
#$ -l h_rt=100:00:00                     
#$ -M n.bulgakova@sheffield.ac.uk
#$ -m bea

cd /data/md1nbu/Density/Validation

module load apps/matlab/2017b          

matlab -nodesktop -r ‘validationkmeanspar’
