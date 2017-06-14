#!/bin/bash
for clone in {CH103,CH103L,VRC26,VRC26int,VRC26L,VRC01_01,VRC01_19,VRC01_13}
do
   for prior in {constant,exponential,logistic}
   do
    echo "#!/bin/bash" > ${clone}_${prior}_process_output.sbatch
    echo "#SBATCH --job-name=${clone}_${prior}_process_output" >> ${clone}_${prior}_process_output.sbatch
    echo "#SBATCH --output=${clone}_${prior}_process_output.out" >> ${clone}_${prior}_process_output.sbatch
    echo "#SBATCH --error=${clone}_${prior}_process_output.err" >> ${clone}_${prior}_process_output.sbatch
    echo "#SBATCH --time=2:00:00" >> ${clone}_${prior}_process_output.sbatch
    echo "#SBATCH --partition=cobey" >> ${clone}_${prior}_process_output.sbatch
    echo "#SBATCH --nodes=1" >> ${clone}_${prior}_process_output.sbatch
    echo "#SBATCH --ntasks-per-node=1" >> ${clone}_${prior}_process_output.sbatch
    echo "#SBATCH --mem-per-cpu=16000" >> ${clone}_${prior}_process_output.sbatch
    
    
    echo ""  >> ${clone}_${prior}_process_output.sbatch
    echo "module load java" >> ${clone}_${prior}_process_output.sbatch
    echo "cd /project/cobey/mvieira/evolution_of_mutability/analyses/BEAST/observed_lineages" >> ${clone}_${prior}_process_output.sbatch
    echo "./process_output.sh $clone $prior" >> ${clone}_${prior}_process_output.sbatch
    
    mv ${clone}_${prior}_process_output.sbatch ${clone}_${prior}
   
   done

done