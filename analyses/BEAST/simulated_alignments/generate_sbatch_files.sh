#!/bin/bash
# GENERATES SBATCH FILES FOR RUNNING CHAINS ON MIDWAY
# For each scenario:
for scenario in {1,2a,2b,2c,2d,2e,3a,3b,3c,3d}
do
    # For each replicate
    for replicate in {1..4}
    do   
    echo "#!/bin/bash" > scenario${scenario}_rep${replicate}.sbatch
    echo "#SBATCH --job-name=scn${scenario}_rep${replicate}" >> scenario${scenario}_rep${replicate}.sbatch
    echo "#SBATCH --output=scn${scenario}_rep${replicate}.out" >> scenario${scenario}_rep${replicate}.sbatch
    echo "#SBATCH --error=scn${scenario}_rep${replicate}.err" >> scenario${scenario}_rep${replicate}.sbatch
    echo "#SBATCH --time=1500:00:00" >> scenario${scenario}_rep${replicate}.sbatch
    echo "#SBATCH --partition=cobey" >> scenario${scenario}_rep${replicate}.sbatch
    echo "#SBATCH --nodes=1" >> scenario${scenario}_rep${replicate}.sbatch
    echo "#SBATCH --ntasks-per-node=1" >> scenario${scenario}_rep${replicate}.sbatch
    
    echo ""  >> scenario${scenario}_rep${replicate}.sbatch
    
    echo "module load beast/1.8" >> scenario${scenario}_rep${replicate}.sbatch

    echo "cd /project/cobey/mvieira/evolution_of_mutability/results/BEAST/simulated_alignments/scenario${scenario}_rep${replicate}/" >> scenario${scenario}_rep${replicate}.sbatch
               
    echo "beast /project/cobey/mvieira/evolution_of_mutability/analyses/BEAST/simulated_alignments/scenario${scenario}_rep${replicate}/scenario${scenario}_rep${replicate}.xml" >> scenario${scenario}_rep${replicate}.sbatch

    mv scenario${scenario}_rep${replicate}.sbatch scenario${scenario}_rep${replicate}/
    done    
done

