#!/bin/bash

# Define arrays for each parameter
ROIs=("dmPFC" "dlPFC")
freqBands=("Theta" "Beta" "HFA")
alignEvents=("S" "R")

# Loop through all combinations of parameters
for ROI in "${ROIs[@]}"; do
    for freqBand in "${freqBands[@]}"; do
        for alignEvent in "${alignEvents[@]}"; do
            # Set the job name and function based on ROI
            if [ "$ROI" == "dmPFC" ]; then
                jobName="Permutation_Test_LMM_dmPFC"
                functionCall="PermTestLMMMPFC"
            else
                jobName="Permutation_Test_LMM_dlPFC"
                functionCall="PermTestLMMLPFC"
            fi

            # Create a unique job name for each combination
            uniqueJobName="${jobName}_${freqBand}_${alignEvent}"

            # Submit the job
            sbatch <<EOT
#!/bin/bash

#SBATCH --job-name=$uniqueJobName
#SBATCH --output=/data/user/anaskhan/PRJ_Stroop/cluster_logs/logfile_%A_%a.out
#SBATCH --error=/data/user/anaskhan/PRJ_Stroop/cluster_errors/results_%A_%a.err
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G

# Load MATLAB module
module load rc/matlab

# Export additional variables
export ROI="$ROI"
export freqBand="$freqBand"
export alignEvent="$alignEvent"

# Run the analysis function
matlab -nodisplay -nosplash -r "$functionCall('$ROI','$freqBand','$alignEvent'); exit"
EOT
        done
    done
done