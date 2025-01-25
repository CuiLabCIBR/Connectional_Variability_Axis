#!/bin/bash
#SBATCH -p q_cn_2
#SBATCH --ntasks=1 # Run a single serial task
#SBATCH --cpus-per-task=1
#SBATCH --job-name=calc_FC

# Define the base directory paths for output and log files
BASE_DIR="/ibmgpfs/cuizaixu_lab/yanghang/projects/fc_variability_axis/hcp"
OUT_DIR="$BASE_DIR/out"
LOG_DIR="$BASE_DIR/log"
SUBLIST="$BASE_DIR/hcp_sublist.txt"

# Create the output and log directories if they don't exist
echo "Checking directories..."
for DIR in "$OUT_DIR" "$LOG_DIR"; do
  if [ ! -d "$DIR" ]; then
    echo "Creating directory: $DIR"
    mkdir -p "$DIR" || { echo "Error creating directory $DIR"; exit 1; }
  else
    echo "Directory exists: $DIR"
  fi
done

# Loop over the subject range and submit jobs with dynamic output and error paths
echo "Submitting jobs..."

# Use cat to read the contents of sublist.txt and pipe it into the while loop
cat "$SUBLIST" | while IFS= read -r subj
do
   # Skip empty lines or comments (lines starting with #)
   if [[ -z "$subj" || "$subj" == \#* ]]; then
      continue
   fi
   
   # Dynamically create the output and log file paths
   OUT_PATH="$OUT_DIR/${subj}.out"
   LOG_PATH="$LOG_DIR/${subj}.log"
   
   # Submit the job with the dynamically created output and log paths
   echo "Submitting job for subject: $subj"
   sbatch -o "$OUT_PATH" -e "$LOG_PATH" "$BASE_DIR/calculate_FC.sh" "$subj" || { echo "Job submission failed for subject $subj"; exit 1; }
done

echo "All jobs submitted successfully."
