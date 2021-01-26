#!/bin/bash
#SBATCH --job-name=GENESIS_results
#SBATCH --cpus-per-task=8
#SBATCH --error=GENESIS_results.stdout
#SBATCH --error=GENESIS_results.stderri
#SBATCH --mem 90000
#SBATCH -p defq

eval #SBATCH --account=barnecaapa=${USER}"

module load R-3.3.1-gcc-6.1.0-tpkh25pakjkfxwqth7kp3gxl4grlf6c2
module load torque

R CMD BATCH genesis_clean_qqman.R