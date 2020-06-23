#!/usr/bin/env bash
#SBATCH --mem=16G

module load ddsclient
ddsclient upload -p ddh-com-data data
