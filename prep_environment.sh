#!/bin/sh
# common variables
export DATA_DIR=/projects/cees2/in_progress/cees_workshop
export REF=/projects/cees2/in_progress/cees_workshop/ref/house_sparrow_genome_assembly-18-11-14_masked.fa
# set up CEES stuff
export PATH=/projects/cees/bin:/projects/cees/scripts:$PATH
umask 0002
module use --append /projects/cees/bin/modules