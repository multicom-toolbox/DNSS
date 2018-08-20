#!/bin/bash
###############################################################################
# Name %n  : generate-pssm.sh
# Desc %d  : Generate and save a PSSM for a query using psi blast and 
#  suggested parameters from BMC Bioinformatics 2008, 9(Suppl 12):S12 
#  doi:10.1186/1471-2105-9-S12-S12
# Input %i : A fasta file and two output file names
# Output %o: The pssm and a report 
#
# Author: Jesse Eickholt
# URL: http://merit.oryo.us
# Date: Wed Jan 25 2012
# Latest Modified Date: 11/02/2016 by Jie Hou
###############################################################################
GLOBAL_PATH='/space1/jh7x3/DNSS_release/';

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <seq.fasta> <out_fname1(pssm)> <out_fname2(report)>";
  exit;
fi

FASTA=$1
PSSM=$2
REPORT=$3

DB="nr90"

# Get path to script, works generally but not if code is sourced
#abspath=$(cd ${0%/*} && echo $PWD/${0##*/})
#BASEPATH=`dirname $abspath`;

# Get paths using PATHFINDER
BLAST_PATH="$GLOBAL_PATH/programs/ncbi-blast-2.2.25+/bin"; #`$PATHFINDER/get-path.sh NCBI_BLAST`/bin";
BLAST_NR_DB="$GLOBAL_PATH/nr_database"; #`$PATHFINDER/get-path.sh BLAST_DB`";
BLASTMAT="$GLOBAL_PATH/programs/ncbi-blast-2.2.25+/matrices"; #`$PATHFINDER/get-path.sh BLAST_MAT`";

echo "Running PSI-Blast with $DB...";
#$BLAST_PATH/psiblast -query "$FASTA" -evalue .001 -inclusion_ethresh .002 -db "$BLAST_NR_DB/$DB" -num_iterations 3 -outfmt "0"\
#  -out $REPORT -num_alignments 2000 -out_ascii_pssm $PSSM

$BLAST_PATH/psiblast -query "$FASTA" -evalue .001 -inclusion_ethresh .002 -db "$BLAST_NR_DB/$DB" -num_iterations 3 -outfmt "0"\
  -out $REPORT -seg yes -out_ascii_pssm $PSSM


