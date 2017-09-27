#!/bin/bash -l

############################### SCRIPT DETAILS ################################
#                                                                             #
# THIS PIPELINE ASSUMES YOU HAVE A DIRECTORY IN WHICH THERE IS A FOLDER       #
# NAMED FAA, WHICH CONTAINS MULTIPLE FAA FILES. THESE FAA FILES CORRESPOND TO #
# ONE GENOME EACH. THE SEQUENCES MUST BE PROTEIN. IN ADDITION, IT ALSO        #
# THAT YOU HAVE EITHER ONE FFN FILE CONTAINING ALL NUCLEOTIDE SEQUENCES FROM  #
# ALL GENOMES, OR AN FFN FOLDER CONTAINING ALL GENOMES IN FFN FILES.          #
#                                                                             #
# PROGRAMS REQUIRED: PYTHON, AFREE, MCL, MUSCLE, MOTHUR                       #
# PYTHON LIBRARIES REQUIRED: BIOPYTHON, NUMPY                                 #
#                                                                             #
# EMAIL OLGAZH@DARTMOUTH.EDU WITH ANY QUESTIONS                    	      #
#                                                                             #
###############################################################################

############################ PARAMETERS FOR SCRIPT ############################
#                                                                             #
# minimum afree score to keep                                                 #
AFREETHRESHOLD=10                                                             #
#                                                                             #
# keep families containing >= x% genomes                                      #
SIZEPERCENT=95                                                                #
#                                                                             #
# discard families with > x avg copies of gene                                #
GENECOPYNUMBER=1.1                                                            #
#                                                                             #
# discard sequences with < x nucleotides aligned                              #
MINALIGNMENT=50                                                               #
#                                                                             #
###############################################################################

# set name as the current directory
NAME=${PWD##*/}

echo "Running pipeline on $NAME"

# check that there is an faa folder
if [ ! -d "faa" ]; then
echo "No faa folder"
exit 1
fi

# check that there are faa files in folder
if [ ! "$(ls faa/*.faa)" ]; then
echo "No faa files found"
exit 1
fi

nGENOMES=(faa/*.faa)
nGENOMES=${#nGENOMES[@]}

# make out folder if not there
if [ ! -d "out" ]; then
mkdir out
fi

# make clusters folder if not there
if [ ! -d "clusters" ]; then
mkdir clusters
fi

# check if there is ffn file with all genes
if [ ! -s $NAME.ffn ]; then
# check if there is ffn folder
if [ -d "ffn" ]; then
# if there is, check if there are ffn files in folder
if [ ! "$(ls ffn/*.ffn)" ]; then
echo "No ffn files"
exit 1
else
# make ffn file with all genes from these ffn files
cat ffn/*.ffn > $NAME.ffn
fi
else
echo "No ffn folder"
fi
fi

cd out

# run afree on all genomes
echo "Running afree on all genomes"
all_vs_all_egm2.py -at $PBS_NUM_PPN -s $AFREETHRESHOLD ../faa/*.faa

if [ ! "$(ls *)" ]; then
echo "No output from afree"
exit 1
fi

cd ..
echo Parsing afree output
parse_afree.py -uo $NAME.mcl out/*

if [ ! -s $NAME.mcl ]; then
echo "MCL did not finish"
exit 1
fi

# run mcl on output from afree
echo "Running MCL on afree output"
mcl $NAME.mcl --abc -I 1.2 -te $PBS_NUM_PPN -o $NAME.clusters

# quickly filter out small clusters
MINSIZE=$(($nGENOMES*$SIZEPERCENT/100))
echo "Filtering out clusters smaller than $MINSIZE genes"
filter_mcl_clusters.py -m $MINSIZE -o $NAME.filtered $NAME.clusters

if [ ! -s $NAME.lengths.pkl.gz ]; then
echo "Building length dictionary"
parse_lengths.py -bl $NAME.lengths.pkl.gz faa/*.faa
fi
if [ ! -s $NAME.lengths.pkl.gz ]; then
echo "Could not create length dictionary"
exit 1
else
echo Filtering out sequences based on protein length
parse_lengths.py -l $NAME.lengths.pkl.gz -o len $NAME.filtered
fi

if [ ! -s $NAME.filtered.len ]; then
echo "No length filtered file"
exit 1
fi

if [ ! -s $NAME.genomes ]; then
cd faa
echo "Building dictionary of gene IDs for each genome"
get_gene_ids.py -o ../$NAME.genomes *.faa
cd ..
fi
if [ ! -s $NAME.genomes ]; then
echo "Can't make genome file"
exit 1
fi

# filter high gene copy count clusters
echo "Filtering out families with high gene copy count"
check_clusters.py -m $NAME.filtered.len -i $NAME.genomes -e $GENECOPYNUMBER

if [ -s $NAME.filtered.len.filtered ]; then
mv $NAME.filtered.len.filtered $NAME.final
else
echo "Failed to filter gene copy number"
exit 1
fi

cd clusters

# get fa files based on filtered clusters
echo "Creating fasta files with nucleotide sequences"
select_sequences_in_fasta.py -i ../$NAME.ffn -o $NAME.fa -m ../$NAME.final

# make sure there were fa files made
if [ ! "$(ls *.fa)" ]; then
echo "No fa files were created"
exit 1
fi

# align with MUSCLE and create mothur batch file
echo "Aligning nucleotide sequences with MUSCLE"
run_mothur_pipeline.py -bmt $PBS_NUM_PPN -i 16 *.fa

if [ ! "$(ls *.aligned.fa)" ]; then
echo "Did not perform alignment"
exit 1
fi

# remove small or non-overlapping sequences
echo "Filtering alignments for < $MINALIGNMENT overlapping nucleotides"
parse_alignment.py -t $MINALIGNMENT -p $PBS_NUM_PPN *.aligned.fa

if [ ! "$(ls *.aligned.fa.out)" ]; then
echo "Did not filter out small alignments"
exit 1
fi

# run mothur on fixed alignments
echo "Running mothur on alignments"
run_mothur_pipeline.py -f *.aligned.fa.out

if [ ! "$(ls *.list)" ]; then
echo "No list files to parse"
exit 1
fi

# parse mothur results
echo "Parsing mothur results"
parse_mothur.py -o ../$NAME.txt *.list

cd ..

# get rid of extra stuff in column labels
sed -i "s/\.aligned\.fa\.filter\.phylip\.fn\.list//g" $NAME.txt

exit 0
