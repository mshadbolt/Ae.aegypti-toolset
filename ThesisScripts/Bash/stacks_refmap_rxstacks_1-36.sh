#!/bin/bash

### Modified Stacks refmap script ###
#
# Originally authored by Rasic, G. and Filipovic, I.
# Published in Rasic et al. 2014 DOI: 10.1186/1471-2164-15-275
#
# Modified by Marion Shadbolt
# email: marion.shadbolt@gmail.com
# Last edited: 27/10/2016
#
## Additions: ##
# 1. Allows input of bam files or incorporates zipping and unzipping of bowtie files
#
# 2. Incorporates rxstacks correction module expected to reduce genotyping error)
#
# Note: Tested to work with Stacks version 1.36, parameters may need editing for later/earlier versions
# By default does not create or populate database used for web interface as current version (1.19) not compatible with
# new database structure
#
# Additions have been highlighted with block comments
#
########################################

#environment and data set parameters
_start_time=$(date +%Y%m%d_%H%M%S)
_number_of_processors=$(nproc)
_user="marion"
_data_set="inputs"; # Name of the dataset (dataset folder).
_scripts_file_path=~/PEARG/$_user/PycharmProjects/MazTrimAlign; # Path to the script files folder.
_input_files_path=$_scripts_file_path/$_data_set; # Path to the input files.
_output_files_path=$_input_files_path/$_start_time
_popmap_file="AllPop.map" # Population map file name
_log_file_path=$_scripts_file_path/log_files; # Path for the log file
_input_file_format="bam" #bam or bowtie files accepted
_stacks_input_files_filter="*."$_input_file_format # file extension of mapped reads

# Stacks refmap parameters
_minimum_stack_depth=5
_batch_number=1
_stacks_binaries_path=~/PEARG/$_user/tools/bin

# Stacks rxstacks parameters. By default peforms haplotype pruning and uses upper bounded SNP error model
# Params can be added/changed below, line 77
_conf_lim=0.25
_err_lim=0.15
_lnl_lim="-10.0"

# Stacks populations parameters - valid commands up to Stacks v.1.36.
# Arguments below are not exhaustive, Add to commands at line 107 as desired
_use_popmap_file="yes" # Use poulation map file (yes / no)
_popmap_file_path=$_input_files_path; # Path for the population map file
_num_pops=3  # number of populations for -p in populations module
_per_pop=75 # percentage of individuals in a population required to process a locus
_snp_output="" # if empty string, outputs every SNP at a locus, can be altered to "--write_single_snp" or "--write_random_snp"
_outputs=( --vcf ) # array of desired output formats, common ones: [ --vcf --vcf_haplotypes --plink ]

# Uncomment if you want to enable database creation
#_database_name=$_data_set"_"$_start_time"_radtags"
#_cor_database="cor_"$_database_name
#_uncor_database="uncor_"$_database_name


clear
read -p "Press [Enter] key to start..."
echo

if [ ! -d "$_log_file_path" ]; then echo $(date +%H:%M:%S) "Creating: "$_log_file_path; mkdir -p $_log_file_path; fi

{ time {

date
if [ ! -d "$_output_files_path" ]; then echo $(date +%H:%M:%S) "Creating: "$_output_files_path; mkdir -p $_output_files_path; fi

if [ $_input_file_format = "bowtie" ];
    then
    cd $_input_files_path
    echo "Unzipping bowtie input files..."
    gunzip *.bowtie.gz
elif [ ! $_input_file_format = "bam" ];
    then
    echo "File format "$_input_file_format" not recognised, must be bam or bowtie"
fi

cd $_output_files_path

# Run Stacks refmap wrapper with user-defined parameters
ref_map.pl \
$( cd $_input_files_path; for file in $_stacks_input_files_filter; do echo -ne "-s "$_input_files_path/$file" "; done) \
$( if [ $_use_popmap_file = "yes" ]; then echo -ne "-O "$_popmap_file_path/$_popmap_file" "; fi) \
-o $_output_files_path \
-T $_number_of_processors \
-m $_minimum_stack_depth \
-b $_batch_number \
-S 

# Re-zip bowties
if [ $_input_file_format = "bowtie" ];
    then
    cd $_input_files_path
    echo "Re-zipping bowtie input files ..."
    gzip -9 *.bowtie
    echo "...done."
fi

# Make directory for corrected Stacks catalogue files
cd $_output_files_path
mkdir corr
_corr_files_path=$_output_files_path/corr

echo "Running rxstacks correction module..."
rxstacks -b $_batch_number -P $_output_files_path -o $_corr_files_path \
--conf_lim $_conf_lim \
--prune_haplo \
--model_type bounded --bound_high $_err_lim \
--lnl_lim $_lnl_lim --lnl_dist \
-t $_number_of_processors --verbose 

# Move original outputs files so they don't interfere with steps below
mkdir outputs
mv batch* outputs

_stacks_output_files_filter="*tags.tsv"

echo "Rematching corrected loci against catalogue..."
cstacks -b $_batch_number -g -p $_number_of_processors -o $_corr_files_path \
$( cd $_corr_files_path; for file in $_stacks_output_files_filter; do echo -ne " -s "$_corr_files_path/${file%.*.*}" "; done)

cd $_corr_files_path

sstacks -b 1 -p $_number_of_processors -c  $_corr_files_path/batch_$_batch_number \
$( cd $_corr_files_path; for file in $_stacks_output_files_filter; do echo -ne " -s "$_corr_files_path/${file%.*.*}" -o "$_corr_files_path" "; done)

####### Remove SNPs in restriction cut site (first four bases) #######

# Filter sumstats file to get non-duplicated list of loci with position > 4 in tag locus
grep "^[^#]" $_output_files_path/outputs/batch_1.sumstats.tsv | \ 
	awk '{if ($5 > 3) {print $2 "\t" $3 "\t" $4 "\t" $5}}' | \
	awk '!x[$0]++{print $1 "\t" $4}' > $_corr_files_path/whitelist.txt

#####################################################################

echo "Running populations module without cut site snps..."
populations -b 1 -P $_corr_files_path -M $_popmap_file_path/$_popmap_file \
	-W $_corr_files_path/whitelist.txt \
	-t $_number_of_processors $_outputs $_snp_output \
	-p $_num_pops -r $_per_pop

# Uncomment if you want to enable database creation - NB: Stacks 1.36 database will not work in web interface of PEARG as of 25/10/16
#mysql -e "create database "$_cor_database
#mysql -e "create database "$_uncor_database_name

#mysql $_cor_database < /home/marion/PEARG/marion/tools/stacks-1.36/sql
#mysql $_uncor_databse < /home/marion/PEARG/marion/tools/stacks-1.36/sql

#load_radtags.pl -D $_cor_database -b 1 -p $_cor_files_path -B \
#-e "Stacks 1.36 m=5 Gv10-Vi14 Corrected data" -c -t population -M 

#-B $_database_name \

#-D ''$_data_set' - Stacks-1.19 m='$_minimum_stack_depth'' \
#-a $(date +%Y-%m-%d) 

#-O $_popmap_file_path/$_data_set.map
#-a YYYY-MM-DD
# -d
# -O $_input_files_path/PopMapEmRo.map \
# -n 1 -m 7

#------------------
#database backup
#echo -ne "Backing up ["$_database_name"] database into "$_output_files_path/$_database_name.sql.gz" ... "
#mysqldump --databases $_database_name \
#| gzip > $_output_files_path/$_database_name.sql.gz
echo "done."
date

} } 2>&1 | tee $_log_file_path/$_start_time"_"$_data_set"_stacks_refmap.log"
