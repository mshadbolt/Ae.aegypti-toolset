#!/bin/bash

### Modified trim_and_align script ###
# Originally authored by Rasic, G. and Filipovic, I.
# Published in Rasic et al. 2014 DOI: 10.1186/1471-2164-15-275
#
# Modified by Marion Shadbolt
# email: marion.shadbolt@gmail.com
# Last edited: 27/10/2016
# 
## Dependencies:
#
# Expects the following programs to be installed and in the path:
#	1. bowtie 0.12.7 or later
#	2. samtools 1.3.1 or later
#	3. python 3.4 or earlier. 
#		- Uses modules: glob, sys, subprocess, collections, os, multiprocessing, time
#
## Additions: ##
# 1. BAM output from bowtie:
# Default alignment behaviour is to output a BAM file for each individual then convert
# to a sorted, indexed bam file that only contains aligned reads. Previous behaviour can be
# enabled by toggling _include_sam_output parameter to "no".
#
# 2. Incorporate merging:
# Automatically merge files so they are ready for stacks input. Merging removes paired and singles files.
# Toggled on and off by _merge_files parameter
#
# 3. Incorporate Male/Female diagnosis
# Automatically performs Male/Female diagnosis by calling the MaleDiag.py script.
# Based on sequences and method published in Fontaine et al. 2016 doi: http://dx.doi.org/10.1101/060061
# Assumes bowtie index of male specific sequences is in the _bowtie_index_path called "Aaegypt-male".
# If reference is named differently, change in MaleDiag.py script.
#
# 4. Perform check of required files before commencing script, does not proceed unless all required files are present.
#
# 5. Automatically detects instrument name so reads from different sequencing runs can done at once.
#
# 6. Creates a tab-delimited stats file that records aligned versus non-aligned reads for each sample against each 
# reference genome used. Saved in input_files_path
#
# Modifications:
# 1. Added -p parameter to gzip and gunzip commands so that _number_of_processors used is not exceeded. Could
# not work out how to adjust this for the trimming stage.
#
# Additions have been highlighted with block comments
########################################

#environment and data set parameters
_start_time=$(date +%Y%m%d_%H%M%S)
_time_stamp=$_start_time
_number_of_processors=8; # or set as $(nproc) to use all available CPUs.
_data_set="inputs"; #Name of the dataset (dataset folder).
_scripts_file_path=~/PEARG/; #Path to the script files folder.
_input_files_path=~/PEARG/; #Path to the input files.
_working_files_path=/home/working/RAMDisk/$_data_set; #Path to the working folder.
_output_files_path=$_input_files_path/$_time_stamp; #Path for all output files.
_log_file_path=$_scripts_file_path/log_files; #Path for the log file

#synchronizing and sorting parameters
_sort_tmp_path=/home/working/RAMDisk/$_data_set/tmp; #location for creating tmp files required by "sort_paired_reads.sh" script.
_verbose_sorting="True"; # Additional parameter for the "sort_paired_reads.sh" script. If it is "True" script will provide stats. 
_read_header_splitter=" "

#trimmer parameters
_skip_trimming="no" # yes / no
_R1_overhang="CATG"
_R2_overhang="AATT"
_R1_first_base_to_keep=1
_R1_last_base_to_keep=90
_R1_quality_threshold=20
_R2_first_base_to_keep=1
_R2_last_base_to_keep=90
_R2_quality_threshold=20
_input_quality_encoding_base=33
_output_quality_encoding_base=33
_quality_encoding_base=$_output_quality_encoding_base

_alignment_software="bowtie" # Choose between bowtie or bowtie2

#bowtie parameters
_bowtie_path=/home/working/bowtie
_bowtie_indexes=( wolbachia AaegL1 ); # Ordered list of indexes to go through ( phiX174 wolbachia Aaeg_mtDNA AaegL1)
_bowtie_index=""
_bowtie_index_path=$_bowtie_path/indexes
_min_size=$_R1_last_base_to_keep
_max_size=400
_maximum_number_of_mismatches_allowed=3
_maximum_number_of_reportable_alignments=1
_include_sam_output="yes" # yes/no. "yes" creates SAM directly from bowtie and converts to sorted, indexed bam file

######New additions######

#samtools parameters
_path_to_reference=/home/working/bowtie/genomes/

#Merge aligned reads?
# If "Y", will merge aligned reads and remove paired and singles files.
_merge_aligned="Y"

#Male/Female diagnosis settings
_sex_diag="Y" # If "Y", will call the sex diagnosis script
_pop_map_name="AllPop.map" #tab-delimited text file linking samples to populations, Should be in your input file path
_single_pop="N" #If only one pop in input files, negates need for pop map
_single_pop_name="PopName"

#check for existence of additional required files:
clear

echo "Checking required files..."
if test -x $_scripts_file_path/sort_paired_reads.sh; then
    echo "Sort paired reads script found."
else
    echo "sort_paired_reads.sh not found or not executable, ensure it is in your scripts directory."
    echo "exiting..."
    exit
fi

if [ $_sex_diag = "Y" ] && [ -x $_scripts_file_path/MaleDiag.py ]; then
    echo "Sex diagnosis script found."
fi

if  [ $_sex_diag = "Y" ] && [ ! -x $_scripts_file_path/MaleDiag.py ]; then
    echo "MaleDiag.py not found or not executable, ensure it is in your scripts directory."
    echo "exiting..."
    exit
fi

if [ $_sex_diag = "Y" ] && [ $_single_pop = "N" ] && [ -f $_input_files_path/$_pop_map_name ]; then
    echo "Pop map file found."
fi

if [ $_sex_diag = "Y" ] && [ $_single_pop = "N" ] && [ ! -f $_input_files_path/$_pop_map_name ]; then
    echo "Pop map file not found, ensure it is in your input directory (or check input dir is correct)."
    echo "exiting..."
    exit
fi

if [ $_single_pop = "Y" ];
	then
	echo "All samples are from a single population, labeling with "$_single_pop_name
fi

if [ ! -d "$_log_file_path" ]; then
    mkdir $_log_file_path
fi

# Create stats file
_stats_file=$_input_files_path/$_start_time"_"$_data_set"_""alignment_stats.tsv"
echo -e "index\tsample\tunaligned_paired\taligned_paired\tunaligned_R1\taligned_R1\tunaligned_R2\taligned_R2" >> $_stats_file

######End additions######


echo "****************************************************"
echo ""
echo "Warning!"
echo "Script for trimming and alignment is about to begin."
echo "" 
read -p "Press [Enter] key to start or [Ctrl]+[c] to stop."
echo
echo

{ time {
date
echo $(date +%H:%M:%S)" Executing script ..."; echo

# Trim R1 and R2 reads:
_trimmed_output_files_path=$_output_files_path/trimmed
if [ ! -d "$_output_files_path" ]; then echo $(date +%H:%M:%S) "Creating: "$_output_files_path; mkdir $_output_files_path; fi
if [ ! -d "$_trimmed_output_files_path" ]; then echo $(date +%H:%M:%S) "Creating: "$_trimmed_output_files_path; mkdir $_trimmed_output_files_path; fi
if [ ! -d "$_working_files_path" ]; then echo $(date +%H:%M:%S) "Creating: "$_working_files_path; mkdir $_working_files_path; fi
if [ ! -d "$_sort_tmp_path" ]; then echo $(date +%H:%M:%S) "Creating: "$_sort_tmp_path; mkdir $_sort_tmp_path; fi
#if [ ! -d "$_working_files_path/trimmed" ]; then echo "Creating: "$_working_files_path/trimmed; mkdir $_working_files_path/trimmed; fi

#Trim and filter all R1 reads
_R1_input_files_filter=*_1.fq.gz
_input_files_filter=$_R1_input_files_filter
_first_base_to_keep=$_R1_first_base_to_keep
_last_base_to_keep=$_R1_last_base_to_keep
_quality_threshold=$_R1_quality_threshold
_trimmed_fragment_size=$(($_last_base_to_keep - $_first_base_to_keep + 1))
_overhang=$_R1_overhang
_sed_overhang_commands=$( for i in $(seq 0 $(( ${#_overhang}-1 )) ); do echo -ne "-e '2~4s/^(.{$i})N(.*)$/\\1"${_overhang:$i:1}\\2/"' "; done )

cd $_input_files_path
for file in $_input_files_filter;
do 
	if [ -r "$file" -a $_skip_trimming = "no" ];
	then
		#echo $_SED_overhang_commands
		# Wait for 1 second if the number of running jobs is greater or equal than a number of availiable processors and then check again
		while [ $(jobs | wc -l) -ge $_number_of_processors ]; do sleep 1; done

		# Since the number of running jobs is less than a number of available processors spawn a job 
		echo $(date +%H:%M:%S) "Trimming "$file" from "$_first_base_to_keep"bp to "$_last_base_to_keep"bp ["$_trimmed_fragment_size"bp] using "$_quality_threshold" as quality threshold  ... "
		( gunzip -d -c -f -k -p $_number_of_processors ./$file \
		| fastq_quality_converter -n -Q$_input_quality_encoding_base \
		| fastq_quality_converter -a -Q$_output_quality_encoding_base \
		| eval sed -r -e '4~4s/\%/\&/g' $_sed_overhang_commands \
		| fastx_trimmer -f $_first_base_to_keep -l $_last_base_to_keep -Q$_quality_encoding_base -v \
		| fastq_quality_trimmer -t $_quality_threshold -l $_trimmed_fragment_size -Q$_quality_encoding_base -v \
		| gzip -f -p $_number_of_processors 1> $_trimmed_output_files_path/${file%%.*}.trimmed.fq.gz ; echo $(date +%H:%M:%S) "Done with "$file; echo) &
	fi
done 



#Trim and filter all R2 reads
_R2_input_files_filter=*_2.fq.gz
_input_files_filter=$_R2_input_files_filter
_first_base_to_keep=$_R2_first_base_to_keep
_last_base_to_keep=$_R2_last_base_to_keep
_quality_threshold=$_R2_quality_threshold
_trimmed_fragment_size=$(($_last_base_to_keep - $_first_base_to_keep + 1))
_overhang=$_R2_overhang
_SED_overhang_commands=$( for i in $(seq 0 $(( ${#_overhang}-1 )) ); do echo -ne "-e '2~4s/^(.{$i})N(.*)$/\\1"${_overhang:$i:1}\\2/"' "; done )

cd $_input_files_path
for file in $_input_files_filter;
do 
	if [ -r "$file" -a $_skip_trimming = "no" ];
	then
		#echo $_SED_overhang_commands
		# Wait for 1 second if the number of running jobs is greater or equal than a number of availiable processors and then check again
		while [ $(jobs | wc -l) -ge $_number_of_processors ]; do sleep 1; done

		# Since the number of running jobs is less than a number of available processors spawn a job
		echo $(date +%H:%M:%S) "Trimming "$file" from "$_first_base_to_keep"bp to "$_last_base_to_keep"bp ["$_trimmed_fragment_size"bp] using "$_quality_threshold" as quality threshold  ... "
		( gunzip -d -c -f -k -p $_number_of_processors ./$file \
		| fastq_quality_converter -n -Q$_input_quality_encoding_base \
		| fastq_quality_converter -a -Q$_output_quality_encoding_base \
		| eval sed -r -e '4~4s/\%/\&/g' $_sed_overhang_commands \
		| fastx_trimmer -f $_first_base_to_keep -l $_last_base_to_keep -Q$_quality_encoding_base -v \
		| fastq_quality_trimmer -t $_quality_threshold -l $_trimmed_fragment_size -Q$_quality_encoding_base -v \
		| gzip -f -p $_number_of_processors 1> $_trimmed_output_files_path/${file%%.*}.trimmed.fq.gz ; echo $(date +%H:%M:%S) "Done with "$file; echo) &
	fi
done 



wait
if [ $_skip_trimming = "no" ]; then echo $(date +%H:%M:%S)" All trimming is done."; fi
echo; echo

#if [ ${#_bowtie_indexes[@] = 0 ]; then echo $(date +%H:%M:%S) "List of indexes is empty. Exiting!"; exit; fi
echo $(date +%H:%M:%S) "List of indexes: "${_bowtie_indexes[@]}
_bowtie_index=( trimmed "${_bowtie_indexes[@]}" )
#for _index in ${_bowtie_indexes[@]}; do echo "Index: "$_index; done
for i in $(seq 1 $((${#_bowtie_index[@]}-1))); 
do 
	echo; date; echo $(date +%H:%M:%S) " Index["$((i))"]: "${_bowtie_index[$((i))]}

	#Prepare trimmed and filtered R1 and R2 reads for aligning
	#read -p "Press [Enter] key to continue ..."

	# Prepare directory structure
	if 	[ ${_bowtie_index[$((i-1))]} = "trimmed" ]; 
	then	_current_input_files_path=$_input_files_path/$_time_stamp/${_bowtie_index[$((i-1))]}
	else	_current_input_files_path=$_input_files_path/$_time_stamp/${_bowtie_index[$((i-1))]}/nonaligned
	fi
	_current_output_files_path=$_input_files_path/$_time_stamp/${_bowtie_index[$((i))]}
	_current_working_output_files_path=$_working_files_path/${_bowtie_index[$((i-1))]}
	_current_aligned_output_files_path=$_working_files_path/${_bowtie_index[$((i))]}

	if [ ! -d "$_current_input_files_path" ]; then echo "No input directory found. Aborting!"; exit; fi
	cd $_current_input_files_path
	_input_files_filter=*_1.*.gz 
	for file in $_input_files_filter;
	do 
		cd $_current_input_files_path
		if [ ! -r "$file" ]; then echo "No input files found. Aborting!"; exit; fi

		_read=${file%%_1.*}
		_read1=${file%%_1.*}_1.${_bowtie_index[$((i-1))]}.fq
		_read2=${file%%_1.*}_2.${_bowtie_index[$((i-1))]}.fq

		if 	[ ${_bowtie_index[$((i-1))]} = "trimmed" ]; 
		then
			echo -ne $(date +%H:%M:%S) "Decompressing "$_read1.gz"  ... "
			gunzip -d -c -f -k -p $_number_of_processors $_current_input_files_path/$_read1.gz 1>$_working_files_path/$_read1
			echo " done."
	
			echo -ne $(date +%H:%M:%S) "Decompressing "$_read2.gz"  ... "
			gunzip -d -c -f -k -p $_number_of_processors $_current_input_files_path/$_read2.gz 1>$_working_files_path/$_read2
			echo " done."
		else
			echo -ne $(date +%H:%M:%S) "Decompressing "${file%%_1.*}_1.fq.gz"  ... "
			gunzip -d -c -f -k -p $_number_of_processors $_current_input_files_path/${file%%_1.*}_1.fq.gz 1>$_working_files_path/${file%%_1.*}_1.${_bowtie_index[$((i-1))]}.fq
			echo " done."
	
			echo -ne $(date +%H:%M:%S) "Decompressing "${file%%_1.*}_1.fq.gz"  ... "
			gunzip -d -c -f -k -p $_number_of_processors $_current_input_files_path/${file%%_1.*}_2.fq.gz 1>$_working_files_path/${file%%_1.*}_2.${_bowtie_index[$((i-1))]}.fq
			echo " done."
		fi	
		date
		
		#### Addition: Automatically detect Instrument name ####		
		echo -ne $(date +%H:%M:%S) "Detecting instrument name  ... "
		_first_line=($(head -n 1 $_working_files_path/${file%%_1.*}_2.${_bowtie_index[$((i-1))]}.fq))	
		_space_split=($(echo $_first_line | tr ":" "\n"))
		_instrument=${_space_split[0]}
		echo "Instrument is "$_instrument"."
		#######################################################
		
		echo $(date +%H:%M:%S) "Matching and sorting reads between "$_read1" and "$_read2"  ... "
		$_scripts_file_path/sort_paired_reads.sh \
		$_working_files_path/$_read1 \
		$_working_files_path/$_read2 \
		$_instrument \
		$_verbose_sorting \
		$_sort_tmp_path \
		$_read_header_splitter
		echo " done.";


		#Aligning pairs...
		#read -p "Press [Enter] key to continue ..."
		_bowtied_output_files_path=$_working_files_path
		mkdir $_bowtied_output_files_path/aligned
		mkdir $_bowtied_output_files_path/nonaligned
		echo $(date +%H:%M:%S) "Aligning "${_read1%.*.*}.${_bowtie_index[$((i-1))]}.matched.fq" and "${_read2%.*.*}.${_bowtie_index[$((i-1))]}.matched.fq" paired reads to "${_bowtie_index[$((i))]}" with segment size between "$_min_size"bp and "$_max_size"bp using ("$_number_of_processors") processors:"
		if [ $_include_sam_output = "yes" ];
		then
		    # Align in sam format
            bowtie \
            -n $_maximum_number_of_mismatches_allowed \
            -m $_maximum_number_of_reportable_alignments \
            -p $_number_of_processors \
            -t \
            --sam \
            --tryhard \
            --chunkmbs 2048 \
            --phred$_quality_encoding_base-quals \
            --minins $_min_size \
            --maxins $_max_size \
            --al $_bowtied_output_files_path/aligned/${file%%_1.*}.fq \
            --un $_bowtied_output_files_path/nonaligned/${file%%_1.*}.fq \
            -1 $_working_files_path/${_read1%.*.*}.${_bowtie_index[$((i-1))]}.matched.fq \
            -2 $_working_files_path/${_read2%.*.*}.${_bowtie_index[$((i-1))]}.matched.fq \
            $_bowtie_index_path/${_bowtie_index[$((i))]} \
            $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.sam

            # Move mapped reads to aligned folder
            if [ -r "$_bowtied_output_files_path/aligned/${file%%_1.*}_1.fq" ];
            then
                mv $_bowtied_output_files_path/aligned/${file%%_1.*}_1.fq $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.paired.fq
                mv $_bowtied_output_files_path/aligned/${file%%_1.*}_2.fq $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.paired.fq
            fi
            if [ -r "$_bowtied_output_files_path/nonaligned/${file%%_1.*}_1.fq" ];
            then
                mv $_bowtied_output_files_path/nonaligned/${file%%_1.*}_1.fq $_bowtied_output_files_path/nonaligned/${file%%_1.*}_1.paired.fq
                mv $_bowtied_output_files_path/nonaligned/${file%%_1.*}_2.fq $_bowtied_output_files_path/nonaligned/${file%%_1.*}_2.paired.fq
            fi

            ####### Addition #######
            # Covert sam to bam, sort and index
			# count all reads in map
			_total_paired=($(samtools view $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.sam | wc -l))					
			# Filter out unmapped reads	
			echo -ne $(date +%H:%M:%S) "Converting sam to bam ..."		
			samtools view -bT $_path_to_reference${_bowtie_index[$((i))]}.fa $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.sam > $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.all.bam
			echo " done."			
			echo -ne $(date +%H:%M:%S) "Filtering bam file..."			
			samtools view -b -F 4 $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.all.bam	> $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.bam
			rm $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.all.bam						
			rm $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.sam
            samtools sort -T $_bowtied_output_files_path/tmp/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.sorted -o $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.sorted.bam $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.bam
            mv $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.sorted.bam $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.bam
			samtools index $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.bam
			_total_paired_aligned=($(samtools view -F 0x4 $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.bam | wc -l))
			echo "done."			
			#########################
        else
            bowtie \
            -n $_maximum_number_of_mismatches_allowed \
            -m $_maximum_number_of_reportable_alignments \
            -p $_number_of_processors \
            -t \
            --tryhard \
            --chunkmbs 2048 \
            --phred$_quality_encoding_base-quals \
            --minins $_min_size \
            --maxins $_max_size \
            --al $_bowtied_output_files_path/aligned/${file%%_1.*}.fq \
            --un $_bowtied_output_files_path/nonaligned/${file%%_1.*}.fq \
            -1 $_working_files_path/${_read1%.*.*}.${_bowtie_index[$((i-1))]}.matched.fq \
            -2 $_working_files_path/${_read2%.*.*}.${_bowtie_index[$((i-1))]}.matched.fq \
            $_bowtie_index_path/${_bowtie_index[$((i))]} \
            $_bowtied_output_files_path/aligned/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.bowtie
            if [ -r "$_bowtied_output_files_path/aligned/${file%%_1.*}_1.fq" ];
            then
                mv $_bowtied_output_files_path/aligned/${file%%_1.*}_1.fq $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.paired.fq
                mv $_bowtied_output_files_path/aligned/${file%%_1.*}_2.fq $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.paired.fq
            fi
            if [ -r "$_bowtied_output_files_path/nonaligned/${file%%_1.*}_1.fq" ];
            then
                mv $_bowtied_output_files_path/nonaligned/${file%%_1.*}_1.fq $_bowtied_output_files_path/nonaligned/${file%%_1.*}_1.paired.fq
                mv $_bowtied_output_files_path/nonaligned/${file%%_1.*}_2.fq $_bowtied_output_files_path/nonaligned/${file%%_1.*}_2.paired.fq
            fi
        fi

		#Aligning singles...
		#read -p "Press [Enter] key to continue ..."
	
		#Combining orphans and unaligned paired reads into a single file...
		echo -ne $(date +%H:%M:%S) "Combining R1 orphans with nonaligned R1 matched sequences into "${file%%_1.*}_1..${_bowtie_index[$((i-1))]}.singles.fq"  ... "
		cat \
		$( if [ -r "$_working_files_path/nonaligned/${file%%_1.*}_1.paired.fq" ]; then echo $_working_files_path/nonaligned/${file%%_1.*}_1.paired.fq; fi) \
		$( if [ -r "$_working_files_path/${file%%_1.*}_1.${_bowtie_index[$((i-1))]}.orphan.fq" ]; then echo $_working_files_path/${file%%_1.*}_1.${_bowtie_index[$((i-1))]}.orphan.fq; fi) \
		> $_working_files_path/${file%%_1.*}_1.${_bowtie_index[$((i-1))]}.singles.fq
		echo " done."

		#Aligning combined orphans and unaligned paired reads as single ended reads...
		#read -p "Press [Enter] key to continue ..."
		
		echo $(date +%H:%M:%S) "Aligning "${file%%_1.*}_1.${_bowtie_index[$((i-1))]}.singles.fq" reads to "${_bowtie_index[$((i))]}" using ("$_number_of_processors") processors:"
		if [ $_include_sam_output = "yes" ];
		then
			# Align in sam format
		    bowtie \
            -n $_maximum_number_of_mismatches_allowed \
            -m $_maximum_number_of_reportable_alignments \
            -p $_number_of_processors \
            -t \
            --sam \
            --tryhard \
            --chunkmbs 2048 \
            --phred$_quality_encoding_base-quals \
            --al $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.fq \
            --un $_bowtied_output_files_path/nonaligned/${file%%_1.*}_1.singles.fq \
            $_bowtie_index_path/${_bowtie_index[$((i))]} \
            $_working_files_path/${file%%_1.*}_1.${_bowtie_index[$((i-1))]}.singles.fq \
            $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.sam

			####### Addition #######
            # Covert sam to bam, sort and index
			# count all reads in map
			echo -ne "Converting sam to sorted, indexed bam file..."
			_total_R1=($(samtools view $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.sam | wc -l))					
			# Filter out unmapped reads			
			samtools view -bT $_path_to_reference${_bowtie_index[$((i))]}.fa $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.sam > $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.all.bam
			echo "done."		
			echo -ne $(date +%H:%M:%S) "Filtering bam ..."		
			samtools view -b -F 4 $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.all.bam	> $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.bam
			rm $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.all.bam			
			rm $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.sam
            samtools sort -T $_bowtied_output_files_path/tmp/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.sorted -o $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.sorted.bam $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.bam
            mv $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.sorted.bam $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.bam
			samtools index $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.bam
			_total_R1_aligned=($(samtools view -F 0x4 $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.bam | wc -l))
			echo "done."
			#########################

		else
            bowtie \
            -n $_maximum_number_of_mismatches_allowed \
            -m $_maximum_number_of_reportable_alignments \
            -p $_number_of_processors \
            -t \
            --tryhard \
            --chunkmbs 2048 \
            --phred$_quality_encoding_base-quals \
            --al $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.fq \
            --un $_bowtied_output_files_path/nonaligned/${file%%_1.*}_1.singles.fq \
            $_bowtie_index_path/${_bowtie_index[$((i))]} \
            $_working_files_path/${file%%_1.*}_1.${_bowtie_index[$((i-1))]}.singles.fq \
            $_bowtied_output_files_path/aligned/${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.bowtie
        fi

		echo -ne $(date +%H:%M:%S) "Combining R2 orphans with nonaligned R2 matched sequences into "${file%%_1.*}_2.${_bowtie_index[$((i-1))]}.singles.fq"  ... "
		cat \
		$( if [ -r "$_working_files_path/nonaligned/${file%%_1.*}_2.paired.fq" ]; then echo $_working_files_path/nonaligned/${file%%_1.*}_2.paired.fq; fi) \
		$( if [ -r "$_working_files_path/${file%%_1.*}_2.${_bowtie_index[$((i-1))]}.orphan.fq" ]; then echo $_working_files_path/${file%%_1.*}_2.${_bowtie_index[$((i-1))]}.orphan.fq; fi) \
		> $_working_files_path/${file%%_1.*}_2.${_bowtie_index[$((i-1))]}.singles.fq
		echo " done."

		echo $(date +%H:%M:%S) "Aligning "${file%%_1.*}_2.${_bowtie_index[$((i-1))]}.singles.fq" reads to "${_bowtie_index[$((i))]}" using ("$_number_of_processors") processors:"
		if [ $_include_sam_output = "yes" ];
		then
			# Align in sam format
		    bowtie \
            -n $_maximum_number_of_mismatches_allowed \
            -m $_maximum_number_of_reportable_alignments \
            -p $_number_of_processors \
            -t \
            --sam \
            --tryhard \
            --chunkmbs 2048 \
            --phred$_quality_encoding_base-quals \
            --al $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.fq \
            --un $_bowtied_output_files_path/nonaligned/${file%%_1.*}_2.singles.fq \
            $_bowtie_index_path/${_bowtie_index[$((i))]} \
            $_working_files_path/${file%%_1.*}_2.${_bowtie_index[$((i-1))]}.singles.fq \
            $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.sam

			####### Addition #######
            # Covert sam to bam, sort and index
			# count all reads in map
			_total_R2=($(samtools view $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.sam | wc -l))					
			# Filter out unmapped reads	
			echo -ne $(date +%H:%M:%S) "Converting sam to bam ..."			
			samtools view -bT $_path_to_reference${_bowtie_index[$((i))]}.fa $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.sam > $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.all.bam
			echo "done."			
			echo -ne $(date +%H:%M:%S) "Filtering bam ..."				
			samtools view -b -F 4 $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.all.bam	> $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.bam
			rm $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.all.bam	
			rm $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.sam
            samtools sort -T $_bowtied_output_files_path/tmp/${file%%_1.*}.${_bowtie_index[$((i))]}.paired.sorted -o $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.sorted.bam $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.bam
            mv $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.sorted.bam $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.bam
			samtools index $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.bam
			_total_R2_aligned=($(samtools view -F 0x4 $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.bam | wc -l))
			echo "done."	
			echo -e ${_bowtie_index[$((i))]}"\t"${file%%_1.*}"\t"$_total_paired"\t"$_total_paired_aligned"\t"$_total_R1"\t"$_total_R1_aligned"\t"$_total_R2"\t"$_total_R2_aligned >> $_stats_file	
			#########################     

		else
            bowtie \
            -n $_maximum_number_of_mismatches_allowed \
            -m $_maximum_number_of_reportable_alignments \
            -p $_number_of_processors \
            -t \
            --tryhard \
            --chunkmbs 2048 \
            --phred$_quality_encoding_base-quals \
            --al $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.fq \
            --un $_bowtied_output_files_path/nonaligned/${file%%_1.*}_2.singles.fq \
            $_bowtie_index_path/${_bowtie_index[$((i))]} \
            $_working_files_path/${file%%_1.*}_2.${_bowtie_index[$((i-1))]}.singles.fq \
            $_bowtied_output_files_path/aligned/${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.bowtie
        fi
		#Gzip all files to storage and cleanup
		#read -p "Press [Enter] key to continue ..."
		echo $(date +%H:%M:%S) "Compressing result files to storage:"	
	
		if [ ! -d "$_current_output_files_path" ]; then echo $(date +%H:%M:%S) "Creating: "$_current_output_files_path; mkdir $_current_output_files_path; fi

		cd $_bowtied_output_files_path/nonaligned
		for file_to_gzip in *.singles.fq;
		do 
			if [ -r "$file_to_gzip" ];
			then
				if [ ! -d "$_current_output_files_path/nonaligned" ]; then echo $(date +%H:%M:%S) "Creating: "$_current_output_files_path/nonaligned; mkdir $_current_output_files_path/nonaligned; fi
				echo -ne $(date +%H:%M:%S) "Compressing nonaligned "$file_to_gzip" into "${file_to_gzip%%.*}.fq.gz"  ... "
				gzip -f -c -9 -p $_number_of_processors $file_to_gzip > $_current_output_files_path/nonaligned/${file_to_gzip%%.*}.fq.gz
				rm $file_to_gzip
				echo " done."
			fi
		done

		#### Addition: step to carry out merging of aligned files ########

		if [ $_merge_aligned=="Y" ];
		    then
		    cd $_bowtied_output_files_path/aligned
			if [ $_include_sam_output = "yes" ];
			then
			    echo "Merging bam files..."
				samtools merge ${file%%_1.*}.${_bowtie_index[$((i))]}.merged.bam ${file%%_1.*}.${_bowtie_index[$((i))]}.paired.bam ${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.bam ${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.bam
				samtools sort -T $_bowtied_output_files_path/tmp/${file%%_1.*}.${_bowtie_index[$((i))]}.merged.sorted -o ${file%%_1.*}.${_bowtie_index[$((i))]}.merged.sorted.bam ${file%%_1.*}.${_bowtie_index[$((i))]}.merged.bam
				mv ${file%%_1.*}.${_bowtie_index[$((i))]}.merged.sorted.bam ${file%%_1.*}.${_bowtie_index[$((i))]}.merged.bam
				samtools index ${file%%_1.*}.${_bowtie_index[$((i))]}.merged.bam
			else
			    echo "Merging bowtie files..."
			    cat ${file%%_1.*}.${_bowtie_index[$((i))]}.paired.bowtie ${file%%_1.*}_1.${_bowtie_index[$((i))]}.singles.bowtie ${file%%_1.*}_2.${_bowtie_index[$((i))]}.singles.bowtie > $_bowtied_output_files_path/aligned ${file%%_1.*}.${_bowtie_index[$((i))]}.merged.bowtie
			fi
			rm *.singles.*
			rm *.paired.*
			echo "done."
		fi
        #################################################################

		### move bam and bai files to aligned path
		echo -ne $(date +%H:%M:%S) "Moving bam & bai files into "$_current_output_files_path/aligned" ..."		
		cd $_bowtied_output_files_path/aligned
		if [ ! -d "$_current_output_files_path/aligned" ]; then echo $(date +%H:%M:%S) "Creating: "$_current_output_files_path/aligned; mkdir $_current_output_files_path/aligned; fi
		mv *.{bam,bai} $_current_output_files_path/aligned
		echo "done."
		
		### gzip files if files outputted as bowtie files
		cd $_bowtied_output_files_path/aligned
		for file_to_gzip in *.*;
		do 
			if [ -r "$file_to_gzip" ];
			then
				if [ ! -d "$_current_output_files_path/aligned" ]; then echo $(date +%H:%M:%S) "Creating: "$_current_output_files_path/aligned; mkdir $_current_output_files_path/aligned; fi
				echo -ne $(date +%H:%M:%S) "Compressing aligned "$file_to_gzip" into "$file_to_gzip.gz"  ... "
				gzip -f -c -9 -p $_number_of_processors $file_to_gzip > $_current_output_files_path/aligned/$file_to_gzip.gz
				rm $file_to_gzip
				echo " done."
			fi
		done

		cd $_bowtied_output_files_path
		for file_to_gzip in *.matched.*;
		do 
			if [ -r "$file_to_gzip" ];
			then
				if [ ! -d "$_current_output_files_path/matched" ]; then echo $(date +%H:%M:%S) "Creating: "$_current_output_files_path/matched; mkdir $_current_output_files_path/matched; fi
				echo -ne $(date +%H:%M:%S) "Compressing matched "$file_to_gzip" into "$file_to_gzip.gz"  ... "
				gzip -f -c -9 -p $_number_of_processors $file_to_gzip > $_current_output_files_path/matched/$file_to_gzip.gz
				rm $file_to_gzip
				echo " done."
			fi
		done

		cd $_bowtied_output_files_path
		for file_to_gzip in *.orphan.*;
		do 
			if [ -r "$file_to_gzip" ];
			then
				if [ ! -d "$_current_output_files_path/orphaned" ]; then echo $(date +%H:%M:%S) "Creating: "$_current_output_files_path/orphaned; mkdir $_current_output_files_path/orphaned; fi
				echo -ne $(date +%H:%M:%S) "Compressing orphaned "$file_to_gzip" into "$file_to_gzip.gz"  ... "
				gzip -f -c -9 -p $_number_of_processors $file_to_gzip > $_current_output_files_path/orphaned/$file_to_gzip.gz
				rm $file_to_gzip
				echo " done."
			fi
		done
		
		rm -r *
		echo
	done

#exit
done

########## Addition: Perform male diagnosis #########################
echo -ne $(date +%H:%M:%S) "Begin Male/Female diagnosis ..."
echo "Input file path "$_input_files_path/$_start_time/${_bowtie_index[$((0))]}
if [ $_sex_diag=="Y" ];
    then
    if [ $_single_pop = "N" ];
        then
        $_scripts_file_path./MaleDiag.py $_input_files_path/$_start_time/${_bowtie_index[$((0))]} $_input_file_path/$_pop_map_name $_number_of_processors
    else
        $_scripts_file_path./MaleDiag.py $_input_files_path/$_start_time/${_bowtie_index[$((0))]} --single $_single_pop_name $_number_of_processors
    fi
	if [ -d "$_input_files_path/$_start_time/${_bowtie_index[$((0))]}/MaleDiag" ];
		then
		mv $_input_files_path/$_start_time/${_bowtie_index[$((0))]}/MaleDiag $_input_files_path/$_start_time
	fi
fi

echo "done."
####################################################################
#exit
date
echo "All done!"
} } 2>&1 | tee $_log_file_path/$_data_set"_"$_start_time"_trim_and_align.log"
