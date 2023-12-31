#!/bin/bash

# requiers pigz !
SRR=$1
GenomRef=$2
threads=16

# incase you have a big cluster make sure to load the conda module
module load tools/conda/23.1.0
source /share/apps/tools/python/3.8/etc/profile.d/conda.sh

#otherwise you have to load the conda config, here are some examples:
#source /home/graujh/miniconda3/etc/profile.d/conda.sh
#source /home/graujh/miniconda3/pkgs/conda-22.9.0-py38h06a4308_0/etc/profile.d/conda.sh



# check if genome exists, otherwise exit #######################################
check_file_exists() {
    local file="$1"

    if [ ! -f "$file" ]; then
        echo "ERROR -> Genome reference $file does not exist. Exiting..."
        exit 1
    fi
}

check_file_exists $GenomRef

# Check and install conda packages #############################################

check_and_install_conda_package() {
    local env_name="$1"
    local package="$2"
    local version="$3"
    local channel="$4"

    # Check if conda environment exists
    if ! conda env list | grep -q "$env_name"; then
        echo "Creating a new conda environment: $env_name"
        conda create -n "$env_name" -y
    fi

    # Activate the conda environment
    echo "Activating conda environment: $env_name"
    conda activate "$env_name"

    # Check if the package is installed in the conda environment
    if ! conda list -n "$env_name" | grep -q "^$package"; then
        echo "$package is not installed in the conda environment. Installing version $version using conda from channel: $channel"
        conda install -n "$env_name" -y -c "$channel" "$package==$version"
    else
        local installed_version=$(conda list -n "$env_name" | awk -v pkg="$package" '$1 == pkg { print $2 }')
        if [[ "$installed_version" == "$version" ]]; then
            echo "$package version $version is already installed in the conda environment."
        else
            echo "Updating $package to version $version using conda from channel: $channel"
            conda install -n "$env_name" -y -c "$channel" "$package==$version"
        fi
    fi
}

install_necessary() {
    check_and_install_conda_package "bwaMP" "bwa-mem2" "2.2.1" "bioconda"
    check_and_install_conda_package "bwaMP" "samtools" "1.17" "bioconda"
    check_and_install_conda_package "sra-toolsMP" "sra-tools" "3.0.3" "bioconda"
    check_and_install_conda_package "fastpMP" "fastp" "0.23.4" "bioconda"
    check_and_install_conda_package "bcftoolsMP" "bcftools" "1.17" "bioconda"
    check_and_install_conda_package "mosdepthMP" "mosdepth" "0.3.3" "bioconda"
    check_and_install_conda_package "mahajrodMP" "mavr" "0.97" "mahajrod"
    check_and_install_conda_package "bedtoolsMP" "bedtools" "2.31.0" "bioconda" 
}

install_necessary


# Check dependencies flag ###################################################
if [[ "$1" == "--check_dependencies" ]]; then
	echo "checking dependencies and exiting thereafter"
	install_necessary
	exit 0
fi

# check if fastqfiles are not already there 
if [ ! -f "$SRR.R2_clipped.fastq.gz" ]; then
        
	#### install FASTERQDUMP #######################################################
	# Specify the installation directory
	install_dir="$PWD"

	# Check if fasterq-dump is available
	if command -v fasterq-dump >/dev/null 2>&1; then
	    echo "fasterq-dump is already installed."
	else
	    echo "fasterq-dump is not found. Downloading and installing..."

	    # remove previos setup and download the latest fasterq-dump binary from NCBI
	    rm -rf sratool* && wget https://ftp.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

	    # Extract the contents of the tarball
	    tar -xf sratoolkit.current-ubuntu64.tar.gz && rm sratoolkit.current-ubuntu64.tar.gz

	    # Move to new file name and delete old one
	    mv sratoolkit.* sratoolkit.current-ubuntu64

	    echo "fasterq-dump is installed in $install_dir."
	fi

	# Set the path to fasterq-dump
	fasterq_dump="./sratoolkit.current-ubuntu64/bin/fasterq-dump"

	# Use the fasterq-dump path in your subsequent commands
	echo "Using fasterq-dump from $fasterq_dump_path"



	## Download and convert to FASTQ ################################################
	$fasterq_dump $SRR --outdir . --threads $threads --split-3 --verbose
	echo "finished fasterq-dump"

	
	# FASTP adapter and quality trimming ###########################################
	conda activate fastpMP
	
	fastp --in1 $SRR'_1.fastq' --in2 $SRR'_2.fastq' \
	 --out1 $SRR.R1_clipped.fastq.gz --out2 $SRR.R2_clipped.fastq.gz \
	 --qualified_quality_phred 15 --unqualified_percent_limit 40 --average_qual 30 -l 50 \
	 -h $SRR'.fastp.html' -j $SRR'.fastp.json' --correction --thread $threads 1>$SRR.fastp.log 2>$SRR.fastp.err \
	 && rm  $SRR'_1.fastq' $SRR'_2.fastq'
	 
	conda deactivate
	echo "finished fastp quality and adapter removal"
	

fi

#### BWA alignment #############################################################
conda activate bwaMP

bwa-mem2 mem -t $threads -M -R '@RG\tID:sample\tLB:sample\tPL:ILLUMINA\tSM:sample' \
	$GenomRef \
	$SRR.R1_clipped.fastq.gz $SRR.R2_clipped.fastq.gz 2>bwa.log | \
	samtools sort -@12 -n -o /dev/stdout | \
	samtools fixmate - /dev/stdout | \
	samtools sort -@ 16 -o /dev/stdout | \
	samtools rmdup -S - $SRR.onref_aligned_sorted_fix_sorted_rmdup.bam \
	&& samtools index $SRR.onref_aligned_sorted_fix_sorted_rmdup.bam

conda deactivate
echo "finished bwa+samtools"



#### bcftools ##################################################################
conda activate bcftoolsMP

bcftools mpileup \
--max-depth 250 --min-MQ 30 --min-BQ 30 \
--adjust-MQ 50 --threads $threads \
--annotate AD,INFO/AD,ADF,INFO/ADF,ADR,INFO/ADR,DP,SP,SCR,INFO/SCR \
--output-type u --fasta-ref $GenomRef $SRR.onref_aligned_sorted_fix_sorted_rmdup.bam \
--output $SRR.onref.bcf


bcftools call --group-samples - --multiallelic-caller \
--output-type u --variants-only \
--format-fields GQ,GP --threads $threads -o $SRR.onref.vcf $SRR.onref.bcf


bcftools filter -o $SRR.onref.filtered.vcf \
--exclude 'QUAL < 20.0 || (FORMAT/SP > 60.0 || FORMAT/DP < 5.0 || FORMAT/GQ < 20.0)' $SRR.onref.vcf

conda deactivate
echo "finished bcftools"



#### mosdepth ##################################################################
conda activate mosdepthMP

mosdepth -t 16 $SRR.mosdepth $SRR.onref_aligned_sorted_fix_sorted_rmdup.bam

conda deactivate
echo "finished mosdepth"


#### mahajrod ##################################################################
conda activate mahajrodMP

get_windows_stats_mosdepth_per_base_file.py -o $SRR.mahajrod -w 1000000 -s 100000 -b 100000000 -i $SRR.mosdepth.per-base.bed.gz
generate_mask_from_coverage_bed.py -c $SRR.mosdepth.per-base.bed.gz -m 42 -o $SRR.mahajrod.masking.bed -x 2.5  -n 0.33

conda deactivate
echo "finished mahajrod"


#### bedtools ##################################################################
conda activate bedtoolsMP

bedtools intersect -v -a $SRR.onref.vcf -b $SRR.mahajrod.masking.bed -header >$SRR.mahajrod.masked.vcf

conda deactivate
echo "finished bedtools"


#### bcftools ##################################################################
conda activate bcftoolsMP

bcftools filter -i TYPE=\"snp\" -O v $SRR.mahajrod.masked.vcf >$SRR.mahajrod.masked.SNP.vcf
bcftools filter -i FMT/GT=\"het\" -O v $SRR.mahajrod.masked.SNP.vcf >$SRR.mahajrod.masked.SNP.hetero.vcf

conda deactivate
echo "finished bcftools"
