#!/bin/sh
#$ -S /bin/sh
set -o errexit

#Analysis directory
AD=`pwd`

# load modules 
module load GATK/4.3.0.0

source ${AD}/${CONFIG_FILE} # Added source config file to upload environmental variables, in particular GENOME 

# Input dir
INPUT_DIR=09-variant-calling

# Output dir
OUTPUT_DIR=10-filtered_variants
if [[ ! -d ${OUTPUT_DIR} ]]; then
	mkdir -m 777 ${OUTPUT_DIR}
fi	


# Number of threads
NT=${SLURM_CPUS_PER_TASK}
MEM=${SLURM_MEM_PER_NODE}
TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`

# Print header
echo ""
echo "#################################"
echo "#" ${SAMPLE}
echo "#" MEMORY PER NODE: ${MEM}
echo "#" CPUS: ${NT}
echo "# Merging gvcf files start:"
echo "#" $TIMESTAMP
echo "#################################"
echo ""

# >>>> NEW CODE

# Search available g.vcf files to merge
MY_GVCF_FILES=()

for file in $(find ${INPUT_DIR} -type f -name "*.g.vcf.gz"); do # Mofified find command. Original command = find ${INPUT_DIR} -type f -name "*.vcf.gz" 
    MY_GVCF_FILES+="-V ${file} "
done

if [[ ! ${MY_GVCF_FILES[@]} ]]; then
    echo
    echo ERROR
    echo No g.vcf files found in ${OUTPUT_DIR}
    echo 
    exit 1
fi

#################################################################
###                     Combine all gvcfs                     ###
#################################################################

OUTPUT_GVCF=${OUTPUT_DIR}/all_samples.vcf.gz

gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' CombineGVCFs \
    -R ${GENOME} \
    ${MY_GVCF_FILES[@]} \
    -O ${OUTPUT_GVCF}

#################################################################
###    Call Variants within chromosomes across all samples    ###
#################################################################

# Print header
echo ""
echo "#################################"
echo "#" ${SAMPLE}
echo "#" MEMORY PER NODE: ${MEM}
echo "#" CPUS: ${NT}
echo "# Calling genotypes start:"
echo "#" $TIMESTAMP
echo "#################################"
echo ""


OUTPUT_GENOTYPED_VCF=${OUTPUT_DIR}/all_samples_genotyped.vcf.gz

gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" GenotypeGVCFs \
    -R ${GENOME} \
    -V ${OUTPUT_GENOTYPED_VCF} \
    --window 35 \
    --cluster 3 \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${FILTERED_SNP_VCF}

# Filter Indels
FILTERED_INDEL_VCF=${OUTPUT_DIR}/all_samples_INDEL_filtered.vcf.gz

gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantFiltration \
    -R ${GENOME} \
    -V ${OUTPUT_GENOTYPED_VCF} \
    --window 35 \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O ${FILTERED_INDEL_VCF}

# Print footnote
TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
echo ""
echo "#################################"
echo "# Finished filtering variants ..."
echo "# $TIMESTAMP"
echo "#################################"
echo ""
