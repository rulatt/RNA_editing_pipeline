# RNA editing pipeline


## DNA-seq Pipeline 

Example to run full DNA seq pipeline

<pre>singularity exec snakemake_9.12.0.sif snakemake -s DNAseq_align_pipeline.sf --cores 4 all </pre>

Available commands : 

<pre>all                           # Run all commands
run_BWA_meme                  # Run BWA meme aligner (output : {DNA_ID}.sam)
run_GATK_SortSAM              # Sort SAM file & convert to BAM file (output : {DNA_ID}.sort.bam) 
run_samtools_process          # Keep only properly paired aligned reads with a MAPQ ≥ 20  (output : {DNA_ID}.processed.bam)
run_GATK_rmdup                # Remove duplicated reads (output : {DNA_ID}.rmdup.bam) 
run_GATK_Base_Recalibrator    # Compute base recalibration (output : {DNA_ID}_recal_data.table) 
run_GATK_Apply_Recalibration  # Apply base recalibration (output : {DNA_ID}.rmdup.bqsr.bam) 
</pre>

Final output {DNA_ID}.rmdup.bqsr.bam is used for RNA editing analysis 
 

### Config file 

A config file must be provided (default : 'config_DNAseq.yaml',  can be overriden with --configfile).

Example : 

<pre># config_DNAseq.yaml
# Configuration file for DNAseq Snakemake pipeline

software_dir: "/path/to/software_dir"  # Path to .sif files and Rscript 
resource_dir: "/path/to/ressource" # Path to ressources files (GATK ressources files)
output_dir: "/path/to/output_dir" # output directory 

threads: "40"  # Number of threads to run BWA meme

reads_R1: "/path/to/DNA_reads_R1.fastq.gz" # Full path to R1 & R2 fastq files
reads_R2: "/path/to/DNA_reads_R2.fastq.gz"

ref_fa: "hg38.genome.fa" # full path to reference genome 

DNA_ID: "DNA_ID_WGS" # DNA_ID, prefix of output .sam & .bam files

reads_LB: "unknown"  # Read group infos 
reads_PL: "ILLUMINA"  # Read group infos 

snps_1000G: "resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf"
known_indels: "resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf"
indels_1000G: "resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf"
</pre>


## RNA-seq Pipeline 

Example to run full RNA seq pipeline 

<pre> singularity exec snakemake_9.12.0.sif snakemake -s RNAseq_align_pipeline.sf --cores 4 all</pre>


Available commands : 

<pre>all                           # Run all commands
run_STAR_genome_generate      # Generates genome index for STAR 
run_STAR_mapping              # Run STAR aligner (output: {RNA_ID}Aligned.sortedByCoord.out.bam)
run_samtools_process          # Keep only properly paired aligned reads with a MAPQ ≥ 20 (output : {RNA_ID}.processed.bam") 
run_GATK_splitN               # Formatting BAM file for next processing steps (output : {RNA_ID}.split.bam")    
run_GATK_rmdup                # Remove duplicated reads (output : {RNA_ID}.rmdup.bam)
run_GATK_Base_Recalibrator    # Compute base recalibration (output : {RNA_ID}_recal_data.table")
run_GATK_Apply_Recalibration  # Apply base recalibration (output : {RNA_ID}.rmdup.split.bqsr.bam ) 
</pre>

Final output {RNA_ID}.rmdup.split.bqsr.bam is used for RNA editing analysis 

### Config file 

A config file must be provided (default : 'config_RNAseq.yaml',  can be overriden with --configfile).

Example : 

<pre># config_RNAseq.yaml
# Configuration file for RNAseq Snakemake pipeline

software_dir: "/path/to/software_dir"
resource_dir: "/path/to/ressource"
output_dir: "/path/to/output_dir"  

genome_dir: "/path/to/STAR_index/" 

threads: "5"
gtf_file: "annotation_file.gtf"  # annotation file, must be located in resource_dir 

reads_R1: "/path/to/RNA_reads_R1.fastq.gz"
reads_R2: "/path/to/RNA_reads_R2.fastq.gz"

ref_fa: "hg38.genome.fa"

RNA_ID: "RNA_Tumor_ID" # prefix of output.bam files

snps_1000G: "resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf"
known_indels: "resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf"
indels_1000G: "resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf"
</pre>  


## Editing Pipeline

Example to run full RNA editing detection pipeline 

<pre> singularity exec snakemake_9.12.0.sif snakemake -s REDItools2_pipeline.sf --cores 4 all</pre>

Available commands : 

<pre>all                           # Run all commands 
run_RNA_detection             # Detect variants in RNA  (output : {RNA_ID}.table.txt )
run_BED_conversion            # Convert RNA variant table in BED table for DNA variant detection  (output : {RNA_ID}.table.bed )
run_DNA_detection             # Detect variants in DNA based on RNA variant positions (output : {DNA_ID}.table.txt )
run_annotation                # Annotate RNA variant table with DNA variant table (output : {RNA_ID}_annotated.table.txt )   
run_r_filtering               # Filter annotated RNA variant table (output : {RNA_ID}_annotated_filtered.table.txt )
</pre>

Final filtered and annotated table : {RNA_ID}_annotated_filtered.table.txt

### Config file 

A config file must be provided (default : 'config_editing.yaml',  can be overriden with --configfile).

Example : 

<pre># config_editing.yaml
# Configuration file for REDItools2.0 Snakemake pipeline

software_dir: "/path/to/software_dir"
resource_dir: "/path/to/ressource"
output_dir: "/path/to/output_dir"

RNA_bam: "/path/to/RNA_alignment_file.bam"  # Full path to RNA BAM file
DNA_bam: "/path/to/DNA_alignment_file.bam" # Full path to DNA BAM file

ref_fa: "hg38.genome.fa"

RNA_ID: "RNA_Tumor_ID" # RNA BAM file prefix
DNA_ID: "DNA_ID_WGS" # DNA BAM file prefix

custom_bed: ""                 # Custom BED file for run_DNA_detection step (if not running run_BED_conversion step) 

coverage_threshold: 10         # Minimum coverage at the edited site
frequency_threshold_min: 0.1   # Minimum Variant Allele Frequency
frequency_threshold_max: 0.95  # Maximum Variant Allele Frequency
KeepOneVar: true               # Keep only monoallelic RNA variant (true or false)
StrictFiltering: true          # Keep only variant with no change in DNA corresponding position  (true or false)
KeepEditing: true              # Keep only A->G & T->C variants (true or false)
</pre>   


Notes : 

You can run only parts of the scripts, for example : 

<pre>singularity exec snakemake_9.12.0.sif snakemake -s REDItools2_pipeline.sf --cores 4 run_annotation run_r_filtering </pre>
