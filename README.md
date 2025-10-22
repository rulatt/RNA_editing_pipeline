# RNA editing pipeline




<pre> bash singularity exec snakemake_9.12.0.sif snakemake -s REDItools2_pipeline.sf --cores 4  </pre>


Example to run full DNA seq pipeline

<pre> singularity exec snakemake_9.12.0.sif snakemake -s DNAseq_align_pipeline.sf --cores 4 run_DNAseq_pipeline </pre>

Available commands : 

<pre>  

run_DNAseq_pipeline           # Run full alignment pipeline
run_BWA_meme                  # Run BWA meme aligner
run_GATK_SortSAM              # Sort SAM file & convert to BAM file
run_samtools_process          # Keep only properly paired aligned reads with a quality of ##"####20 
run_GATK_rmdup                # remove duplicated reads
run_GATK_Base_Recalibrator    # Compute base recalibration
run_GATK_Apply_Recalibration  # Apply base recalibration
  
</pre>

Final *.rmdup.bqsr.bam BAM File is used for RNA editing analysis 
 

'config_DNAseq.yaml' 

<pre> # config_DNAseq.yaml
# Configuration file for DNAseq Snakemake pipeline

software_dir: "/path/to/software_dir"  # Path to .sif files and Rscript 
resource_dir: "/path/to/ressource" # Path to ressources files (GATK ressources files)
output_dir: "/path/to/output_dir" # output directory 

threads: "40"  # Number of threads to run BWA meme

reads_R1: "/path/to/DNA_reads_R1.fastq.gz" # Full path to R1 & R2 fastq files
reads_R2: "/path/to/DNA_reads_R2.fastq.gz"

ref_fa: "/path/to/hg38.genome.fa" # full path to reference genome 

DNA_ID: "DNA_ID_WGS" # DNA_ID, prefix of output .sam & .bam files

reads_LB: "unknown"  # Read group infos 
reads_PL: "ILLUMINA"  # Read group infos 
  
</pre>



Example to run full RNA seq pipeline 


<pre> singularity exec snakemake_9.12.0.sif snakemake -s RNAseq_align_pipeline.sf --cores 4 run_STAR_pipeline </pre>


Available commands : 

<pre>  
  
run_STAR_pipeline             # Run full alignment pipeline
run_STAR_genome_generate      # Generates genome index for STAR 
run_STAR_mapping              # Run STAR aligner 
run_samtools_process          # Keep only properly paired aligned reads with a quality of ##"####20
run_GATK_splitN               # Formatting BAM file for next processing steps    
run_GATK_rmdup                # remove duplicated reads
run_GATK_Base_Recalibrator    # Compute base recalibration
run_GATK_Apply_Recalibration  # Apply base recalibration
  
</pre>

Final *rmdup.split.bqsr.bam BAM file is used for RNA editing analysis 


'config_RNAseq.yaml'

<pre> # config_RNAseq.yaml
# Configuration file for RNAseq Snakemake pipeline

software_dir: "/path/to/software_dir"
resource_dir: "/path/to/ressource"
output_dir: "/path/to/output_dir"  

genome_dir: "/mnt/iribhm/genomes/hg38-gatk/STAR_hg38_gatk_index_GTF_TEST/"

threads: "5"
gtf_file: "annotation_file.gtf"  # annotation file, must be located in resource_dir 

reads_R1: "/path/to/RNA_reads_R1.fastq.gz"
reads_R2: "/path/to/RNA_reads_R2.fastq.gz"

ref_fa: "/path/to/hg38.genome.fa"

RNA_ID: "RNA_Tumor_ID" # prefix of output.bam files
  
</pre>  


Example to run full RNA editing detection pipeline 

<pre> singularity exec snakemake_9.12.0.sif snakemake -s REDItools2_pipeline.sf --cores 4 run_editing_pipeline </pre>

Available commands : 

<pre>  
  
run_editing_pipeline          # Run full editing pipeline
run_RNA_detection             # Detects variants in RNA  
run_BED_conversion            # Converts RNA variant table in BED table for DNA variant detection 
run_DNA_detection             # Detects variants in DNA based on RNA variant positions
run_annotation                # Annotate RNA variant table with DNA variant table    
run_r_filtering               # Filter annotated RNA variant table
 
</pre>


Final file for downstream analysis : 
Intermediate files : 


'config_editing.yaml' 

<pre> # config_editing.yaml
# Configuration file for REDItools2.0 Snakemake pipeline

software_dir: "/path/to/software_dir"
resource_dir: "/path/to/ressource"
output_dir: "/path/to/output_dir"

RNA_bam: "/path/to/RNA_alignment_file.bam"  # Full path to RNA BAM file
DNA_bam: "/path/to/DNA_alignment_file.bam" # Full path to DNA BAM file

ref_fa: "hg38.genome.fa"

RNA_ID: "RNA_Tumor_ID" # RNA BAM file prefix
DNA_ID: "DNA_ID_WGS" # DNA BAM file prefix

custom_bed: ""                 # Custom bed use for Reditools on DNA steps, none by default (generated after)

coverage_threshold: 10         # Minimum coverage at the edited site
frequency_threshold_min: 0.1   # Minimum Variant Allele Frequency
frequency_threshold_max: 0.95  # Maximum Variant Allele Frequency
KeepOneVar: true               # Keep only monoallelic RNA variant (true or false)
StrictFiltering: true          # Keep only variant with no change in DNA corresponding position  (true or false)
KeepEditing: true              # Keep only A->G & T->C variants (true or false)

</pre>   


Notes : 
Default config file names are defined in each snakemake file and can be  overriden with --configfile

You can Run only part of the scripts, for example : 
