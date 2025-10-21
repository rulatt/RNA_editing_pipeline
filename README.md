# RNA editing pipeline




<pre> bash singularity exec snakemake_9.12.0.sif snakemake -s REDItools2_pipeline.sf --cores 4  </pre>



<pre> 

  


</pre>




'config_DNAseq.yaml' 

<pre> 
# config_DNAseq.yaml
# Configuration file for DNAseq Snakemake pipeline

software_dir: "/path/to/software_dir"  # Path to .sif files and Rscript 
resource_dir: "/path/to/ressource" # Path to ressources files (GATK ressources files)
output_dir: "/path/to/output_dir" # output directory 

threads: "40"  # Number of threads to run BWA 

reads_R1: "/path/to/DNA_reads_R1.fastq.gz" # Full path to R1 & R2 fastq files
reads_R2: "/path/to/DNA_reads_R2.fastq.gz"

ref_fa: "/path/to/hg38.genome.fa" # full path to reference genome 

DNA_ID: "DNA_ID_WGS" # DNA_ID, prefix of output .sam & .bam files

reads_LB: "unknown"  # Read group infos 
reads_PL: "ILLUMINA"  # Read group infos 
  
</pre>


'config_RNAseq.yaml'

<pre> 

# config_RNAseq.yaml
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


'config_editing.yaml' 

<pre> 
  
# config_editing.yaml
# Configuration file for REDItools2.0 Snakemake pipeline

software_dir: "/path/to/software_dir"
resource_dir: "/path/to/ressource"
output_dir: "/mnt/iribhm/homes/rulattuc/P_Mathieu/Reditools_RES/A549_editing_RES/DLTR00064_rmdup_bqsr"

RNA_bam: "/mnt/iribhm/homes/rulattuc/P_Mathieu/Reditools_RES/ADAR/DLTR00064/DLTR00064.rmdup.split.bqsr.bam"
DNA_bam: "/mnt/iribhm/homes/rulattuc/P_Mathieu/Reditools_RES/A549_Genome_Mapped/A549_WGS.rmdup.bqsr.bam"

ref_fa: "Homo_sapiens_assembly38.fa"

RNA_ID: "DLTR00064_rmdup_bqsr"
DNA_ID: "A549_WGS"

custom_bed: ""

coverage_threshold: 10         # Minimum coverage at the edited site
frequency_threshold_min: 0.1   # Minimum Variant Allele Frequency
frequency_threshold_max: 0.95  # Maximum Variant Allele Frequency
KeepOneVar: true               # Keep only monoallelic RNA variant (true or false)
StrictFiltering: true          # Keep only variant with no change in DNA corresponding position  (true or false)
KeepEditing: true              # Keep only A->G & T->C variants (true or false)

</pre>   


Notes : 
Default config file names are defined in each snakemake file and can be  overriden with --configfile
