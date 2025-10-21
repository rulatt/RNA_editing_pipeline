# RNA editing pipeline




<pre> bash singularity exec snakemake_9.12.0.sif snakemake -s REDItools2_pipeline.sf --cores 4  </pre>



<pre> 

  


</pre>




'config_DNAseq.yaml' 

<pre> 

software_dir: "/path/to/software_dir"  # Path to .sif files and Rscript 
resource_dir: "/path/to/ressource" # Path to ressources files (GATK ressources files)
output_dir: "/path/to/output_dir" # output directory 

threads: "40"  # Number of threads to run BWA 

reads_R1: "/path/to/reads_R1.fastq.gz" # Full path to R1 & R2 fastq files
reads_R2: "/path/to/reads_R2.fastq.gz"

ref_fa: "/path/to/hg38.genome.fa" # full path to reference genome 

DNA_ID: "DNA_ID_WGS" # DNA_ID, prefix of output .sam & .bam files

reads_LB: "unknown"  # Read group infos 
reads_PL: "ILLUMINA"  # Read group infos 
  
</pre>



Config file name = config_editing.yaml, can be overriden with --configfile
