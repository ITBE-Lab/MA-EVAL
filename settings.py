# folder temporary files shall be stored in (best use a fast drive)
TEMP_FILE_FOLDER = "/mnt/ssd0/temp/"

# folder that contains all the git projects of the analyzed aligners
WORK_SPACE_FOLDER = "~/workspace/"

# the reference prefix in packed form
PACK_PREFIX = "/MAdata/genome/GRCh38.p12"
# the reference as fasta file
REFERENCE_FASTA = "/MAdata/chrom/human/GCA_000001405.27_GRCh38.p12_genomic.fna"

# fastaq file for ultralong oxford nanopore reads
UON_READS = "/mnt/hdd0/giab/HG002/Ultralong_OxfordNanopore/combined_2018-08-10.part.0.fastq"
# fasta file for pacBio reads
PAC_BIO_READS = "/mnt/hdd0/giab/HG002/subreads/*.fasta"
# fasta file for illumina reads
ILLUMINA_READS = "/mnt/hdd0/giab/HG002/NIST_HiSeq_HG002_Homogeneity-10953946/paired-reads/*.fastq"

# in SV analysis
# filter out all the small indels that are created due to the read simulation...
MIN_SV_INDEL_LEN = 25 