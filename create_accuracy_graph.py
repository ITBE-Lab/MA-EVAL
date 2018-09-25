from check_accuracy import *
from read_simulation import *
import glob
from load_contigs import load_contigs

def create_accuracy_graph(
        file_regex, ref, ref_fasta, out_prefix, 
        limit_sample_files_to=10, end=350, step=10, aligners=[], sampler="bwa",
        single_core=False, num_reads=1000, with_secondary=True
    ):
    files_list = []
    for file in glob.glob(file_regex)[:limit_sample_files_to]:
        files_list.append(file)

    # load all Chromosomes
    contigs = load_contigs(ref_fasta)

    read_simulator = ReadSimulator(ref, ref_fasta, contigs)
    if not read_simulator.load_sampled_from_file(out_prefix + "_sampled_distrib"):
        read_simulator.sample_distrib_from_fasta(fasta_files=files_list, pick=sampler)
        read_simulator.save_sampled_to_file( out_prefix + "_sampled_distrib")

    if len(aligners) > 0:
        read_simulator.accuracy_coverage_graph(
                out_file_name=out_prefix+"_acc_cov", end=end, step=step, aligners=aligners,
                single_core=single_core, num_reads=num_reads, allow_secondary=with_secondary
            )
    read_simulator.read_distribution_graph(
            out_file_name=out_prefix + "_read_distrib",
            sample_files=files_list,
            sampler=sampler
        )

create_accuracy_graph(
    PAC_BIO_READS,
    PACK_PREFIX,
    REFERENCE_FASTA,
    "pacBio",
    sampler="mm",
    aligners=[
        Ngmlr(REFERENCE_FASTA, threads=32),
        Blasr(PACK_PREFIX, REFERENCE_FASTA, threads=32),
        MA(PACK_PREFIX, fast="pacBio", threads=32),
        Minimap2(PACK_PREFIX, presetting="map-pb", threads=32),
        G_MAP(REFERENCE_FASTA, threads=32),
    ]
)

create_accuracy_graph(
    ILLUMINA_READS,
    PACK_PREFIX,
    REFERENCE_FASTA,
    "Illumina",
    sampler="bwa",
    end=100,
    num_reads=100000,
    aligners=[
        Bowtie2(PACK_PREFIX, threads=32),
        MA(PACK_PREFIX, fast="acc", threads=32),
        BWA_MEM(PACK_PREFIX, threads=32),
    ]
)

create_accuracy_graph(
    UON_READS,
    PACK_PREFIX,
    REFERENCE_FASTA,
    "oxfNano",
    sampler="mm-ont",
    limit_sample_files_to=5,
    aligners=[
        Blasr(PACK_PREFIX, REFERENCE_FASTA, threads=32),
        MA(PACK_PREFIX, fast="pacBio", threads=32),
        Minimap2(PACK_PREFIX, presetting="map-ont", threads=32),
        G_MAP(REFERENCE_FASTA, threads=32),
        Ngmlr(REFERENCE_FASTA, threads=32),
    ]
)
