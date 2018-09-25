from bokeh.plotting import figure, show
from Bio import SeqIO
import os
import subprocess
from command_line_aligner import *

def show_read_length_fasta(file_name):
    buckets = []
    x = []

    for i in range(0, 60000, 500):
        x.append(i)
        buckets.append(0)
    
    for record in SeqIO.parse(file_name, "fasta"):
        read_length = int(len(record) / 500)
        if read_length >= len(buckets):
            print("read of length", read_length*500,"nt")
        else:
            buckets[read_length] += 1

    plot = figure()

    plot.vbar(x, 500, buckets)

    show(plot)

def split_fasta_q_file(file_name, num_sequences=100000):
    if file_name.endswith(".fasta"):
        cmd = "awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%"+ str(num_sequences) +"==0){file=sprintf(\""
        cmd += file_name
        cmd += ".part.%d.fasta\",n_seq);}"
        cmd += "print >> file; n_seq++; next;} { print >> file; }' <"
        cmd += file_name
        os.system(cmd)
    
    elif file_name.endswith(".fastq"):
        cmd = "awk 'BEGIN {n_seq=0;} /^@/ {if(n_seq%"+ str(num_sequences) +"==0){file=sprintf(\""
        cmd += file_name
        cmd += ".part.%d.fastq\",n_seq);}"
        cmd += "print >> file; n_seq++; next;} { print >> file; }' <"
        cmd += file_name
        os.system(cmd)
    else:
        print("unknown file ending.")

def calc_mem_usage(cmd_str, out_file):
    pre = "valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out " 
    post = " > /dev/null 2>&1"
    mem_usage_cmd = pre + cmd_str + post
    os.system(mem_usage_cmd)
    os.system("grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\\1/' | sort -g | tail -n 1 > " + out_file)

def get_mem_usages(fasta_file_long, fasta_file_short, reference, ref_fasta):
    ngmlr_cmd = "~/workspace/ngmlr/bin/ngmlr-0.2.8/ngmlr -r " + ref_fasta
    ngmlr_cmd += " -q " + fasta_file_long + "> /dev/null"
    calc_mem_usage(ngmlr_cmd, "max_mem/ngmlr.txt")
    return

    ma_cmd = "~/workspace/aligner/ma -t 1 -m pacBio -x " + reference + " -i " + fasta_file_long 
    calc_mem_usage(ma_cmd, "max_mem/MA.txt")

    mm_cmd = "~/workspace/minimap2/minimap2 --MD -c -a " + reference + ".mmi " + fasta_file_long 
    calc_mem_usage(mm_cmd, "max_mem/minimap.txt")

    bwa_cmd = "~/workspace/bwa/bwa mem " + reference + "bwa " + fasta_file_short 
    calc_mem_usage(bwa_cmd, "max_mem/bwa.txt")

    bowtie_cmd = "~/workspace/bowtie2/bowtie2-2.3.3.1/bowtie2 -x " 
    bowtie_cmd += reference + "bowtie2 -f -U " + fasta_file_short 
    calc_mem_usage(bowtie_cmd, "max_mem/bowtie.txt")

    g_map_cmd = "~/workspace/graphmap/bin/Linux-x64/graphmap align -v 0 -r " + ref_fasta 
    g_map_cmd += " -d " + fasta_file_long + " -o /dev/null"
    calc_mem_usage(g_map_cmd, "max_mem/graph_map.txt")

    blasr_cmd = "~/workspace/blasr/build/bin/blasr " + fasta_file_long + " " + ref_fasta
    blasr_cmd += " --sa " + reference + "blasr --sam --out /dev/null --cigarUseSeqMatch"
    calc_mem_usage(blasr_cmd, "max_mem/blasr.txt")

def start_up_times():
    reference = "/MAdata/genome/GRCh38.p12"
    ref_fasta = "/MAdata/chrom/human/GCA_000001405.27_GRCh38.p12_genomic.fna"
    aligners1 = [
            Ngmlr(ref_fasta, threads=32),
            #MA(reference, fast="pacBio", threads=32, silent=True),
            #Minimap2(reference, presetting="map-pb", threads=32, silent=True),
            #G_MAP(ref_fasta, threads=32),
            #Bowtie2(reference, threads=32),
            #BWA_MEM(reference, threads=32),
            #Blasr(reference, ref_fasta, threads=32),
        ]


    for aligner in aligners1[0:1]:
        print(aligner.get_name(), aligner.get_start_up_time(), "s")

start_up_times()

get_mem_usages(
        "/mnt/hdd0/giab/HG002/subreads/m150210_062505_42163R_c100779642550000001823165208251542_s1_p0.2.subreads.fasta",
        "/mnt/hdd0/giab/HG002/NIST_HiSeq_HG002_Homogeneity-10953946/paired-reads/2A1_CGATGT_L001_R1_002.fastq",
        "/MAdata/genome/GRCh38.p12",
        "/MAdata/chrom/human/GCA_000001405.27_GRCh38.p12_genomic.fna"
    )

