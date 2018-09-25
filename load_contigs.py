from Bio import SeqIO

def load_contigs(ref_fasta):
    contigs = []
    start = 0
    for contig in SeqIO.parse(ref_fasta, "fasta"):
        #if len(contigs) > 2:
        #    break
        name = contig.name
        length = len(contig)
        #merely load chromosomes
        if name[:2] != "CM":
            start += length
            continue
        print("contig:", name,":", start, "-", start+length)
        seq = str(contig.seq)
        assert len(seq) == length
        contigs.append( (name, seq) )
        start += length
    return contigs