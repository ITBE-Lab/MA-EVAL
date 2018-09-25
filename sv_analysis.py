from check_accuracy import *
import glob
from read_simulation import *
import pysam
from bokeh.core.properties import value
from bokeh.io import show, output_file
from bokeh.plotting import figure
from bokeh.models import FactorRange
from load_contigs import load_contigs
from settings import *


def alignments_to_coverage(alignment_list, contig_name):
    cov_list = []
    for index, alignment in enumerate(alignment_list):
        if alignment.contig_name != contig_name:
            continue
        start = alignment.start
        length = 0
        for amount, char in alignment.cigar:
            if char in ['M', '=', 'X']: # match / mismatch
                length += amount
            if char in ['S', 'H', 'P']: # soft & hard -clipping & padding
                pass
            if char in ['I']: # insertion
                if amount < MIN_SV_INDEL_LEN:
                    pass
                else:
                    cov_list.append( (start, start+length, index) )
                    start += length
                    length = 0
            if char in ['D']: # deletion
                if amount < MIN_SV_INDEL_LEN:
                    length += amount
                else:
                    cov_list.append( (start, start+length, index) )
                    start += length + amount
                    length = 0
        cov_list.append( (start, start+length, index) )
        assert start + length == alignment.start + alignment.length
    #print(cov_list)
    return cov_list

def close(a, b, max_dist=25):
    return abs(a-b) < max_dist

def exact(coverage_list, point, check_start=True, check_end=True):
    for start, end, index in coverage_list:
        if check_start and close(start, point):
            return index
        if check_end and close(end, point):
            return index
    return None

def approx(coverage_list, point, check_start=True, check_end=True):
    for start, end, index in coverage_list:
        if check_start and close(start, point, max_dist=100):
            return index
        if check_end and close(end, point, max_dist=100):
            return index
    return None

#
# there must be at least two fragments between start - 100 and end + 100
#
def forced_del(coverage_list, start, end):
    count = 0
    for start_c, end_c, index in coverage_list:
        if start_c < start - 10 and end_c > end + 10:
            return True
    return False

#
# exact hits must be present for start and end
# gap must be covered by the same alignment
#
def precise_del(coverage_list, start, end):
    a = exact(coverage_list, start)
    b = exact(coverage_list, end)
    if a is None or b is None:
        return False
    return a == b
#
# approx hits must be present for start and end
# gap must be covered by the same alignment
#
def indicated_del(coverage_list, start, end):
    a = approx(coverage_list, start)
    b = approx(coverage_list, end)
    if a is None or b is None:
        return False
    return a == b
#
# exact hits must be present for start and end
# gap can be covered by different alignments
#
def split_del(coverage_list, start, end):
    a = exact(coverage_list, start)
    b = exact(coverage_list, end)
    if a is None or b is None:
        return False
    return True
#
# exact hits must be present for start or end
#
def trimmed_del(coverage_list, start, end):
    a = approx(coverage_list, start)
    b = approx(coverage_list, end)
    if a is None and b is None:
        return False
    return True

#
# there must be at least two fragments between pos - 100 and pos + 100
#
def forced_ins(coverage_list, pos):
    for start_c, end_c, index in coverage_list:
        if start_c < pos - 100 and end_c > pos + 100:
            return True
    return False
#
# exact hits must be present for start and end
# gap must be covered by the same alignment
#
def precise_ins(coverage_list, pos):
    a = exact(coverage_list, pos, check_start=False, check_end=True)
    b = exact(coverage_list, pos, check_start=True, check_end=False)
    if a is None or b is None:
        return False
    return a == b
#
# approx hits must be present for start and end
# gap must be covered by the same alignment
#
def indicated_ins(coverage_list, pos):
    a = approx(coverage_list, pos, check_start=False, check_end=True)
    b = approx(coverage_list, pos, check_start=True, check_end=False)
    if a is None or b is None:
        return False
    return a == b
#
# exact hits must be present for start and end
# gap can be covered by different alignments
#
def split_ins(coverage_list, pos):
    a = exact(coverage_list, pos, check_start=False, check_end=True)
    b = exact(coverage_list, pos, check_start=True, check_end=False)
    if a is None or b is None:
        return False
    return True
#
# exact hits must be present for start or end
#
def trimmed_ins(coverage_list, pos):
    a = approx(coverage_list, pos, check_start=False, check_end=True)
    b = approx(coverage_list, pos, check_start=True, check_end=False)
    if a is None and b is None:
        return False
    return True

# p, s, t, f
def analyze(coverage_list, length, origin, sv_type, sv_size, sv_pos):
    if sv_type == "ins":
        #print(origin + sv_pos)
        if forced_ins(coverage_list, origin + sv_pos ):
            return 1, 0, 0, 0, 0, 0
        if precise_ins(coverage_list, origin + sv_pos ):
            return 0, 1, 0, 0, 0, 0
        if split_ins(coverage_list, origin + sv_pos ):
            return 0, 0, 1, 0, 0, 0
        if indicated_ins(coverage_list, origin + sv_pos ):
            return 0, 0, 0, 1, 0, 0
        if trimmed_ins(coverage_list, origin + sv_pos ):
            return 0, 0, 0, 0, 1, 0
    elif sv_type == "del":
        if forced_del(coverage_list, origin + sv_pos, origin + sv_pos + sv_size ):
            return 1, 0, 0, 0, 0, 0
        if precise_del(coverage_list, origin + sv_pos, origin + sv_pos + sv_size ):
            return 0, 1, 0, 0, 0, 0
        if split_del(coverage_list, origin + sv_pos, origin + sv_pos + sv_size ):
            return 0, 0, 1, 0, 0, 0
        if indicated_del(coverage_list, origin + sv_pos, origin + sv_pos + sv_size ):
            return 0, 0, 0, 1, 0, 0
        if trimmed_del(coverage_list, origin + sv_pos, origin + sv_pos + sv_size ):
            return 0, 0, 0, 0, 1, 0
    return 0, 0, 0, 0, 0, 1

def sv_analysis(
            file_regex, ref, out_prefix, limit_sample_files_to=10, aligners=[],
            lengths=[100, 250, 500, 1000, 5000, 10000, 50000], num_reads=1000,
            #lengths=[500, 5000, 50000], num_reads=1000,
            sv_types=["del", "ins"]
        ):
    files_list = []
    for file in glob.glob(file_regex)[:limit_sample_files_to]:
        files_list.append(file)
        #break
    # load all Chromosomes
    contigs = load_contigs(ref_fasta)
    
    read_simulator = ReadSimulator(ref, ref_fasta, contigs)
    if not read_simulator.load_sampled_from_file(out_prefix + "_sampled_distrib"):
        read_simulator.sample_distrib_from_fasta(fasta_files=files_list)
        read_simulator.save_sampled_to_file( out_prefix + "_sampled_distrib")
    
    # filter the indels, so that there are only those <= 10
    read_simulator.ins_length_distrib = [l for l in read_simulator.ins_length_distrib if l <= 10]
    read_simulator.del_length_distrib = [l for l in read_simulator.del_length_distrib if l <= 10]

    colors = ["#548235", "#eab200", "#ff6600", "#7030a0", "#aa0000", "#d0cece"]

    x_axis = []
    for x in lengths:
        for a in aligners:
            x_axis.append( (str(x) + "nt", a.get_name()) )
    y_stack = ["precise", "split", "indicated", "forced", "trimmed", "unaligned"]
    print(
            "sv_type", "sv_size", "forced", "precise", "split", "indic.", "trimmed", "failed",
            "tries", "runtime", "aligner", sep="\t"
        )
    output_file(out_prefix + "_split_reads.html")
    for sv_type in sv_types:
        plot = figure(title="SV: " + sv_type, width=1200, x_range=FactorRange(*x_axis))
        data = {"sv_sizes": x_axis}
        for y in y_stack:
            data[y] = []
        for sv_size in lengths:
            if sv_type == "del":
                read_simulator.adjust_read_lengths(sv_size, sv_size + 1000)
            if sv_type == "ins":
                read_simulator.adjust_read_lengths(-sv_size, 1000)
            read_simulator.set_sv_type(sv_size, sv_type)
            read_list = read_simulator.generate_reads(num_reads)
            query_list = []
            for sample_id, sample in enumerate(read_list):
                _, _, _, sequence, _ = sample
                query_list.append( (str(sample_id), sequence) )

            for aligner in aligners:
                runtime, alignments = aligner.align(query_list, prefix="split_reads", taskset=False)
                # group all alignments by sample id
                alignments_dict = {}
                for alignment in alignments:
                    sample_id = int(alignment.name)
                    if not sample_id in alignments_dict:
                        alignments_dict[sample_id] = []
                    alignments_dict[sample_id].append(alignment)

                num_forced = 0
                num_precise = 0
                num_indicated = 0
                num_split = 0
                num_trimmed = 0
                num_failed = 0
                num_tries = 0

                for sample_id, sample in enumerate(read_list):
                    if not sample_id in alignments_dict:
                        num_failed += 1
                        continue
                    contig, origin, length, sequence, sv_pos = sample
                    alignment_list = alignments_dict[sample_id]
                    coverage_list = alignments_to_coverage(alignment_list, contig)

                    c, p, s, i, t, f = analyze(coverage_list, length, origin, sv_type, sv_size, sv_pos)
                    if False: # indicated
                        if i == 1 and aligner.get_name()[:2] == "MA":
                            with open(".temp_ma_indicated", "a") as txt_file:
                                txt_file.write( str(origin + sv_pos) )
                                txt_file.write("\n")
                                txt_file.write(sample[-2])
                                txt_file.write("\n")
                            #print("wrote indicated read to file...")
                    if False: # forced
                        if c == 1 and aligner.get_name()[:2] == "MA":
                            with open(".temp_ma_forced", "a") as txt_file:
                                txt_file.write( str(origin + sv_pos) )
                                txt_file.write("\n")
                                txt_file.write(sample[-2])
                                txt_file.write("\n")
                            #print("wrote forced read to file...")
                    if False: # trimmed
                        if t == 1 and aligner.get_name()[:2] == "MA":
                            with open(".temp_ma_trimmed", "a") as txt_file:
                                txt_file.write( str(origin + sv_pos) )
                                txt_file.write("\n")
                                txt_file.write(sample[-2])
                                txt_file.write("\n")
                            #print("wrote trimmed read to file...")
                    if False: # failed
                        if f == 1 and aligner.get_name()[:2] == "MA":
                            with open(".temp_ma_failed", "a") as txt_file:
                                txt_file.write( str(origin + sv_pos) )
                                txt_file.write("\n")
                                txt_file.write(sample[-2])
                                txt_file.write("\n")
                            #print("wrote failed read to file...")
                    assert c + p + s + t + f + i == 1
                    num_forced += c
                    num_precise += p
                    num_split += s
                    num_indicated += i
                    num_trimmed += t
                    num_failed += f
                    num_tries += len(alignment_list)

                data["precise"].append(num_precise)
                data["forced"].append(num_forced)
                data["split"].append(num_split)
                data["indicated"].append(num_indicated)
                data["trimmed"].append(num_trimmed)
                data["unaligned"].append(num_failed)
                print(
                        sv_type, sv_size, num_forced, num_precise, num_split,
                        num_indicated, num_trimmed, num_failed, num_tries,
                        runtime, aligner.get_name(), sep="\t"
                    )
            read_simulator.adjust_read_lengths(0, 0)
            print()
        plot.vbar_stack(y_stack, x='sv_sizes', width=0.9, color=colors, source=data,
             legend=[value(x) for x in y_stack])
        show(plot)
        
        with open(out_prefix + "_sv_" + sv_type + ".json", "w") as f:
            json.dump([y_stack, data], f)
        print()

sv_analysis(
    UON_READS,
    PACK_PREFIX,
    "oxfNano",
    aligners=[
        MA(PACK_PREFIX, fast="pacBio", threads=32),
        Minimap2(PACK_PREFIX, presetting="map-ont", threads=32),
        Blasr(PACK_PREFIX, REFERENCE_FASTA, threads=32),
        G_MAP(REFERENCE_FASTA, threads=32),
        Ngmlr(REFERENCE_FASTA, threads=32),
    ]
)

sv_analysis(
    PAC_BIO_READS,
    PACK_PREFIX,
    "pacBio",
    aligners=[
        MA(PACK_PREFIX, fast="pacBio", threads=32),
        Minimap2(PACK_PREFIX, presetting="map-pb", threads=32),
        Blasr(reference, REFERENCE_FASTA, threads=32),
        G_MAP(REFERENCE_FASTA, threads=32),
        Ngmlr(REFERENCE_FASTA, threads=32),
    ]
)