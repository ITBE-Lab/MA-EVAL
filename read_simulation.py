import random
import math
from Bio import SeqIO
import pysam
import os
import numpy
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column, layout
from check_accuracy import *
from command_line_aligner import *
import json
from pathlib import Path
import sys
from settings import *

UPPERCASE_LETTERS = set(map(chr, range(65, 91)))
NUMBERS = set(map(str, range(0, 10)))

def write_read_list_to_fasta(reads, file_name):
    with open(file_name, "w") as file:
        for index, read in enumerate(reads):
            file.write(">read_" + str(index) + "_from_" + str(read[0]) + "_len_" + str(read[1]) + "\n")
            file.write(read[2])
            file.write("\n")

class ReadSimulator:
    def __init__(self, reference, ref_fasta, contigs):
        self.mut_prob = .0075
        self.ins_prob = .1
        self.ins_length_distrib = list(range(1, 16))
        self.del_prob = .03
        self.del_length_distrib = list(range(1, 16))
        self.read_length_distrib = list(range(1000, 10000, 1000))
        self.from_regions = None
        self.contigs = contigs
        self.adjusted_size = 0
        self.min_read_size = 0

        self.contig_prob_list = []
        tot_len = 0
        for name, seq in contigs:
            self.contig_prob_list.append(len(seq))
            tot_len += len(seq)
        for index in range(len(self.contig_prob_list)):
            self.contig_prob_list[index] /= tot_len

        self.reference = reference
        self.ref_fasta = ref_fasta

        self.structural_variant_applier = None
        self.sv_size = 0

    def get_random_contig(self):
        f = random.random()
        i = 0
        while f > self.contig_prob_list[i] and i < len(self.contigs):
            f -= self.contig_prob_list[i]
            i += 1
        return self.contigs[i]

    def get_cigars(self, fasta_file, pick="mm"):
        file_name = "/mnt/ssd0/temp/" + fasta_file.replace("/", "_")
        if pick == "ma":
            ma_cmd = WORK_SPACE_FOLDER + "aligner/ma -m acc -x " + self.reference 
            ma_cmd += " -i " + fasta_file + " -o " + file_name + ".sam"
            #print(ma_cmd)
            #exit()
            os.system(ma_cmd)
        if pick == "bwa":
            bwa_cmd = WORK_SPACE_FOLDER + "bwa/bwa mem -t 32 " + self.reference + "bwa " + fasta_file
            bwa_cmd += " > " + file_name + ".sam"
            # silent:
            #bwa_cmd += " 2> /dev/null"
            os.system(bwa_cmd)
        if pick == "mm":
            mm_cmd = WORK_SPACE_FOLDER + "minimap2/minimap2 -t 32 -x map-pb --MD -c -a "
            mm_cmd += self.ref_fasta + " "
            mm_cmd += fasta_file + " > " + file_name + ".sam"
            os.system(mm_cmd)
        if pick == "mm-ont":
            mm_cmd = WORK_SPACE_FOLDER + "minimap2/minimap2 -t 32 -x map-ont --MD -c -a "
            mm_cmd += self.ref_fasta + " "
            mm_cmd += fasta_file + " > " + file_name + ".sam"
            os.system(mm_cmd)

        for alignment in pysam.AlignmentFile(file_name + ".sam"):
            if not alignment.cigartuples is None:
                tag = None
                if alignment.has_tag("MD"):
                    tag = alignment.get_tag("MD")
                yield (alignment.query_name, alignment.cigartuples, tag)
        os.remove(file_name + ".sam")

    def evict_mq_0_reads(self, reads):
        query_list = []
        for sample_id, sample in enumerate(reads):
            _, _, _, sequence = sample
            query_list.append( (str(sample_id), sequence) )
        ret_list = []
        runtime, alignments = Minimap2(self.reference, threads=32).align(query_list, False)
        for alignment in alignments:
            if alignment.mapping_quality > 0:
                ret_list.append(reads[int(alignment.name)])
        print(
                "evicted ", len(query_list) - len(ret_list), 
                "out of", len(query_list), 
                "seeds due to mapping quality 0."
            )
        return ret_list

    def get_characteristics_from_cigar(self, cigar, tag):
        num_nuc = 0
        num_ins = 0
        num_del = 0
        num_ins_nuc = 0
        num_del_nuc = 0
        num_mm = 0
        bHaveM = False
        for operation, length in cigar:
                ## 0 -> M  BAM_CMATCH
                ## 1 -> I  BAM_CINS
                ## 2 -> D  BAM_CDEL
                ## 3 -> N  BAM_CREF_SKIP
                ## 4 -> S  BAM_CSOFT_CLIP
                ## 5 -> H  BAM_CHARD_CLIP
                ## 6 -> P  BAM_CPAD
                ## 7 -> =  BAM_CEQUAL
                ## 8 -> X  BAM_CDIFF
                ## 9 -> B  BAM_CBACK
                
                # match (=)
                if operation == 7: 
                    num_nuc += length
                # insertion (I)
                elif operation == 1: 
                    self.ins_length_distrib.append(length)
                    num_ins += 1
                    num_ins_nuc += length
                # deletion (D)
                elif operation == 2: 
                    self.del_length_distrib.append(length)
                    num_del += 1
                    num_del_nuc += length
                    num_nuc += length
                # mismatch (X)
                elif operation == 8: 
                    num_mm += length
                    num_nuc += length
                elif operation == 0:
                    bHaveM = True
                    num_nuc += length
                elif operation != 4 and operation != 5:
                    print("got symbol:", operation, "in sam output")
        if bHaveM:
            if not tag is None and num_mm == 0:
                for index in range(len(tag)):
                    if tag[index] in UPPERCASE_LETTERS:
                        if index == 0 or tag[index - 1] in NUMBERS or tag[index - 1] in UPPERCASE_LETTERS:
                            num_mm += 1
            else:
                print("cigar symbol M in alignment but", tag is None, num_mm == 0)
        return num_nuc, num_ins, num_del, num_mm, num_ins_nuc, num_del_nuc

    def sample_distrib_from_fasta(
                self,
                fasta_files = ["/mnt/ssd0/arne/3_C01/m54015_171229_224813.subreads.partial.fasta"],
                pick="mm"
            ):
        print("Sampling read characteristics from fasta...")
        self.read_length_distrib = []
        print("\tsetting up read length distribution...")
        for fasta_file in fasta_files:
            if fasta_file.endswith(".fasta"):
                for read in SeqIO.parse(fasta_file, "fasta"):
                    self.read_length_distrib.append(len(read))
            if fasta_file.endswith(".fastq"):
                for read in SeqIO.parse(fasta_file, "fastq"):
                    self.read_length_distrib.append(len(read))
        print("\tdone")
        print("\tcomputing cigars...")
        num_nuc = 0
        num_ins = 0
        num_del = 0
        num_ins_nuc = 0
        num_del_nuc = 0
        num_mm = 0
        self.ins_length_distrib = []
        self.del_length_distrib = []
        cigar_list = []
        for fasta_file in fasta_files:
            print("\t" + fasta_file)
            for _, cigar, tag in self.get_cigars(fasta_file, pick=pick):
                num_nuc_, num_ins_, num_del_, num_mm_, num_ins_nuc_, num_del_nuc_ = self.get_characteristics_from_cigar(cigar, tag)
                num_nuc += num_nuc_
                num_ins += num_ins_
                num_del += num_del_
                num_ins_nuc += num_ins_nuc_
                num_del_nuc += num_del_nuc_
                num_mm += num_mm_
        self.ins_prob = num_ins / num_nuc
        self.del_prob = num_del / num_nuc
        self.mut_prob = num_mm / num_nuc
        print("done")
        print("probabilities:")
        print("\tinsertion:", self.ins_prob)
        print("\tdeletion:", self.del_prob)
        print("\tmutation:", self.mut_prob)
        if self.mut_prob == 0:
            print("WARNING: aligner did not output any mutations...")
        print("indel sizes:")
        print(
                "\tinsertion [50% 98% 100%]:", 
                [numpy.percentile(self.ins_length_distrib, x) for x in [50, 98, 100]],
                "adjusted by length:",
                num_ins_nuc / num_nuc
            )
        print(
                "\tdeletion [50% 98% 100%]:", 
                [numpy.percentile(self.del_length_distrib, x) for x in [50, 98, 100]],
                "adjusted by length:",
                num_del_nuc / num_nuc
            )
        print(
                "read length [50% 98% 100%]:", 
                [numpy.percentile(self.read_length_distrib, x) for x in [50, 98, 100]]
            )
        print("\tdone")

    def save_sampled_to_file(self, file):
        with open(file + ".json", "w") as f:
            json.dump([
                self.ins_prob,
                self.del_prob,
                self.mut_prob,
                self.ins_length_distrib,
                self.del_length_distrib,
                self.read_length_distrib
                ], f)

    def load_sampled_from_file(self, file):
        file = file + ".json"
        if not Path(file).is_file():
            print(file, "does not exist")
            return False
        print("loading distributions from file")
        def _decode(o):
            if isinstance(o, str) or isinstance(o, unicode):
                try:
                    return float(o)
                except ValueError:
                    return o
            elif isinstance(o, dict):
                return {_decode(k): _decode(v) for k, v in o.items()}
            elif isinstance(o, list):
                return [_decode(v) for v in o]
            else:
                return o
        with open(file, "r") as f:
            json_file = json.loads(f.read(), object_hook=_decode)
            self.ins_prob, self.del_prob, self.mut_prob, self.ins_length_distrib, self.del_length_distrib, self.read_length_distrib = json_file
        print("probabilities:")
        print("\tinsertion:", self.ins_prob)
        print("\tdeletion:", self.del_prob)
        print("\tmutation:", self.mut_prob)
        if self.mut_prob == 0:
            print("WARNING: aligner did not output any mutations...")
        print("indel sizes:")
        print(
                "\tinsertion [50% 98% 100%]:", 
                [numpy.percentile(self.ins_length_distrib, x) for x in [50, 98, 100]]
            )
        print(
                "\tdeletion [50% 98% 100%]:", 
                [numpy.percentile(self.del_length_distrib, x) for x in [50, 98, 100]]
            )
        print(
                "read length [50% 98% 100%]:", 
                [numpy.percentile(self.read_length_distrib, x) for x in [50, 98, 100]]
            )
        return True


    def load_specific_regions_from_tsv(
                self, tsv_file, col_w_chr, col_w_start, col_w_end, 
                col_w_rev_comp, one_spot_only=False, min_identity=0, min_size=0
            ):
        self.from_regions = []
        with open(tsv_file, "r") as tsv:
            for line in tsv.readlines()[1:]:
                row = line.split("\t")
                if float(row[26]) < min_identity:
                    continue
                if int(row[col_w_end]) - int(row[col_w_start]) < min_size:
                    continue
                chr_id = row[col_w_chr]
                if chr_id in ["chrM"]:
                    continue
                is_rev = (row[col_w_rev_comp] == "-")
                if is_rev: # ignore all rev comp positions for now
                    continue
                name_trans_dict = {
                    "chr1"  : "CM000663.2",
                    "chr2"  : "CM000664.2",
                    "chr3"  : "CM000665.2",
                    "chr4"  : "CM000666.2",
                    "chr5"  : "CM000667.2",
                    "chr6"  : "CM000668.2",
                    "chr7"  : "CM000669.2",
                    "chr8"  : "CM000670.2",
                    "chr9"  : "CM000671.2",
                    "chr10" : "CM000672.2",
                    "chr11" : "CM000673.2",
                    "chr12" : "CM000674.2",
                    "chr13" : "CM000675.2",
                    "chr14" : "CM000676.2",
                    "chr15" : "CM000677.2",
                    "chr16" : "CM000678.2",
                    "chr17" : "CM000679.2",
                    "chr18" : "CM000680.2",
                    "chr19" : "CM000681.2",
                    "chr20" : "CM000682.2",
                    "chr21" : "CM000683.2",
                    "chr22" : "CM000684.2",
                    "chrX"  : "CM000685.2",
                    "chrY"  : "CM000686.2"
                }
                if chr_id in name_trans_dict:
                    chr_id = name_trans_dict[chr_id]
                if chr_id.startswith("chr"):
                    chr_id = chr_id[chr_id.find("_")+1:]
                if chr_id.endswith("_random"):
                    chr_id = chr_id[:len(chr_id) - len("_random")]
                if chr_id[-2] == "v":
                    chr_id = chr_id[:-2] + "." + chr_id[-1]

                self.from_regions.append( 
                        (chr_id, is_rev, int(row[col_w_start]), int(row[col_w_end]) ) 
                    )
                if one_spot_only:
                    print(
                            "generating for region:",
                            row[col_w_chr] + ":" + row[col_w_start] + "-" + row[col_w_end],
                            "other region is at:",
                            row[7] + ":" + row[8] + "-" + row[9]
                        )
                    # break here to force the generation of all data in one spot
                    break

    def mutate(self, char):
        if char in ["C", "c"]:
            return random.choice(['a', 't', 'g'])
        elif char in ["G", "g"]:
            return random.choice(['a', 't', 'c'])
        elif char in ["A", "a"]:
            return random.choice(['c', 't', 'g'])
        elif char in ["T", "t"]:
            return random.choice(['a', 'c', 'g'])
        else:
            assert False

    def draw_from_dist(self, dist):
        return dist[random.randint(0, len(dist)-1)]

    def repl_n(self, q):
        pos = 0
        while pos < len(q):
            if not q[pos].lower() in ["a", "c", "g", "t"]:
                q = q[:pos] + random.choice(["a", "c", "g", "t"]) + q[pos+1:]
            pos += 1
        return q

    ##
    # @brief apply modifications
    # @details
    # revcomp, deletions, mutations, insertions
    # the order is important.
    # code makes sure that the same nucleotide is never modified twice
    # also indels must be at least one nuc apart from each other...
    #
    def disfigure(self, q, prob_modifier):
        #deletions
        pos = 0
        while pos < len(q):
            if random.random() < self.del_prob * prob_modifier:
                q = q[:pos] + q[pos + self.draw_from_dist(self.del_length_distrib):]
                pos += 1
            pos += 1

        # mutations
        pos = 0
        while pos < len(q):
            if random.random() < self.mut_prob * prob_modifier:
                q = q[:pos] + self.mutate(q[pos]) + q[pos+1:]
            pos += 1

        # insertions
        pos = 0
        while pos < len(q):
            if random.random() < self.ins_prob * prob_modifier:
                ins_size = self.draw_from_dist(self.ins_length_distrib)
                ins_str = ""
                for _ in range(ins_size):
                    ins_str += random.choice(["a", "c", "g", "t"])
                q = q[:pos] + ins_str + q[pos:]
                pos += ins_size + 1
            else:
                pos += 1
        return q

    def set_sv_type(self, size=0, t="none"):
        def add_insertion(ins_len, q_len, q):
            pos = int(q_len/2)
            rand_seq = ""
            for _ in range(ins_len):
                rand_seq += random.choice(["a", "c", "g", "t"])
            q = q[:pos] + rand_seq + q[pos:]
            q_len += ins_len
            return q_len, q, pos + 1

        def add_deletion(del_len, q_len, q):
            pos = int(q_len/2)
            q = q[:pos] + q[pos + del_len:]
            q_len -= del_len
            return q_len, q, pos + 1
        #
        # actually set the variant applier
        #
        if t == "none":
            self.structural_variant_applier = None
        if t == "ins":
            self.structural_variant_applier = add_insertion
        if t == "del":
            self.structural_variant_applier = add_deletion
        self.sv_size = size

    def generate_read(self, prob_modifier, reverse=None):
        q_len = self.draw_from_dist(self.read_length_distrib) + self.adjusted_size
        if q_len < self.min_read_size:
            q_en = self.min_read_size
        q_len_ex = q_len + 100

        contig_name, contig = self.get_random_contig()
        #shorten the read if it does not fit into the contig...
        if len(contig) < q_len_ex:
            print("Warning read longer than contig")
            q_len_ex = len(contig) - 1
            q_len = q_len_ex - 100

        q_from = random.randint(0, len(contig)-q_len_ex)
        q_to = q_from + q_len_ex
        q = contig[q_from : q_to]
        q = self.repl_n(q)
        
        if not self.structural_variant_applier is None:
            q_len, q, sv_pos = self.structural_variant_applier(self.sv_size, q_len, q)

        # disfigure and cut to size
        q = self.disfigure(q, prob_modifier)[:q_len]

        #reverse complement
        is_rev = False
        if (reverse is None and random.randint(0,1) == 1) or (reverse == True):
            comp = {
                'A' : 'T',
                'T' : 'A',

                'G' : 'C',
                'C' : 'G'
            }# dict
            q_ = ""
            for nuc in reversed(q):
                q_ += comp[nuc.upper()]
            q = q_
            is_rev = True

        # code for structural variations
        if not self.structural_variant_applier is None:
            return contig_name, q_from, q_len, q, sv_pos

        # return read position, length and sequence
        return contig_name, q_from, q_len, q

    def adjust_read_lengths(self, adjust_size, min_read_size):
        self.adjusted_size = adjust_size
        self.min_read_size = min_read_size

    def generate_reads(self, num, prob_modifier=1):
        reads = []
        for _ in range(num):
            reads.append(self.generate_read(prob_modifier))
        return reads

    def get_mm_indel_buckets(self, fasta_files, sampler):
        mm_buckets  = []
        ins_buckets = []
        del_buckets = []

        for fasta_file in fasta_files:
            for _, cigar, tag in self.get_cigars(fasta_file, pick=sampler):
                num_nuc, num_ins, num_del, num_mm, num_ins_nuc, num_del_nuc = self.get_characteristics_from_cigar(cigar, tag)

                mm_ratio  = num_mm  / num_nuc
                ins_ratio = num_ins / num_nuc
                del_ratio = num_del / num_nuc

                mm_percent_change  = ( (mm_ratio  / self.mut_prob) - 1 ) * 100
                ins_percent_change = ( (ins_ratio / self.ins_prob) - 1 ) * 100
                del_percent_change = ( (del_ratio / self.del_prob) - 1 ) * 100

                mm_buckets .append(mm_percent_change)
                ins_buckets.append(ins_percent_change)
                del_buckets.append(del_percent_change)

        return mm_buckets, ins_buckets, del_buckets


    def make_bar_plot(self, name, l, num_buckets=35):
        plot = figure(title=name, width=1200, y_axis_type="log")

        l_sorted = sorted(l)

        start = l_sorted[0]
        end = l_sorted[-1]

        if end - start == 0:
            plot.quad(left=start-1, right=start+2, top=1, bottom=0)
            return plot

        x = []
        y = []
        zeros = []
        for i in range(num_buckets):
            # left border of each bucket
            x.append( start+(end-start)*i/num_buckets )
            # bucket height
            y.append(0)
        # right border of last bucket
        x.append(end)

        ## x = []
        ## start = -100
        ## end = 400
        ## num_buckets = int(500 / 25)-1
        ## for i in range(-100, 400, 25):
        ##     x.append(i)

        # dummy bucket height
        y.append(0)
        
        for val in l_sorted:
            x_percent = (val - start) / (end - start)

            ## if x_percent > 1:
            ##     x_percent = 1

            assert(x_percent <= 1)
            x_pos = int( num_buckets * x_percent )
            y[x_pos] += 1

        # all buckets but the last shall exclude elements on their right border
        # so there is one dummy bucket on the very right that we now add the second last bucket
        y[-2] += y[-1]
        # we then remove the dummy bucket
        y = y[:-1]

        # transform heights into percentages:
        sum_y = 0
        for y_ele in y:
            sum_y += y_ele
        y = [y_ele / sum_y for y_ele in y]

        min_y= min(0.1 if y_ele <= 0 else y_ele for y_ele in y)

        for i in range(num_buckets):
            # bucket bottom (will always be zero...)
            zeros.append(min_y/2)

        plot.quad(left=x[:-1], right=x[1:], top=y, bottom=zeros)

        return plot

    def accuracy_coverage_graph(
                self,
                out_file_name = "accuracy_coverage_graph",
                sample_files = ["/mnt/ssd0/arne/3_C01/m54015_171229_224813.subreads.partial.fasta"],
                num_reads = 1000,
                start = -100,
                end = 350,
                step = 10,
                evict_mapq_0 = False,
                allow_secondary = True,
                aligners = [],
                single_core=False
            ):
        colors = [ "red", "blue", "green", "purple", "yellow", "black", "magenta" ]

        output_file(out_file_name + ".html")
        plot1 = figure(title="accuracy graph",     width=1200)
        plot2 = figure(title="coverage graph",     width=1200)
        plot3 = figure(title="runtime graph",      width=1200)
        plot4 = figure(title="abandon rate graph", width=1200)
        plot5 = figure(title="num alignments",     width=1200)

        coverages = []
        accuracies = []
        runtimes = []
        #startup_times = []
        abandons = []
        num_al = []
        #failed_lists = []
        names = []
        for aligner in aligners:
            names.append(aligner.get_name())
            coverages.append( [] )
            accuracies.append( [] )
            runtimes.append( [] )
            abandons.append( [] )
            num_al.append( [] )
            #failed_lists.append( [] )
            #if num_reads > 100 and end > 0:
            #    startup_times.append( aligner[1].get_start_up_time() )
            #else:
            #    startup_times.append( 0 )
            #print(aligner[1].get_name(), "'s startup time is:", startup_times[-1], "seconds")
        x_pos_list = []
        while start <= end:
            x_pos_list.append(start)
            start += step

        """
        for name, one, none, num, tries, aligner.elapsed_time in result_list:
            print(name, "found \t:", one, "(average coverage", cov/one, ") missed:", none,
                "of", num, " simulated pac_bio reads using", tries, "alignments and",
                aligner.elapsed_time, "seconds")
        """

        for a, x_pos in enumerate(x_pos_list):
            prob_mod = (x_pos + 100) / 100
            print("\t", a, "/", len(x_pos_list), "prob modifier:", prob_mod)
            read_list = self.generate_reads(num_reads, prob_mod)

            if evict_mapq_0:
                read_list = self.evict_mq_0_reads(read_list)

            results = check_accuracy(
                    read_list, l=aligners, allow_secondary=allow_secondary,
                    allow_supplementary=False, prefix=out_file_name, single_core=single_core
                )
            for index, tup in enumerate(results):
                name, one, none, num, tries, elapsed_time, cov, abandon, failed_list = tup
                #failed_lists[index].extend(failed_list)
                if one == 0:
                    coverages[index].append(float('NaN'))
                else:
                    coverages[index].append(cov / num_reads)
                accuracies[index].append(one / num_reads)
                abandons[index].append(abandon / num_reads)
                runtimes[index].append(elapsed_time) #- startup_times[index])
                num_al[index].append( tries / num_reads )
                print("accuracy", aligners[index].get_name(), "=", one / num_reads)
                print("runtime", aligners[index].get_name(), "=", elapsed_time)
                # - startup_times[index])

        for aligner, coverage, accuracy, runtime, color, abandon, tries_l in zip(
                aligners, coverages, accuracies, runtimes, colors, abandons, num_al
            ):
            plot1.line(x_pos_list, accuracy, color=color, legend=aligner.get_name())
            plot2.line(x_pos_list, coverage, color=color, legend=aligner.get_name())
            plot3.line(x_pos_list, runtime,  color=color, legend=aligner.get_name())
            plot4.line(x_pos_list, abandon,  color=color, legend=aligner.get_name())
            plot5.line(x_pos_list, tries_l,  color=color, legend=aligner.get_name())

        
        with open(out_file_name + ".json", "w") as f:
            json.dump([names, x_pos_list, accuracies, coverages, runtimes, abandons], f)

        show(column(plot1, plot2, plot3, plot4, plot5))

    def read_distribution_graph(
                self,
                out_file_name = "read_distribution_graph",
                sample_files = ["/mnt/ssd0/arne/3_C01/m54015_171229_224813.subreads.partial.fasta"],
                sampler="mm"
            ):
        output_file(out_file_name + ".html")

        mm_buckets, ins_buckets, del_buckets = self.get_mm_indel_buckets(sample_files, sampler)

        plot1 = self.make_bar_plot("mismatch distribution", mm_buckets)
        plot2 = self.make_bar_plot("insertion distribution", ins_buckets)
        plot3 = self.make_bar_plot("deletion distribution", del_buckets)
        # plot4 = self.make_bar_plot(
        #         "sampled read length distribution", 
        #         self.read_length_distrib
        #     )
        # plot5 = self.make_bar_plot(
        #         "sampled insertion length distribution", 
        #         self.ins_length_distrib
        #     )
        # plot6 = self.make_bar_plot(
        #         "sampled deletion length distribution", 
        #         self.del_length_distrib
        #     )
        with open(out_file_name + ".json", "w") as f:
            json.dump([mm_buckets, ins_buckets, del_buckets], f)

        show(column(plot1, plot2, plot3))#, plot4, plot5, plot6))