import os
import sys
import time
import statistics
from settings import *



class CommandLine():
    class Alignment():
        def __init__(self):
            self.start = 0
            self.length = 0
            self.name = ""
            self.contig_name = ""
            self.is_forw_starand = True
            self.mapping_quality = 0
            self.secondary = False
            self.supplementary = False
            self.cigar = []

    def get_name(self):
        raise Exception("This function shall be overwritten by al children")

    def create_command(self, in_file, out_file):
        raise Exception("This function shall be overwritten by al children")

    def __get_sam(self, queries, taskset, prefix, allow_error=False):
        in_file = TEMP_FILE_FOLDER + prefix + self.get_name().replace(" ", "_") + ".fasta"
        out_file = TEMP_FILE_FOLDER + prefix + self.get_name().replace(" ", "_") + ".sam"

        with open(in_file, "w") as f:
            for name, sequence in queries:
                if len(sequence) == 0:
                    continue
                assert not "/" in name
                # write in pacbio bla format so that blasr does not cause trouble...
                f.write(">" + name + "/0/0_" + str(len(sequence)) + "\n" )
                while len(sequence) >= 60:
                    line = sequence[:60]
                    sequence = sequence[60:]
                    f.write(line + "\n")
                f.write(sequence + "\n")

        #assemble the shell command
        cmd_str = self.create_command(in_file, out_file)
        if taskset:
            cmd_str = "taskset 1 " + cmd_str

        start_time = time.time()
        returncode = os.system(cmd_str)
        elapsed_time = time.time() - start_time

        if returncode != 0:
            if allow_error:
                os.remove(in_file)
                return elapsed_time, ""
            print("call command:", cmd_str)
            print("subprocess returned with ERROR")
            exit()

        os.remove(in_file)

        sam_file = []
        with open(out_file, "r") as f:
            sam_file = f.readlines()

        os.remove(out_file)

        return elapsed_time, sam_file

    ##
    # @origin https://stackoverflow.com/questions/10321978/integer-to-bitfield-as-a-list
    def __bitfield(self, num):
        return list(reversed([True if digit=='1' else False for digit in bin(num)[2:]]))

    def __check_flag(self, string, flag_bit):
        bits = self.__bitfield(int(string))
        if len(bits) <= flag_bit:
            return False
        return bits[flag_bit]

    def secondary(self, string):
        return self.__check_flag(string, 8)

    def supplementary(self, string):
        return self.__check_flag(string, 11)

    ##
    # @details
    # queries is a list of (name, sequence) tuples (str, str)
    # returns (elapsed_time, alignments)
    # where alignments is a list of Alignment objects
    def align(self, queries, taskset=True, prefix=""):
        elapsed_time, lines = self.__get_sam(queries, taskset, prefix)
        #print(sam)

        #transform sam file into list data structure
        alignments = []

        if self.output_type() == "SAM":
            # remove the header of the sam format
            while len(lines) > 0 and ( len(lines[0]) == 0 or lines[0][0] is '@' ):
                lines = lines[1:]

            for line in lines:
                # ignore empty lines
                if len(line) == 0:
                    continue
                # ignore samples where the name starts with dummy...
                if line[:6] == "dummy_":
                    continue
                # split into columns
                columns = line.split("\t")

                try:
                    align_length = 0
                    ##
                    # brief helper function to read cigars...
                    def read_cigar(cigar):
                        if cigar[0] == '*':
                            return
                        number_start = 0
                        while number_start < len(cigar):
                            symbol_start = number_start
                            #increase symbol_start while cigar[symbol_start] is a number...
                            while cigar[symbol_start] in [str(x) for x in range(10)]:
                                symbol_start += 1
                            yield int(cigar[number_start:symbol_start]), cigar[symbol_start]
                            number_start = symbol_start + 1

                    assert(len(columns) >= 5)
                    qLen = 0
                    rLen = 0
                    cigar_list = []
                    for amount, char in read_cigar(columns[5]):
                        cigar_list.append( (amount, char) )
                        if char in ['M', 'I', '=', 'X', 'S']:
                            qLen += amount
                        if char in ['M', 'D', '=', 'X', 'N']:
                            rLen += amount
                        if char in ['M', 'X', '=', 'D']:
                            align_length += amount
                        # sanity check...
                        elif not char in ['M', 'I', 'D', '=', 'X', 'S', 'N', 'H', 'P']:
                            print("Error: got wierd cigar symbol", char, "in cigar", columns[5])
                            exit()
                    alignment = self.Alignment()

                    alignment.start = int(columns[3])
                    alignment.cigar = cigar_list
                    alignment.contig_name = columns[2]
                    alignment.length = rLen
                    alignment.mapping_quality = int(columns[4])/255
                    alignment.name = columns[0].split("/")[0]
                    alignment.secondary = self.secondary(columns[1])
                    alignment.supplementary = self.supplementary(columns[1])
                    alignments.append(alignment)
                except Exception as e:
                    print("Error:", e)
                    print(line)
                    exit()
        else:
            raise Exception(self.output_type() + " output is unsupported.")

        return elapsed_time, alignments

    def output_type(self):
        return "SAM"

    def get_start_up_time(self, num_tries=10):
        start_up_times = []
        for _ in range(num_tries):
            elapsed_time, _ = self.__get_sam([("empty_query", "")], False, ".temp", True)
            start_up_times.append(elapsed_time)
        return statistics.mean(start_up_times)

###
# The list of aligners that can be analyzed:
###

class Minimap2(CommandLine):
    def __init__(self, index_str, num_results=None, z_drop=None, presetting=None, silent=False, threads=None):
        self.minimap2_home = WORK_SPACE_FOLDER + "minimap2/"
        self.num_results = num_results
        self.presetting = presetting
        if presetting is None:
            self.index_str = index_str + ".mmi"
        else:
            self.index_str = index_str + presetting + ".mmi"
        if not z_drop is None:
            self.z_drop = " -z " + str(z_drop)
        else:
            self.z_drop = ""
        if silent:
            self.silent = " 2> /dev/null"
        else:
            self.silent = ""
        if threads is None:
            self.threads = ""
        else:
            self.threads = " -t " + str(threads)

    def create_command(self, in_file, out_file):
        cmd_str = self.minimap2_home + "minimap2 -c -a "
        cmd_str += self.threads
        if not self.num_results is None:
            in_file += " --secondary=yes -N " + self.num_results
        if not self.presetting is None:
            cmd_str += " -x " + self.presetting
        return cmd_str + " " + self.index_str + " " + in_file + self.z_drop + " > " + out_file + self.silent

    def get_name(self):
        return "Minimap 2"

class MA(CommandLine):
    def __init__(self, index_str, num_results=None, fast=True, finder_mode=False, soc_width=None, threads=None, silent=False):
        self.ma_home = WORK_SPACE_FOLDER + "aligner/"
        self.index_str = index_str
        self.num_results = num_results
        self.mode = "acc"
        if threads is None:
            self.threads = " -t 1 "
        else:
            self.threads = " -t " + str(threads)
        if isinstance(fast, bool):
            if fast:
                self.mode = "fast"
        elif isinstance(fast, str):
            self.mode = fast
        else:
            print(fast, "is neither string nor bool; aborting!")
            exit()
        self.finder_mode = finder_mode
        self.soc_width = soc_width
        if silent:
            self.silent = " 2> /dev/null"
        else:
            self.silent = ""

    def create_command(self, in_file, out_file):
        cmd_str = self.ma_home + "ma -m " + self.mode + self.threads
        if self.finder_mode:
            cmd_str += " -d"
        if not self.soc_width is None:
            cmd_str += " --SoCWidth " + self.soc_width
        if not self.num_results is None:
            cmd_str += " -n " + self.num_results
        return cmd_str + " -x " + self.index_str + " -i " + in_file + " -o " + out_file +self.silent
    
    def get_name(self):
        return "MA " + self.mode

    def output_type(self):
        if self.finder_mode:
            return "FINDER"
        else:
            return "SAM"

class BWA_MEM(CommandLine):
    def __init__(self, index_str, num_results=None, z_drop=None, presetting=None, silent=True, threads=None):
        self.bwa_home = WORK_SPACE_FOLDER + "bwa/"
        self.index_str = index_str + "bwa"
        self.num_results = num_results
        self.presetting = presetting
        if threads is None:
            self.threads = ""
        else:
            self.threads = " -t " + str(threads)
        if not z_drop is None:
            self.z_drop = " -d " + str(z_drop)
        else:
            self.z_drop = ""
        if silent:
            self.silent = " 2> /dev/null"
        else:
            self.silent = ""

    def create_command(self, in_file, out_file):
        cmd_str = self.bwa_home + "bwa mem " + self.threads
        if not self.num_results is None:
            cmd_str += " -a"
        if not self.presetting is None:
            cmd_str += " -x " + self.presetting
        return cmd_str + " " + self.index_str + " " + in_file + self.z_drop + " > " + out_file + self.silent

    def get_name(self):
        return "BWA-MEM"

class Ngmlr(CommandLine):
    def __init__(self, genome_file, threads=None):
        self.ngmlr_home = WORK_SPACE_FOLDER + "ngmlr/bin/ngmlr-0.2.8/"
        self.genome_file = genome_file
        if threads is None:
            self.threads = ""
        else:
            self.threads = " -t " + str(threads)

    def create_command(self, in_file, out_file):
        cmd_str = self.ngmlr_home + "ngmlr -r " + self.genome_file + " -q " + in_file + self.threads
        #print(cmd_str)
        return cmd_str + " > " + out_file

    def get_name(self):
        return "Ngmlr"

class Blasr(CommandLine):
    def __init__(self, index_str, genome_str, threads=None):
        self.blasr_home = WORK_SPACE_FOLDER + "blasr/build/bin/"
        self.index_str = index_str + "blasr"
        self.genome_str = genome_str
        if threads is None:
            self.threads = ""
        else:
            self.threads = " --nproc " + str(threads)

    def create_command(self, in_file, out_file):
        cmd_str = self.blasr_home + "blasr " + self.threads + " " + in_file
        # --nproc 32
        return cmd_str + " " + self.genome_str + " --sa " + self.index_str + " --sam --out " + out_file + " --cigarUseSeqMatch"

    def get_name(self):
        return "Blasr"

class G_MAP(CommandLine):
    def __init__(self, genome_str, threads=None):
        self.g_home = WORK_SPACE_FOLDER + "graphmap/bin/Linux-x64/"
        self.genome_str = genome_str
        if threads is None:
            self.threads = ""
        else:
            self.threads = " -t " + str(threads)

    def create_command(self, in_file, out_file):
        cmd_str = self.g_home + "graphmap align -v 0 -r " + self.genome_str + self.threads
        return cmd_str + " -d " + in_file + " -o " + out_file

    def get_name(self):
        return "GraphMap"

class Bowtie2(CommandLine):
    def __init__(self, index_str, threads=None):
        self.bowtie2_home = WORK_SPACE_FOLDER + "bowtie2/bowtie2-2.3.3.1/"
        self.index_str = index_str + "bowtie2"
        if threads is None:
            self.threads = ""
        else:
            self.threads = " -p " + str(threads)

    def create_command(self, in_file, out_file):
        cmd_str = self.bowtie2_home + "bowtie2 " + self.threads
        index_str = " -x " + self.index_str
        input_str = " -f -U " + in_file
        return cmd_str + " " + index_str + " " + input_str + " > " + out_file

    def do_checks(self):
        return False

    def get_name(self):
        return "Bowtie2"