from command_line_aligner import *

def near(start_align, start_orig, end_align, end_orig, score_by_coverage=False):
    if end_align >= start_orig and start_align <= end_orig:
        if score_by_coverage:
            return (min(end_align, end_orig) - max(start_align, start_orig)) / ( end_orig - start_orig )
        else:
            return 1
    return 0

def check_accuracy(read_list, reference=None, ref_seq=None, l=None, allow_secondary=False, allow_supplementary=False, silent=True, prefix="", single_core=False):
    if l is None:
        l = [
            MA(reference, fast="pacBio"),
            BWA_MEM(reference),
            Minimap2(reference),
            Ngmlr(ref_seq),
        ]

    result_list = []
    query_list = []

    for sample_id, sample in enumerate(read_list):
        _, _, _, sequence = sample
        query_list.append( (str(sample_id), sequence) )

    for aligner in l:
        runtime, alignments = aligner.align(query_list, prefix=prefix, taskset=single_core)

        tries = 0
        found_list = []
        tried_list = []
        cov_list = []
        for _ in range(len(read_list)):
            found_list.append(False)
            tried_list.append(False)
            cov_list.append(0.0)

        for alignment in alignments:
            if alignment.secondary and not allow_secondary:
                continue
            if alignment.supplementary and not allow_supplementary:
                continue
            sample_id = int(alignment.name)

            contig, origin, length, sequence = read_list[sample_id]

            tried_list[sample_id] = True

            if alignment.contig_name == contig:
                if near(
                        alignment.start, origin, 
                        alignment.start + alignment.length, origin + length
                    ):
                    found_list[sample_id] = True
                    cov = near(
                            alignment.start, origin, 
                            alignment.start + alignment.length, origin + length, 
                            True
                        )
                    if cov > cov_list[sample_id]:
                        cov_list[sample_id] = cov
            tries += 1

        one = 0
        cov = 0.0
        none = 0
        abandon = 0
        failed_list = []
        for sample_id, tried in enumerate(tried_list):
            if not tried:
                abandon += 1
        for sample_id, found in enumerate(found_list):
            if found:
                one += 1
                cov += cov_list[sample_id]
            else:
                none += 1
                failed_list.append(query_list[sample_id][1])
                #print("WARNING: aligner", name, "found none for", sample_id, read_list[sample_id])

        result_list.append( (aligner.get_name(), one, none, len(read_list), tries, runtime, cov, abandon, failed_list) )

    if not silent:
        for name, one, none, num, tries, runtime, cov, abandon, failed_list in result_list:
            print(name, "found \t:", one, "(average coverage", cov/max(one,1), ") missed:", none,
                "of", num, " simulated pac_bio reads using", tries, "alignments and",
                runtime, "seconds", abandon, "queries were abandoned")
    return result_list
