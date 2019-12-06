# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 11:07:17 2015
@author: peter
updated by ibrahim
"""

#### invocation: compare.py <reference.vcf> <output.phase>
####   - reference.vcf is the truth set (e.g. GIAB phase)
####   - output.phase  is the output of HapCut2

from collections import defaultdict
import os, sys
import random


def parse_hapblock_file(hapblock_file, use_SNP_index=True):
    blocklist = (
        []
    )  # data will be a list of blocks, each with a list tying SNP indexes to haplotypes
    try:
        with open(hapblock_file, "r") as hbf:

            for line in hbf:
                if len(line) < 3:  # empty line
                    continue
                if "BLOCK" in line:
                    blocklist.append([])
                    continue

                elements = line.strip().split("\t")
                if len(elements) < 3:  # not enough elements to have a haplotype
                    continue

                pos = (
                    int(elements[4]) - 1 if not use_SNP_index else int(elements[0]) - 1
                )
                allele1 = elements[1]
                allele2 = elements[2]

                blocklist[-1].append((pos, allele1, allele2))

    except FileNotFoundError:
        # most of the time, this should mean that the program timed out and therefore didn't produce a phase.
        print("File was missing. Error calculation will continue without it.")
        pass

    return blocklist


def parse_vcf_phase(vcf_file, use_SNP_index=False):
    block = []
    with open(vcf_file, "r") as vcf:
        i = 0
        for line in vcf:
            if line[0] == "#":
                continue

            elements = line.strip().split("\t")
            if len(elements) < 10:
                continue

            phase_data = elements[9]
            pos = int(elements[1]) - 1 if not use_SNP_index else i
            if phase_data[0:3] == "1|0" or phase_data[0:3] == "0|1":
                block.append((pos, phase_data[0:1], phase_data[2:3]))

            i += 1

    return [
        block
    ]  # we return a list containing the single block so format consistent with hapblock file format


def parse_runtime_file(runtime_file):
    with open(runtime_file, "r") as rf:
        return float(rf.readline().strip())


def count_SNPs(vcf_file):
    count = 0
    with open(vcf_file, "r") as infile:
        for line in infile:
            if line[:1] == "#":
                continue
            if len(line.strip().split("\t")) < 5:
                continue
            count += 1
    return count


# given a VCF file from a single chromosome return the name of that chromosome
# will almost always be from a single chromosome but don't assume that
def get_ref_name(vcf_file):
    with open(vcf_file, "r") as infile:
        for line in infile:
            if line[0] == "#":
                continue
            elif len(line.strip().split("\t")) < 5:
                continue
            return line.strip().split("\t")[0]
    # default
    print("ERROR")
    exit(1)


# assume hapblock list only uses genomic index
# add on a SNP index to it for the computation of AN50
not_in = 0


def create_SNP_ix(hapblock_list, vcf_file):
    global not_in
    snp_ix = 0
    approx_len = 0
    vcf_dict = dict()
    with open(vcf_file, "r") as infile:
        for line in infile:
            if line[:1] == "#":
                continue
            el = line.strip().split("\t")
            if len(el) < 5:
                continue

            genomic_pos = int(el[1]) - 1
            vcf_dict[genomic_pos] = snp_ix
            snp_ix += 1
            approx_len = genomic_pos

    new_hapblock_list = []

    not_in = 0
    for blk in hapblock_list:
        new_blk = []
        for (genomic_pos, a1, a2) in blk:
            if genomic_pos not in vcf_dict:
                not_in += 1
                continue
            new_blk.append((vcf_dict[genomic_pos], genomic_pos, a1, a2))
        new_hapblock_list.append(new_blk)

    return new_hapblock_list, approx_len


# returns a binary list representing which positions are covered or not.
def find_covered_positions(frag_file, num_snps):

    covered = [0] * num_snps

    with open(frag_file, "r") as infile:
        for line in infile:
            el = line.strip().split()
            num_blks_oldformat = int((len(el) - 3) / 2)
            num_blks_newformat = int((len(el) - 6) / 2)

            if num_blks_oldformat == int(el[0]):
                # old format
                for i in range(0, num_blks_oldformat):
                    pos = int(el[2 * i + 2]) - 1  # SNPs are 1-indexed
                    seq = el[2 * i + 3]
                    for j, seq_base in enumerate(seq):
                        snp_ix = pos + j
                        covered[snp_ix] = 1
                        # move along on qual string
            elif num_blks_newformat == int(el[0]):
                # new format

                for i in range(0, num_blks_newformat):
                    pos = int(el[2 * i + 5]) - 1  # SNPs are 1-indexed
                    seq = el[2 * i + 6]
                    for j, seq_base in enumerate(seq):
                        snp_ix = pos + j
                        covered[snp_ix] = 1
                        # move along on qual string

    return covered


# this function is needed for "counting ahead" at beginning of blocks.
# error_rate() needs to properly combine switch errors into mismatches
# such that switch errors are minimized (basically, if a block begins
# with an odd number of consecutive switch errors, it should assume
# that this is a sequence of all mismatches and not a "1-less" sequence of mismatches
# with a switch error at the end of it.
def count_consecutive_switches(t_dict, hap, allele):
    count = 0
    first_SNP = True
    switched = False

    for i, a1, a2 in hap:
        x = t_dict[i]  # base in true haplotype
        y = a1 if allele == 0 else a2  # base in assembled haplotype
        if x == "-" or y == "-":
            if first_SNP:
                continue
            else:
                break
        elif first_SNP:
            switched = t_dict[i] != y
            first_SNP = False
        elif (x != y and not switched) or (x == y and switched):
            count += 1
            switched = not switched
        else:
            break
    return count


# combine two dicts
def merge_dicts(d1, d2):
    d3 = d2.copy()
    for k, v in list(d1.items()):
        assert k not in d3
        d3[k] = v
    return d3


# the "error_result" abstraction and its overloaded addition operator are handy
# for combining results for the same chromosome across blocks (when the "ground truth"
# is a set of blocks rather than trio), and combining results across different chromosomes
# into genome wide stats
class error_result:
    def __init__(
        self,
        tool_name=None,
        dataset_name=None,
        ref=None,
        switch_count=None,
        poss_sw=None,
        mismatch_count=None,
        poss_mm=None,
        flat_count=None,
        poss_flat=None,
        phased_count=None,
        num_covered=None,
        num_snps=None,
        maxblk_snps=None,
        approx_len=None,
        runtime=None,
        AN50_spanlst=None,
        N50_spanlst=None,
        switch_loc=None,
        mismatch_loc=None,
        missing_loc=None,
    ):
        def create_dict(val, d_type):
            new_dict = defaultdict(d_type)
            if ref != None and val != None:
                new_dict[ref] = val
            return new_dict

        self.ref = set()  # set of references in this result (e.g. all chromosomes)
        if ref != None:
            self.ref.add(ref)

        self.tool_name = tool_name
        self.dataset_name = dataset_name

        # these are things that can be summed for the same reference,
        # e.g. switch counts for separate blocks are additive
        self.switch_count = create_dict(switch_count, int)
        self.poss_sw = create_dict(poss_sw, int)
        self.mismatch_count = create_dict(mismatch_count, int)
        self.poss_mm = create_dict(poss_mm, int)
        self.flat_count = create_dict(flat_count, int)
        self.poss_flat = create_dict(poss_flat, int)
        self.phased_count = create_dict(phased_count, int)
        self.AN50_spanlst = create_dict(AN50_spanlst, list)
        self.N50_spanlst = create_dict(N50_spanlst, list)

        # these are things that are non-additive properties, because they
        # refer to the whole reference and would be double-counted
        # e.g. if we combine errors for two blocks, on same chromosome, we add their errors
        # but we can't just add "num_snps", their chromosomes' total snp counts
        # so we use dictionaries to make sure these properties aren't duplicated

        self.num_covered = create_dict(num_covered, int)
        self.num_snps = create_dict(num_snps, int)
        self.maxblk_snps = create_dict(maxblk_snps, int)
        self.approx_lens = create_dict(approx_len, int)
        self.runtime = create_dict(runtime, int)

        self.switch_loc = create_dict(switch_loc, list)
        self.mismatch_loc = create_dict(mismatch_loc, list)
        self.missing_loc = create_dict(missing_loc, list)

    # combine two error rate results
    def __add__(self, other):
        new_err = error_result()

        if self.tool_name == None and other.tool_name == None:
            new_err.tool_name = None
        elif self.tool_name == None and other.tool_name != None:
            new_err.tool_name = other.tool_name
        elif self.tool_name != None and other.tool_name == None:
            new_err.tool_name = self.tool_name
        elif self.tool_name != None and other.tool_name != None:
            assert self.tool_name == other.tool_name
            new_err.tool_name = self.tool_name

        if self.dataset_name == None and other.dataset_name == None:
            new_err.dataset_name = None
        elif self.dataset_name == None and other.dataset_name != None:
            new_err.dataset_name = other.dataset_name
        elif self.dataset_name != None and other.dataset_name == None:
            new_err.dataset_name = self.dataset_name
        elif self.dataset_name != None and other.dataset_name != None:
            assert self.dataset_name == other.dataset_name
            new_err.dataset_name = self.dataset_name

        new_err.ref = self.ref.union(other.ref)
        new_err.switch_count = merge_dicts(self.switch_count, other.switch_count)
        new_err.poss_sw = merge_dicts(self.poss_sw, other.poss_sw)
        new_err.mismatch_count = merge_dicts(self.mismatch_count, other.mismatch_count)
        new_err.poss_mm = merge_dicts(self.poss_mm, other.poss_mm)
        new_err.flat_count = merge_dicts(self.flat_count, other.flat_count)
        new_err.poss_flat = merge_dicts(self.poss_flat, other.poss_flat)
        new_err.phased_count = merge_dicts(self.phased_count, other.phased_count)
        new_err.AN50_spanlst = merge_dicts(self.AN50_spanlst, other.AN50_spanlst)
        new_err.N50_spanlst = merge_dicts(self.N50_spanlst, other.N50_spanlst)
        new_err.num_covered = merge_dicts(self.num_covered, other.num_covered)
        new_err.num_snps = merge_dicts(self.num_snps, other.num_snps)
        new_err.maxblk_snps = merge_dicts(self.maxblk_snps, other.maxblk_snps)
        new_err.approx_lens = merge_dicts(self.approx_lens, other.approx_lens)
        new_err.runtime = merge_dicts(self.runtime, other.runtime)
        new_err.switch_loc = merge_dicts(self.switch_loc, other.switch_loc)
        new_err.mismatch_loc = merge_dicts(self.mismatch_loc, other.mismatch_loc)
        new_err.missing_loc = merge_dicts(self.missing_loc, other.missing_loc)
        return new_err

    def update_runtime(self, ref, runtime_file):
        if os.path.isfile(runtime_file):
            self.runtime[ref] = parse_runtime_file(runtime_file)
        else:
            self.runtime[ref] = 0

    def get_num_covered(self):
        return sum(self.num_covered.values())

    def get_switch_count(self):
        return sum(self.switch_count.values())

    def get_mismatch_count(self):
        return sum(self.mismatch_count.values())

    def get_flat_count(self):
        return sum(self.flat_count.values())

    def get_poss_sw(self):
        return sum(self.poss_sw.values())

    def get_poss_mm(self):
        return sum(self.poss_mm.values())

    def get_poss_flat(self):
        return sum(self.poss_flat.values())

    def get_num_snps(self):
        return sum(self.num_snps.values())

    def get_approx_len(self):
        return sum(self.approx_lens.values())

    def get_phased_count(self):
        return sum(self.phased_count.values())

    # error rate accessor functions
    def get_switch_rate(self):
        switch_count = self.get_switch_count()
        poss_sw = self.get_poss_sw()

        if poss_sw > 0:
            return switch_count / poss_sw
        else:
            return 0

    def get_mismatch_rate(self):
        mismatch_count = self.get_mismatch_count()
        poss_mm = self.get_poss_mm()

        if poss_mm > 0:
            return mismatch_count / poss_mm
        else:
            return 0

    def get_switch_mismatch_rate(self):
        poss_mm = self.get_poss_mm()

        if poss_mm > 0:
            return (self.get_switch_count() + self.get_mismatch_count()) / poss_mm
        else:
            return 0

    def get_flat_error_rate(self):
        flat_count = self.get_flat_count()
        poss_flat = self.get_poss_flat()
        if poss_flat > 0:
            return flat_count / poss_flat
        else:
            return 0

    def get_missing_rate(self):
        num_cov = self.get_num_covered()
        if num_cov > 0:
            return 1 - sum(self.phased_count.values()) / num_cov
        else:
            return 0

    def get_AN50(self):
        AN50 = 0
        AN50_spanlst = sum(list(self.AN50_spanlst.values()), [])
        AN50_spanlst.sort(reverse=True)
        phased_sum = 0
        for span, phased in AN50_spanlst:
            phased_sum += phased
            if phased_sum > self.get_num_snps() / 2:
                AN50 = span
                break
        return AN50

    def get_N50(self):
        N50 = 0
        N50_spanlst = sum(list(self.N50_spanlst.values()), [])
        N50_spanlst.sort(reverse=True)

        total = 0
        for span in N50_spanlst:
            total += span
            if total > self.get_approx_len() / 2:
                N50 = span
                break
        return N50

    def get_max_blk_snp_percent(self):
        snps_in_max_blks = sum(self.maxblk_snps.values())
        sum_all_snps = self.get_num_snps()

        if sum_all_snps > 0:
            return snps_in_max_blks / sum_all_snps
        else:
            return 0

    def get_runtime(self):
        return sum(self.runtime.values())

    def __str__(self):
        s = """
tool:            {}
dataset:         {}
switch rate:     {}
mismatch rate:   {}
flat rate:       {}
missing rate:    {}
switch errors:   {}
poss. switch:    {}
mismatch errors: {}
poss. mismatch:  {}
flat errors:     {}
poss. flat:      {}
phased count:    {}
num covered:     {}
AN50:            {}
N50:             {}
max blk snp %:   {}
runtime:         {}
missed vcf:      {}
        """.format(
            self.tool_name,
            self.dataset_name,
            self.get_switch_rate(),
            self.get_mismatch_rate(),
            self.get_flat_error_rate(),
            self.get_missing_rate(),
            self.get_switch_count(),
            self.get_poss_sw(),
            self.get_mismatch_count(),
            self.get_poss_mm(),
            self.get_flat_count(),
            self.get_poss_flat(),
            self.get_phased_count(),
            self.get_num_covered(),
            self.get_AN50(),
            self.get_N50(),
            self.get_max_blk_snp_percent(),
            self.get_runtime(),
            not_in,
        )
        return s


# compute error rates by using phase data in a VCF as ground truth
# requires VCF to have trio phase information
def hapblock_vcf_error_rate(
    assembly_file,
    frag_file,
    vcf_file,
    runtime_file=None,
    use_SNP_index=True,
    tool_name=None,
    dataset_name=None,
    largest_blk_only=False,
):
    # parse and get stuff to compute error rates
    t_blocklist = parse_vcf_phase(vcf_file, use_SNP_index)
    a_blocklist = parse_hapblock_file(assembly_file, use_SNP_index)
    # compute error result object, update the runtime and AN50 / completeness

    if largest_blk_only:
        largest_blk = []
        for blk in a_blocklist:
            if len(blk) > len(largest_blk):
                largest_blk = blk

        a_blocklist = [largest_blk]

    err = error_rate_calc(
        t_blocklist,
        a_blocklist,
        vcf_file,
        frag_file,
        runtime_file,
        use_SNP_index,
        tool_name=tool_name,
        dataset_name=dataset_name,
    )
    return err


# num_covered should be the number of SNPs with coverage in the fragment matrix file
# debug returns extra lists with indexes of errors
# t_ = truth
# a_ = hapcut
def error_rate_calc(
    t_blocklist,
    a_blocklist,
    vcf_file,
    frag_file=None,
    runtime_file=None,
    use_SNP_index=True,
    phase_set=None,
    tool_name=None,
    dataset_name=None,
):

    ref_name = get_ref_name(vcf_file)
    num_snps = count_SNPs(vcf_file)
    num_covered = (
        sum(find_covered_positions(frag_file, num_snps)) if frag_file != None else 0
    )

    if use_SNP_index:
        a_blocklist_double_index, approx_len = create_genomic_ix(a_blocklist, vcf_file)
    else:
        a_blocklist_double_index, approx_len = create_SNP_ix(a_blocklist, vcf_file)

    switch_count = 0
    mismatch_count = 0
    poss_sw = 0  # count of possible positions for switch errors
    poss_mm = 0  # count of possible positions for mismatches
    phased_count = 0
    maxblk_snps = 0
    switch_loc = []
    mismatch_loc = []
    missing_loc = []
    AN50_spanlst = []
    N50_spanlst = []

    for blk in a_blocklist_double_index:

        first_pos = -1
        last_pos = -1
        first_SNP = -1
        last_SNP = -1
        pos_snp = {}
        blk_phased = 0

        for snp_ix, pos, a1, a2 in blk:

            if a1 == "-" or (phase_set != None and snp_ix not in phase_set):
                missing_loc.append(pos)
            else:
                phased_count += 1
                pos_snp[pos] = snp_ix
                blk_phased += 1
                # if first_pos == -1:
                #     first_pos = pos
                #     first_SNP = snp_ix
                # last_pos = pos
                # last_SNP = snp_ix

        if blk_phased == 0:
            AN50_spanlst.append((0, blk_phased))
            N50_spanlst.append((0))
        else:
            first_pos = min(pos_snp)
            last_pos = max(pos_snp)
            first_SNP = pos_snp[first_pos]
            last_SNP = pos_snp[last_pos]
            blk_total = last_SNP - first_SNP + 1
            AN50_spanlst.append(
                (
                    (last_pos - first_pos) * (blk_phased / (blk_total + 0.0001)),
                    blk_phased,
                )
            )
            N50_spanlst.append((last_pos - first_pos))

        if blk_phased > maxblk_snps:
            maxblk_snps = blk_phased

    for t_block in t_blocklist:

        switched = False
        last_base_was_switch = False

        # convert t_block to a dict for convenience
        t_dict = defaultdict(lambda: "-")
        for i, a1, a2 in t_block:
            t_dict[i] = a1

        # iterate over SNPs in the true and assembled haplotypes in parallel
        # i is the index of the current base. x is the current base in the true haplotype. y is the current base in the assembled haplotype.
        for a_block in a_blocklist:

            blk_switches = [0, 0]
            blk_mismatches = [0, 0]
            blk_switchlist = [[], []]
            blk_mmlist = [[], []]
            for a in [
                0,
                1,
            ]:  # choose which allele to score. this only makes a difference for minimizing switch errors vs mismatches in corner cases.

                first_SNP = True
                for blk_ix, (i, a1, a2) in enumerate(a_block):
                    y = a1 if a == 0 else a2
                    x = t_dict[i]

                    if (
                        x == "-"
                        or y == "-"
                        or (phase_set != None and i not in phase_set)
                    ):
                        continue

                    if first_SNP:
                        switched = x != y
                        if (
                            count_consecutive_switches(t_dict, a_block[blk_ix:], a) % 2
                            == 1
                        ):
                            last_base_was_switch = True
                        else:
                            last_base_was_switch = False
                        first_SNP = False
                        continue

                    # if there is a mismatch against the true haplotype and we are in a normal state,
                    # or if there is a "match" that isn't a match because we are in a switched state,
                    # then we need to flip the state again and iterate the count
                    if (x != y and not switched) or (
                        x == y and switched
                    ):  # current base is mismatched, implying a switch
                        switched = not switched  # flip the "switched" status

                        if (
                            last_base_was_switch
                        ):  # if last base was a switch then this is actually a single-base mismatch
                            # count the 2 switches as a single-base mismatch instead
                            blk_mismatches[a] += 1
                            blk_mmlist[a].append(i)
                            blk_switches[a] -= 1  # undo count from last base switch
                            if len(blk_switchlist[a]) > 0:
                                blk_switchlist[a].pop()
                            if blk_switches[a] < 0:
                                blk_switches[a] = 0
                            last_base_was_switch = False

                        else:

                            blk_switches[a] += 1
                            blk_switchlist[a].append(i)
                            last_base_was_switch = True

                    else:  # current base is not mismatched
                        last_base_was_switch = False

                # special case for switch on last base of previous a_block; should count as a mismatch
                if last_base_was_switch:
                    # count the switch as a single-base mismatch instead
                    blk_mismatches[a] += 1
                    blk_mmlist[a].append(i)
                    blk_switches[a] -= 1
                    if len(blk_switchlist[a]) > 0:
                        blk_switchlist[a].pop()

                    if blk_switches[a] < 0:
                        blk_switches[a] = 0

            if blk_switches[0] < blk_switches[1]:
                switch_count += blk_switches[0]
                mismatch_count += blk_mismatches[0]
                switch_loc += blk_switchlist[0]
                mismatch_loc += blk_mmlist[0]

            else:
                switch_count += blk_switches[1]
                mismatch_count += blk_mismatches[1]
                switch_loc += blk_switchlist[1]
                mismatch_loc += blk_mmlist[1]

        assert len(switch_loc) == switch_count
        assert len(mismatch_loc) == mismatch_count

        # tally up how many possible positions there are for switch errors and mismatches
        # count how many phased SNPs there are so we can calculate a rate of pruned SNPs

        for blk in a_blocklist:
            phased_known = 0
            for i, a1, a2 in blk:
                if (
                    t_dict[i] != "-"
                    and a1 != "-"
                    and (phase_set == None or i in phase_set)
                ):
                    phased_known += 1

            # a switch error is only possible in blocks len 4 or greater
            # this is because switches on the ends are counted as mismatches.
            # the -3 term: -1 because only between SNPs counts, and -2 for the two ends.
            if phased_known >= 4:
                poss_sw += phased_known - 3
            # a mismatch can happen in any block length 2 or more, in any position.
            if phased_known >= 2:
                poss_mm += phased_known

        poss_flat = poss_mm
        flat_count = 0

        # iterate over SNPs in the true and assembled haplotypes in parallel
        # i is the index of the current base. x is the current base in the true haplotype. y is the current base in the assembled haplotype.
        for a_block in a_blocklist:

            flat_count1 = 0
            flat_count2 = 0

            for i, a1, a2 in a_block:

                if (
                    a1 == "-"
                    or a2 == "-"
                    or t_dict[i] == "-"
                    or (phase_set != None and i not in phase_set)
                ):
                    continue

                if a1 != t_dict[i]:
                    flat_count1 += 1
                if a2 != t_dict[i]:
                    flat_count2 += 1

            if flat_count1 < flat_count2:
                flat_count += flat_count1
            else:
                flat_count += flat_count2

    runtime = -1
    if runtime_file != None:
        runtime = parse_runtime_file(runtime_file)

    total_error = error_result(
        ref=ref_name,
        tool_name=tool_name,
        dataset_name=dataset_name,
        switch_count=switch_count,
        poss_sw=poss_sw,
        mismatch_count=mismatch_count,
        poss_mm=poss_mm,
        flat_count=flat_count,
        poss_flat=poss_flat,
        phased_count=phased_count,
        num_covered=num_covered,
        num_snps=num_snps,
        maxblk_snps=maxblk_snps,
        approx_len=approx_len,
        runtime=runtime,
        AN50_spanlst=AN50_spanlst,
        N50_spanlst=N50_spanlst,
        switch_loc=switch_loc,
        mismatch_loc=mismatch_loc,
        missing_loc=missing_loc,
    )

    return total_error


if __name__ == "__main__":
    dataset_name = os.path.basename(sys.argv[2])
    result = hapblock_vcf_error_rate(
        sys.argv[2], None, sys.argv[1], dataset_name=dataset_name, use_SNP_index=False
    )
print(result)
