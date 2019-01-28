#!/usr/bin/env python
"""Script to automate genome assembly from BLASR alignments."""

import argparse
import os
import sys
import pysam

def arg_file(in_file):
    """Checks if the specified file exists."""
    if not os.path.isfile(in_file):
        raise argparse.ArgumentTypeError(
            'Error: Unable to locate file "%s"' % in_file
        )
    else:
        return in_file

def arg_frac(value):
    """Checks that the value provided is numeric and between 0 and 1."""
    if not value.replace('.', '').isdigit() or float(value) > 1:
        raise argparse.ArgumentTypeError(
            'Value must be between 0 and 1. User supplied: %s' % value
        )
    return float(value)

def arg_value(value):
    """Checks that the value provided is numeric and >= 1."""
    if not value.isdigit() or int(value) < 1:
        raise argparse.ArgumentTypeError(
            'Value must be an integer >= 1. User supplied: %s' % value
        )
    return int(value)

def parse_args():
    """Handle user-supplied and default arguments for this script."""
    parser = argparse.ArgumentParser(
        description='Script to automate joining sequences in a genome '
                    'assembly based on the alignments of contigs at the ends '
                    'of those respective sequences.'
    )
    parser.add_argument(
        '-b', '--bams', type=arg_file, nargs='+', required=True,
        help='BAM files to process.'
    )
    parser.add_argument(
        '-l', '--length', type=arg_value, nargs='?', required=True,
        help='The minimum mapped length of a contig used to join two '
             'sequences.'
    )
    parser.add_argument(
        '-d', '--distance', type=arg_value, nargs='?', required=True,
        help='Maximum distance to look for bridging contigs at the ends of '
             'sequences.'
    )
    parser.add_argument(
        '-f', '--fraction', nargs='?', type=arg_frac, default=1.0,
        help='Maximum fraction of the sequence length to look for '
             'bridging contigs at the ends of sequences. If this value is '
             'less than that specified for "-d/--distance" for an individual '
             'sequence the fractional value will be used instead. (Default '
             '1.0)'
    )
    parser.add_argument(
        '-r', '--breakpoints', type=arg_file, nargs='?', default=None,
        help='Tab-delimited file of reference sequences and coordinates to '
             'split them at.'
    )
    parser.add_argument(
        '-c', '--cutoff', nargs='?', default=None,
        help='Sequence identifier after which no further sequences should be '
             'considered for joining.'
    )
    parser.add_argument(
        '-k', '--blacklist', nargs='+', default=[],
        help='Sequences to exclude from analysis.'
    )
    parser.add_argument(
        '-t', '--table', action='store_true',
        help='Produce a detailed table of joins and the strains supporting '
             'them.'
    )
    args = parser.parse_args()
    return args

def load_breakpoints(args):
    """Parses breakpoints into a dictionary."""
    break_dict = {}
    if args.breakpoints:
        with open(args.breakpoints, 'r') as in_handle:
            for line in in_handle:
                if line[0] != '#':
                    line = line.strip().split('\t')
                    seq = line[0]
                    coords = [int(coord)-1 for coord in line[1:]] # 0-based
                    break_dict[seq] = coords
    return break_dict

def parse_seq(args, in_bam, seq_id_ord, seq_join_dict, seq_info):
    """Parses sequences from a BAM file to a sequence join dictionary."""
    seq, subseq_num, seq_end, start, end = seq_info
    if subseq_num != None:
        seq_id = '%s.%s' % (seq, subseq_num)
    else:
        seq_id = seq
    if seq_end != 'end':
        seq_id_ord.append(seq_id)
    for contig in in_bam.fetch(seq, start, end):
        contig_name = contig.query_name.split('/')[0]
        mapped_len = contig.query_length
        rev = contig.is_reverse
        if mapped_len >= args.length:
            if contig_name not in seq_join_dict:
                seq_join_dict[contig_name] = {}
            if seq_id not in seq_join_dict[contig_name]:
                seq_join_dict[contig_name][seq_id] = [seq_end, mapped_len, rev]
            else:
                if seq_end != seq_join_dict[contig_name][seq_id][0]:
                    seq_join_dict[contig_name][seq_id][0] = 'all'
                if mapped_len > seq_join_dict[contig_name][seq_id][1]:
                    seq_join_dict[contig_name][seq_id][1] = mapped_len
    return seq_id_ord

def end_swap(seq_str):
    """Swaps the suffixes for sequence ends."""
    if 'beg' in seq_str:
        return seq_str.replace('beg', 'end')
    else:
        return seq_str.replace('end', 'beg')

def parse_bam_output(seq_id_ord, seq_join_dict):
    """Parses the output of the BAM alignments to identify joins."""
    parsed_output = []
    for seq in seq_id_ord:
        for contig_name in seq_join_dict:
            if seq in seq_join_dict[contig_name]:
                seq_end, seq_len, seq_rev = seq_join_dict[contig_name][seq]
                base_seq = seq.split('.')[0]
                for other_seq in seq_id_ord:
                    other_base_seq = other_seq.split('.')[0]
                    if base_seq != other_base_seq:
                        if other_seq in seq_join_dict[contig_name]:
                            other_seq_end, other_seq_len, other_seq_rev = (
                                seq_join_dict[contig_name][other_seq])
                            if (
                                    seq_end != 'all' and
                                    other_seq_end != 'all' and
                                    len(set([seq_end, seq_rev, other_seq_end,
                                             other_seq_rev])) == 3
                            ):
                                seq_id = '%s_%s' % (seq, seq_end)
                                other_seq_id = '%s_%s' % (other_seq,
                                                          other_seq_end)
                                parsed_output.append(
                                    [seq_id, other_seq_id, contig_name,
                                     seq_len, other_seq_len]
                                )
                            elif (seq_end == 'all' and
                                  seq_end != other_seq_end):
                                other_seq_id = '%s_%s' % (other_seq,
                                                          other_seq_end)
                                if seq_rev != other_seq_rev:
                                    seq_id = '%s_%s' % (seq, other_seq_end)
                                else:
                                    seq_id = (
                                        '%s_%s' %
                                        (seq, end_swap(other_seq_end))
                                    )
                                parsed_output.append(
                                    [seq_id, other_seq_id, contig_name,
                                     seq_len, other_seq_len]
                                )
                            elif (other_seq_end == 'all' and
                                  seq_end != other_seq_end):
                                seq_id = '%s_%s' % (seq, seq_end)
                                if seq_rev != other_seq_rev:
                                    other_seq_id = '%s_%s' % (other_seq,
                                                              seq_end)
                                else:
                                    other_seq_id = (
                                        '%s_%s' %
                                        (other_seq, end_swap(seq_end))
                                    )
                                parsed_output.append(
                                    [seq_id, other_seq_id, contig_name,
                                     seq_len, other_seq_len]
                                )
    return parsed_output

def parse_bam(args, bam, break_dict):
    """Parses BAMs to a verbose list for later superscaffold construction."""
    seq_list = []
    with pysam.AlignmentFile(bam, 'rb') as in_bam:
        cutoff = False
        seq_ord = []
        seq_end_dict = {}
        for seq_desc in in_bam.header['SQ']:
            seq = seq_desc['SN']
            seq_end = seq_desc['LN']-1 # 0-based
            seq_ord.append(seq)
            seq_end_dict[seq] = seq_end
        seq_id_ord = []
        seq_join_dict = {}
        for seq in seq_ord:
            if seq in break_dict and not cutoff:
                seq_end = seq_end_dict[seq]
                subseq_num = 1
                sub_start = 0
                for sub_end in break_dict[seq]:
                    subseq_id = '%s.%s' % (seq, subseq_num)
                    subseq_len = (sub_end-sub_start)+1
                    if subseq_id not in args.blacklist:
                        cutoff_dist = min(int(subseq_len*args.fraction),
                                          args.distance)
                        seq_id_ord = parse_seq(
                            args, in_bam, seq_id_ord, seq_join_dict,
                            [seq, subseq_num, 'beg', sub_start,
                             (sub_start+cutoff_dist)-1]
                        )
                        seq_id_ord = parse_seq(
                            args, in_bam, seq_id_ord, seq_join_dict,
                            [seq, subseq_num, 'end', (sub_end-cutoff_dist)+1,
                             sub_end]
                        )
                        seq_list.append(subseq_id)
                    subseq_num += 1
                    sub_start = sub_end+1
                subseq_len = (seq_end-sub_start)+1
                cutoff_dist = min(int(subseq_len*args.fraction),
                                  args.distance)
                seq_id_ord = parse_seq(
                    args, in_bam, seq_id_ord, seq_join_dict,
                    [seq, subseq_num, 'beg', sub_start,
                     (sub_start+cutoff_dist)-1]
                )
                seq_id_ord = parse_seq(
                    args, in_bam, seq_id_ord, seq_join_dict,
                    [seq, subseq_num, 'end', (seq_end-cutoff_dist)+1, seq_end]
                )
            elif not cutoff and seq not in args.blacklist:
                seq_end = seq_end_dict[seq]
                seq_len = seq_end_dict[seq]+1
                cutoff_dist = min(int(seq_len*args.fraction),
                                  args.distance)
                seq_id_ord = parse_seq(
                    args, in_bam, seq_id_ord, seq_join_dict,
                    [seq, None, 'beg', 0, cutoff_dist-1]
                )
                seq_id_ord = parse_seq(
                    args, in_bam, seq_id_ord, seq_join_dict,
                    [seq, None, 'end', (seq_end-cutoff_dist)+1, seq_end]
                )
                seq_list.append(seq)
            if seq == args.cutoff:
                cutoff = True
        parsed_output = parse_bam_output(seq_id_ord, seq_join_dict)
    return parsed_output, seq_list

def parse_bams(args, break_dict):
    """Parses all BAMs supplied to this script into a verbose list."""
    concat_parsed_output = []
    for bam in args.bams:
        parsed_output, seq_list = parse_bam(args, bam, break_dict)
        concat_parsed_output += parsed_output
    return concat_parsed_output, seq_list

def add_strain_to_dict(seq_dict, seq_1, seq_2, strain, contig):
    """Adds strain and contig information to a dictionary of sequence joins."""
    if seq_1 not in seq_dict:
        seq_dict[seq_1] = {}
    if seq_2 not in seq_dict[seq_1]:
        seq_dict[seq_1][seq_2] = {}
    seq_dict[seq_1][seq_2][strain] = contig

def seq_report(seq_str):
    """Reformats sequence strings to report relative sequence orientation."""
    seq_str = seq_str.split('_')
    if seq_str[-1] == 'beg':
        return '%s_fwd' % '_'.join(seq_str[:-1])
    else:
        return '%s_rev' % '_'.join(seq_str[:-1])

def anchor_seq_joins(concat_parsed_output):
    """Determines the best contig for joining two sequences."""
    seq_anchor_dict = {}
    for join_list in concat_parsed_output:
        seq_1, seq_2, contig, seq_1_len, seq_2_len = join_list
        strain = contig.split('_contig_')[0]
        anchor = seq_1_len + seq_2_len
        if strain not in seq_anchor_dict:
            seq_anchor_dict[strain] = {}
        if seq_1 not in seq_anchor_dict[strain]:
            seq_anchor_dict[strain][seq_1] = anchor
        else:
            if anchor > seq_anchor_dict[strain][seq_1]:
                seq_anchor_dict[strain][seq_1] = anchor
        if seq_2 not in seq_anchor_dict[strain]:
            seq_anchor_dict[strain][seq_2] = anchor
            if anchor > seq_anchor_dict[strain][seq_2]:
                seq_anchor_dict[strain][seq_2] = anchor
    return seq_anchor_dict

def build_matched_dict(seq_anchor_dict, concat_parsed_output):
    """Builds a sequence dictionary of matched ends based on anchor scores."""
    seq_set = set()
    seq_dict = {}
    score_dict = {}
    matched_dict = {}
    for join_list in concat_parsed_output:
        seq_1, seq_2, contig, seq_1_len, seq_2_len = join_list
        strain, contig_id = contig.split('_contig_')
        anchor = seq_1_len + seq_2_len
        if (anchor == seq_anchor_dict[strain][seq_1] and
                anchor == seq_anchor_dict[strain][seq_2]):
            if seq_2 not in seq_dict:
                add_strain_to_dict(seq_dict, seq_1, seq_2, strain,
                                   contig_id)
            elif seq_1 not in seq_dict[seq_2]:
                add_strain_to_dict(seq_dict, seq_1, seq_2, strain,
                                   contig_id)
            for seq in [seq_1, seq_2]:
                seq_set.add(seq)
                seq_set.add(end_swap(seq))
    for seq in seq_dict:
        for other_seq in seq_dict[seq]:
            score = len(seq_dict[seq][other_seq])
            if score not in score_dict:
                score_dict[score] = set()
            score_dict[score].add((seq, other_seq))
    for score in sorted(score_dict, reverse=True):
        for match_pair in score_dict[score]:
            seq_1, seq_2 = match_pair
            if seq_1 not in matched_dict and seq_2 not in matched_dict:
                matched_dict[seq_1] = seq_2
                matched_dict[seq_2] = seq_1
    return seq_dict, matched_dict, seq_set

def build_super_dict(matched_dict, seq_set):
    """Joins sequences together where possible based on strain joins."""
    super_seq_dict = {}
    matched_super_seq_ends = set()
    matched_seq_set = set(matched_dict)
    super_seq_ends = seq_set - matched_seq_set
    for super_seq_end in super_seq_ends:
        if super_seq_end not in matched_super_seq_ends:
            super_seq_dict[super_seq_end] = []
            seq_end = super_seq_end
            while end_swap(seq_end) in matched_dict:
                super_seq_dict[super_seq_end].append(seq_report(seq_end))
                seq_end = matched_dict[end_swap(seq_end)]
            super_seq_dict[super_seq_end].append(seq_report(seq_end))
            matched_super_seq_ends.add(end_swap(seq_end))
    return super_seq_dict

def build_superscaffolds(concat_reformatted_output):
    """Builds superscaffolds from the concatenated parsed BAM output."""
    seq_anchor_dict = anchor_seq_joins(concat_reformatted_output)
    seq_dict, matched_dict, seq_set = build_matched_dict(
        seq_anchor_dict, concat_reformatted_output
    )
    super_seq_dict = build_super_dict(matched_dict, seq_set)
    return seq_dict, super_seq_dict

def report_superscaffolds(super_seq_dict, seq_list):
    """Prints a report of sequences joined into superscaffolds and orphans."""
    super_seq_counter = 0
    matched_seq_ids = set()
    for super_seq in sorted(super_seq_dict, key=lambda s:
                            len(super_seq_dict[s]), reverse=True):
        if len(super_seq_dict[super_seq]) > 1:
            for seq in super_seq_dict[super_seq]:
                seq_id = '_'.join(seq.split('_')[:-1])
                matched_seq_ids.add(seq_id)
            seq_id = '_'.join(super_seq_dict[super_seq][-1].split('_')[:-1])
            matched_seq_ids.add(seq_id)
            super_seq_counter += 1
            sys.stdout.write(
                'Superscaffold %s:\t%s\n' %
                (super_seq_counter, '\t'.join(super_seq_dict[super_seq]))
            )
    for seq in seq_list:
        if seq not in matched_seq_ids:
            sys.stdout.write('Orphan: %s\n' % seq)
            matched_seq_ids.add(seq)

def report_table(seq_dict):
    """Prints a detailed table of joins and the strains supporting them."""
    for seq in sorted(seq_dict, key=lambda s:
                      (float(s.split('_')[1]), s.split('_')[2])):
        for match in sorted(seq_dict[seq], key=lambda m:
                            (float(m.split('_')[1]), m.split('_')[2])):
            strain_list = []
            for strain in sorted(seq_dict[seq][match]):
                strain_list.append('%s (%s)' %
                                   (strain, seq_dict[seq][match][strain]))
            sys.stdout.write('%s\t%s\t%s\n' %
                             (seq, match, ', '.join(strain_list)))

def main():
    """The main portion of the script."""
    args = parse_args()
    break_dict = load_breakpoints(args)
    concat_parsed_output, seq_list = parse_bams(args, break_dict)
    seq_dict, super_seq_dict = build_superscaffolds(concat_parsed_output)
    report_superscaffolds(super_seq_dict, seq_list)
    if args.table:
        report_table(seq_dict)

if __name__ == "__main__":
    main()
