import numpy as np
import json
import pysam

def read_fasta_file(filename):
    sequences = {}
    current_sequence = ""

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences[sequence_name] = current_sequence
                    current_sequence = ""
                sequence_name = line[1:]
            else:
                current_sequence += line

        if current_sequence:
            sequences[sequence_name] = current_sequence

    return sequences


def read_fasta_file2(filename):
    with open(filename, 'r') as file:
        sequence_name = ""
        current_sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    yield sequence_name, current_sequence
                    current_sequence = ""
                sequence_name = line[1:]
            else:
                current_sequence += line
        if current_sequence:
            yield sequence_name, current_sequence


def read_sam_file(filename):
    alignments = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('@'):  # Ignore lines starting with '@' (header information)
                fields = line.split('\t')
                alignment = {
                    'QNAME': fields[0],
                    'FLAG': int(fields[1]),
                    'RNAME': fields[2],
                    'POS': int(fields[3]),
                    'MAPQ': int(fields[4]),
                    'CIGAR': fields[5],
                    'RNEXT': fields[6],
                    'PNEXT': int(fields[7]),
                    'TLEN': int(fields[8]),
                    'SEQ': fields[9],
                    'QUAL': fields[10],
                }
                alignments.append(alignment)

    return alignments


def read_bam_file(filename):
    alignments = []

    with pysam.AlignmentFile(filename, "r") as file:
        for read in file:
            alignment = {
                'QNAME': read.query_name,
                'FLAG': read.flag,
                'RNAME': read.reference_name,
                'POS': read.reference_start + 1,  # pysam是0-based，SAM是1-based
                'MAPQ': read.mapping_quality,
                'CIGAR': read.cigarstring,
                'RNEXT': read.next_reference_name,
                'PNEXT': read.next_reference_start + 1 if read.next_reference_name else 0,
                'TLEN': read.template_length,
                'SEQ': read.query_sequence,
                'QUAL': read.qual if read.qual else "*",  # pysam可能返回None
            }
            alignments.append(alignment)

    return alignments


def read_chr_len_file(file_path):
    result = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                parts = line.split()
                if len(parts) == 2:
                    chrom, length = parts
                    result.append([chrom, int(length)])

    return result


def calcu_bam_dep(chr_len_list, filename, read_len):

    sample_depth = {}
    for chro in chr_len_list:
        chro_name = chro[0]
        chro_len = chro[1]
        sample_depth[chro_name] = np.zeros(chro_len)
    alignments = read_bam_file(filename)
    for idx, sam in enumerate(alignments):
        if sam['RNAME'] != '*':
            chro_name = sam['RNAME']
            pos = sam['POS']
            for j in range(pos - 1, pos + read_len - 1):
                sample_depth[chro_name][j] += 1
    return sample_depth


def read_config(file_path):
    paths = {}
    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip().strip("'")
                paths[key] = value
    return paths


def load_npz_file(file_path):
    data = np.load(file_path)
    loaded_data = {key: data[key] for key in data}
    data.close()
    return loaded_data


def load_from_json(input_file):
    with open(input_file, 'r') as jsonfile:
        standardized_depths_serializable = json.load(jsonfile)

    standardized_depths = {}
    for chro_name, standardized_depth in standardized_depths_serializable.items():
        standardized_depths[chro_name] = np.array(standardized_depth)

    return standardized_depths

