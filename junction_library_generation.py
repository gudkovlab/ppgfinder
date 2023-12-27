import sys
import pandas as pd
import numpy as np
import subprocess
import tempfile
import GTF
import time

def pull_sequences(path_to_genome, coordinates_list, path_to_samtools=''):
    sequences = {}
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        for coord in coordinates_list:
            temp_file.write(coord + '\n')
        temp_file_path = temp_file.name

    bash_command = f"{path_to_samtools}samtools faidx {path_to_genome} -r {temp_file_path}"
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, _ = process.communicate()

    output = output.decode()
    blocks = output.strip().split('>')[1:]
    for block in blocks:
        lines = block.strip().split('\n')
        header = lines[0]
        seq = ''.join(lines[1:])
        coord_key = header.split(' ')[0]
        sequences[coord_key] = seq

    return sequences

# def construct_transcript_sequence(group, sequences):
#     sequence = ''
#     junction_positions = []
#     for _, exon in group.iterrows():
#         coord = f"{exon['seqname']}:{exon['start']}-{exon['end']}"
#         exon_sequence = sequences.get(coord, '')
#         sequence += exon_sequence
#         if sequence:
#             junction_positions.append(len(sequence))

#     return sequence, junction_positions
    
def construct_transcript_sequence(group, sequences):
    sequence = ''
    junction_positions = []
    intron_lengths = []  # List to hold the lengths of introns

    prev_exon_end = None  # Variable to hold the end position of the previous exon

    for _, exon in group.iterrows():
        coord = f"{exon['seqname']}:{exon['start']}-{exon['end']}"
        exon_sequence = sequences.get(coord, '')
        sequence += exon_sequence

        # Calculate intron length if this is not the first exon
        if prev_exon_end is not None and sequence:
            intron_length = exon['start'] - prev_exon_end - 1
            intron_lengths.append(intron_length)
        
        # Update junction position after adding exon sequence
        if sequence:
            junction_positions.append(len(sequence))

        prev_exon_end = exon['end']

    return sequence, junction_positions, intron_lengths


def create_exon_junction_library(path_to_gtf, path_to_genome, path_to_samtools, overhang):
    print("Reading GTF data...")
    gtf_data = GTF.dataframe(path_to_gtf)

    print("Processing exons...")
    # read first column of path_to_genome+'.fai' to learn which chromosomes are in the genome. Filter only down to those chromosomes
    fai_chomosomes = pd.read_csv(path_to_genome+'.fai', sep='\t', header=None, usecols=[0])
    fai_chomosomes = fai_chomosomes[0].tolist()
    gtf_data = gtf_data[gtf_data['seqname'].isin(fai_chomosomes)]
    
    gtf_data = gtf_data[gtf_data['feature'] == 'exon']
    gtf_data[['start', 'end']] = gtf_data[['start', 'end']].apply(pd.to_numeric)
    gtf_data.sort_values(['transcript_id', 'start'], inplace=True)

    print("Fetching sequences...")
    coordinates_list = [f"{row['seqname']}:{row['start']}-{row['end']}" for _, row in gtf_data.iterrows()]
    sequences = pull_sequences(path_to_genome, coordinates_list, path_to_samtools)

    total_transcripts = len(gtf_data['transcript_id'].unique())
    processed_transcripts = 0
    next_percent_transcripts = 10

    junctions = {}
    for transcript_id, group in gtf_data.groupby('transcript_id'):
        gene_id = group['gene_id'].iloc[0]
        sequence, junction_positions, intron_lengths = construct_transcript_sequence(group, sequences)
        if len(junction_positions) < 2:
            continue
        else:
            junction_positions = junction_positions[:-1]

        for tjid, junction in enumerate(junction_positions):
            intron_len  = intron_lengths[tjid]
            left = max(0, junction - overhang)
            right = min(len(sequence), junction + overhang)
            subseq = sequence[left:right]
            junction_key = (gene_id, subseq)
            junctions.setdefault(junction_key, []).append((transcript_id, tjid, intron_len))

        processed_transcripts += 1
        current_percent_transcripts = int((processed_transcripts / total_transcripts) * 100)
        if current_percent_transcripts >= next_percent_transcripts:
            print(f"Processed {processed_transcripts}/{total_transcripts} transcripts ({current_percent_transcripts}%)")
            next_percent_transcripts += 10

    print("Consolidating junction sequences...")
    consolidated_junctions = []
    total_junctions = len(junctions)
    processed_junctions = 0
    next_percent_junctions = 10

    for (gene_id, sequence), transcript_data in junctions.items():
        junction_name = f"{gene_id}|{'|'.join([f'{tr_name}|{tj_id+1}|{intron_len}' for tr_name, tj_id,intron_len in transcript_data])}"
        consolidated_junctions.append({'junction_name': junction_name, 'sequence': sequence})
        processed_junctions += 1
        current_percent_junctions = int((processed_junctions / total_junctions) * 100)
        if current_percent_junctions >= next_percent_junctions:
            print(f"Consolidated {processed_junctions}/{total_junctions} junctions ({current_percent_junctions}%)")
            next_percent_junctions += 10

    return pd.DataFrame(consolidated_junctions)

def save_as_pickle(df, file_path):
    df.to_pickle(file_path)

def save_as_fasta(df, file_path):
    with open(file_path, 'w') as fasta_file:
        for _, row in df.iterrows():
            fasta_file.write(f">{row['junction_name']}\n{row['sequence']}\n")

def main():
    if len(sys.argv) != 6:
        print(sys.argv)
        print("Usage: python junction_library_generation.py <gtf_file> <genome_fasta_file> <output_pkl_file> <output_fasta_file> <overhang>")
        sys.exit(1)

    path_to_gtf = sys.argv[1]
    path_to_genome = sys.argv[2]
    output_pkl_file = sys.argv[3]
    output_fasta_file = sys.argv[4]
    overhang = int(sys.argv[5])

    junction_library = create_exon_junction_library(path_to_gtf, path_to_genome, '', overhang)

    save_as_pickle(junction_library, output_pkl_file)
    save_as_fasta(junction_library, output_fasta_file)

if __name__ == "__main__":
    main()
