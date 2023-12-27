#!/bin/bash

# Function to perform BLAST for a single species genome
blast_species_genome() {
    local genome_dir="$1"
    local junction_library="$2"
    local output_file="$3"

    # Extract the directory path of the junction library
    junction_library_dir=$(dirname "$junction_library")

    # Enter the genome directory
    cd "$genome_dir" || exit 1

    # Find the genome file (modify the pattern if needed)
    genome_file=$(find . -maxdepth 1 -name "*.fa" -print -quit)

    if [ -z "$genome_file" ]; then
        echo "No genome file found in $genome_dir. Skipping BLAST."
        return
    fi

    # Extract filename without extension for database name
    db_name="${genome_file%.fa}.blast.db"

    # Check if the .blast.db file already exists
    if [ -f "${db_name}.nsq" ] || [ -f "${db_name}.nin" ] || [ -f "${db_name}.nhr" ]; then
        echo "BLAST database $db_name already exists. Skipping makeblastdb."
    else
        # Run makeblastdb
        makeblastdb -in "$genome_file" -input_type fasta -dbtype nucl -out "$db_name"
        echo "makeblastdb completed for $genome_file"
    fi

    # Define the default output file in the junction library directory
    if [ -z "$output_file" ]; then
        output_file="${junction_library_dir}/$(basename ${genome_file%.fa})_blast"
    fi

    # Run blastn
    blastn -query "$junction_library" -db "$db_name" -outfmt 6 -num_alignments 1 -num_descriptions 1 -out "$output_file"
    echo "blastn completed for $genome_file. Output saved as $output_file"
}

# Main script starts here
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <genome_directory> <junction_library_path> [output_filename]"
    exit 1
fi

genome_dir="$1"
junction_library="$2"
output_filename="$3"

blast_species_genome "$genome_dir" "$junction_library" "$output_filename"
