import argparse
from Bio import SeqIO
import pandas as pd

def count_repeats(sequence, mono_nucleotide, min_length, max_length):
    """
    Counts the occurrences of mononucleotide repeats in a DNA sequence.

    Parameters:
        sequence (str): The DNA sequence to analyze.
        mono_nucleotide (str): The nucleotide to search for (e.g., 'A', 'T', 'G', 'C').
        min_length (int): The minimum repeat length to consider.
        max_length (int): The maximum repeat length to consider.

    Returns:
        dict: A dictionary mapping repeat lengths to their counts in the sequence.
    """
    counts = {}
    sequence = str(sequence)

    # Loop from max_length down to min_length
    for length in range(max_length, min_length - 1, -1):
        # Create a repeat string (e.g., 'AAAA' for length=4 and mono_nucleotide='A')
        repeat = mono_nucleotide * length
        # Count the number of times the repeat appears in the sequence
        count = sequence.count(repeat)
        counts[length] = count

        # Remove the counted repeat sequences from the original sequence
        sequence = sequence.replace(repeat, '')

    return counts

def analyze_fasta(file_name, min_length, max_length):
    """
    Analyzes a FASTA file to count mononucleotide repeats and calculate their percentages.

    Parameters:
        file_name (str): The path to the input FASTA file.
        min_length (int): The minimum repeat length to consider.
        max_length (int): The maximum repeat length to consider.

    Returns:
        tuple: A tuple containing two DataFrames:
            - df_counts (pandas.DataFrame): Counts of mononucleotide repeats.
            - df_percentages (pandas.DataFrame): Percentages of mononucleotide repeats.
    """
    # The mononucleotides to check
    mono_nucleotides = ['A', 'T', 'G', 'C']

    # Create column names for the DataFrames
    columns = ['Sample']
    for mono_nucleotide in mono_nucleotides:
        for length in range(min_length, max_length + 1):
            columns.append(f'{mono_nucleotide}-{length}')

    # Create empty DataFrames for counts and percentages
    df_counts = pd.DataFrame(columns=columns)
    df_percentages = pd.DataFrame(columns=columns)

    # Parse the input FASTA file and analyze each sequence
    with open(file_name, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sample_name = record.id
            sequence = record.seq
            total_length = len(sequence)

            data_counts = {'Sample': sample_name}
            data_percentages = {'Sample': sample_name}

            # Count repeats for each mononucleotide
            for mono_nucleotide in mono_nucleotides:
                counts = count_repeats(sequence, mono_nucleotide, min_length, max_length)

                for length, count in counts.items():
                    column_name = f'{mono_nucleotide}-{length}'
                    data_counts[column_name] = count
                    data_percentages[column_name] = (count * length / total_length) * 100

            # Append the results to the DataFrames
            df_counts = df_counts.append(data_counts, ignore_index=True)
            df_percentages = df_percentages.append(data_percentages, ignore_index=True)

    return df_counts, df_percentages

def main():
    """
    Main function to parse command-line arguments and analyze the input FASTA file.
    Outputs two CSV files: counts_output.csv and percentages_output.csv.
    """
    parser = argparse.ArgumentParser(description='Count mononucleotide SSRs in a FASTA file.')
    parser.add_argument('fasta_file', type=str, help='The input FASTA file.')
    parser.add_argument('min_length', type=int, help='The minimum SSR length to consider.')
    parser.add_argument('max_length', type=int, help='The maximum SSR length to consider.')

    args = parser.parse_args()

    # Analyze the FASTA file and get the counts and percentages DataFrames
    df_counts, df_percentages = analyze_fasta(args.fasta_file, args.min_length, args.max_length)

    # Save the outputs to CSV files
    df_counts.to_csv('counts_output.csv', index=False)
    df_percentages.to_csv('percentages_output.csv', index=False)

if __name__ == "__main__":
    main()