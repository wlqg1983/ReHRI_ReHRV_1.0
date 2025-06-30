from Bio import SeqIO
import os, sys, shutil

def split_fasta(fasta_file, output_dir='.', project_id=None):
    """
    Splits a FASTA file into multiple files, each containing one sequence.

    Parameters:
    fasta_file (str): Path to the input FASTA file.
    output_dir (str): Directory where the output files will be saved. Default is the current directory.
    project_id (str): A project identifier used to name the copied FASTA file.

    Output:
    The sequences will be saved in separate FASTA files with the format chr{i}.fasta.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Copy the FASTA file to the current directory with a new name if project_id is provided
    if project_id:
        multi_sequence_file = f'multi_sequence_{project_id}.fasta'
        shutil.copy(fasta_file, multi_sequence_file)
    else:
        multi_sequence_file = fasta_file

    # Iterate over each sequence in the FASTA file
    for i, record in enumerate(SeqIO.parse(multi_sequence_file, "fasta"), start=1):
        # Modify the record ID by prepending "chr{i}"
        record.id = f'chr{i}'
        # Ensure the description only has the original header after "chr{i}"
        record.description = record.description.split(' ', 1)[-1]
        
        # Create the filename based on the index
        filename = os.path.join(output_dir, f'chr{i}.fasta')
        # Write the sequence to the new file
        with open(filename, 'w', encoding='utf-8') as f:
            SeqIO.write(record, f, "fasta")
            
    os.remove(f'multi_sequence_{project_id}.fasta') if os.path.isfile(f'multi_sequence_{project_id}.fasta') else None

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python split_fasta.py <fasta_file> <project_id> [output_dir]")
        sys.exit(1)

    fasta_file = sys.argv[1]
    project_id = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else '.'
    
    split_fasta(fasta_file, output_dir, project_id)
    
    
    
    