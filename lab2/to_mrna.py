import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

input_file = "data.fa"


filename, extension = os.path.splitext(input_file)
output_file = f"{filename}.rna{extension}"

records = SeqIO.parse(input_file, "fasta")

for record in records:
    rna_seq = record.seq.transcribe()
    rna_record = SeqRecord(rna_seq, id=record.id, description=record.description)
    SeqIO.write(rna_record, output_file, "fasta")
