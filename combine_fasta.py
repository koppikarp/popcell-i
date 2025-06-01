import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_fasta_header(header):
    refseq_match = re.search(r'([NX][MR]_\d+\.\d+)', header)
    range_match = re.search(r'range=(\S+)', header)
    
    if refseq_match and range_match:
        return refseq_match.group(1), range_match.group(1)
    elif refseq_match:
        print(f"Warning: Found RefSeq ID but no range for {header}")
        return refseq_match.group(1), None
    elif range_match:
        print(f"Warning: Found range but no RefSeq ID for {header}")
        return None, range_match.group(1)
    else:
        print(f"Warning: Could not parse FASTA header: {header}")
        return None, None

def find_coding_region(sequence):
    start = re.search(r'[A-Z]', sequence)
    end = re.search(r'[A-Z]', sequence[::-1])
    
    if start and end:
        return start.start(), len(sequence) - end.start()
    else:
        return 0, len(sequence)

def combine_sequences(exon_file, coding_file, output_path):
    exon_records = list(SeqIO.parse(exon_file, "fasta"))
    coding_records = list(SeqIO.parse(coding_file, "fasta"))
    
    print(f"Loaded {len(exon_records)} records from intron file")
    print(f"Loaded {len(coding_records)} records from full sequence file")
    
    coding_dict = {}
    for record in coding_records:
        refseq_id, range_info = parse_fasta_header(record.description)
        if refseq_id and range_info:
            coding_dict[(refseq_id, range_info)] = str(record.seq)
    
    print(f"Created dictionary with {len(coding_dict)} entries from full sequence file")

    combined_records = []
    for exon_record in exon_records:
        refseq_id, range_info = parse_fasta_header(exon_record.description)
        
        if (refseq_id, range_info) in coding_dict:
            exon_seq = str(exon_record.seq)
            coding = coding_dict[(refseq_id, range_info)]
            
            # Find coding region in intron sequence
            coding_start, coding_end = find_coding_region(exon_seq)

            # Find the coding region directly from the full sequence
            coding_start_match = re.search(r'[A-Z]', coding)
            coding_end_match = re.search(r'[A-Z]', coding[::-1])
            
            if not coding_start_match or not coding_end_match:
                print(f"Warning: Could not find coding region in full sequence for {refseq_id}, {range_info}")
                continue
            
            # Extract UTR regions from full sequence
            upstream = exon_seq[:coding_start].lower()
            downstream = exon_seq[coding_end:].lower()

            # Find the coding region directly from the full sequence
            coding_start_full = re.search(r'[A-Z]', coding).start()
            coding_end_full = len(coding) - re.search(r'[A-Z]', coding[::-1]).start()

            # Adjust UTR regions extraction
            utr_5 = coding[len(upstream):coding_start_full].lower()
            utr_3 = coding[coding_end_full:len(coding) - len(downstream)].lower()

            
            # Combine sequences
            combined_seq = upstream + utr_5 + exon_seq[coding_start:coding_end] + utr_3 + downstream
            
            combined_records.append(SeqRecord(Seq(combined_seq), id=exon_record.id, description=exon_record.description))
        else:
            print("No match found in full sequence file")
    
    print(f"\nWriting {len(combined_records)} combined records to output file")
    SeqIO.write(combined_records, output_path, "fasta")
    return output_path

exon_file = 'SLC7A10_exons.fa'
coding_file = 'SLC7A10_cds.fa'
output_path = 'SLC7A10_combined.fa'
combine_sequences(exon_file, coding_file, output_path)