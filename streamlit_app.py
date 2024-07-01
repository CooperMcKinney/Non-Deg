import streamlit as st
import os
import re
from collections import Counter
from io import BytesIO

def reverse_complement_sequence(sequence):
    # Define the complement mapping
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    # Reverse complement the sequence
    return sequence.translate(complement)[::-1]

def extract_pattern(sequence, pattern):
    return re.findall(pattern, sequence)

def sort_and_count_sequences(sequences):
    sorted_sequences = sorted(sequences)
    sequence_counts = Counter(sorted_sequences)
    sorted_counts = sorted(sequence_counts.items(), key=lambda item: item[1], reverse=True)
    return sorted_counts

st.title("UMI Sequence Analysis")

uploaded_file = st.file_uploader("Upload a FASTQ or FASTQSanger file", type=["fastq", "fastqsanger", "txt"])

if uploaded_file is not None:
    sequence = uploaded_file.getvalue().decode("utf-8").strip()
    
    reverse_complement = reverse_complement_sequence(sequence)
    
    pattern = "GTGCAC...........CTAA.A.A......TACC.GA."
    matches = extract_pattern(reverse_complement, pattern)
    
    degseq_pattern = "CTAA.A.A......TACC"
    degseq_matches = extract_pattern("\n".join(matches), degseq_pattern)
    
    sorted_counts = sort_and_count_sequences(degseq_matches)
    
    st.write("### Sorted Unique Sequences by Frequency")
    sorted_counts_text = "\n".join([f"{count} {seq}" for seq, count in sorted_counts])
    st.text(sorted_counts_text)
    
    # Prepare the sorted counts for download
    sorted_counts_file = BytesIO()
    sorted_counts_file.write(sorted_counts_text.encode())
    sorted_counts_file.seek(0)
    
    st.download_button(
        label="Download Sorted Unique Sequences",
        data=sorted_counts_file,
        file_name="sorted_unique_sequences.txt",
        mime="text/plain"
    )