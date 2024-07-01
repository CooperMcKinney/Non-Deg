import streamlit as st
import re
from collections import Counter
from io import BytesIO
import matplotlib.pyplot as plt
import numpy as np

def reverse_complement_sequence(sequence):
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return sequence.translate(complement)[::-1]

def extract_pattern(sequence, pattern):
    return re.findall(pattern, sequence)

def sort_and_count_sequences(sequences):
    sorted_sequences = sorted(sequences)
    sequence_counts = Counter(sorted_sequences)
    sorted_counts = sorted(sequence_counts.items(), key=lambda item: item[1], reverse=True)
    return sorted_counts

def normalize_counts(counts):
    max_count = max(counts, key=lambda item: item[1])[1]
    normalized = [(seq, count / max_count) for seq, count in counts]
    return normalized

def plot_normalized_counts(normalized_counts, label=None):
    ranks = list(range(1, len(normalized_counts) + 1))
    normalized_values = [count for seq, count in normalized_counts]
    plt.plot(ranks, normalized_values, marker='o', linestyle='None', label=label)

def sequence_analysis_page(page_title, key_prefix, session_state_key):
    st.title(page_title)

    uploaded_file = st.file_uploader("Upload a FASTQ or FASTQSanger file", type=["fastq", "fastqsanger", "txt"], key=key_prefix)

    if uploaded_file is not None:
        sequence = uploaded_file.getvalue().decode("utf-8").strip()
        search_pattern = st.text_input("Enter the pattern to search for (use '.' for any character):", "GTGCAC...........CTAA.A.A......TACC.GA.", key=f"{key_prefix}_search")
        trimmed_pattern = st.text_input("Trimmed sequences", "CTAA.A.A......TACC", key=f"{key_prefix}_trimmed")

        if st.button('Search Patterns', key=f"{key_prefix}_search_button"):
            reverse_complement = reverse_complement_sequence(sequence)
            matches = extract_pattern(reverse_complement, search_pattern)
            st.write("### Pattern Extraction Results")
            st.write(f"Number of sequences found: {len(matches)}")

            if trimmed_pattern:
                degseq_matches = extract_pattern("\n".join(matches), trimmed_pattern)
                sorted_counts = sort_and_count_sequences(degseq_matches)

                st.write("### Sorted Unique Sequences by Frequency")
                col1, col2 = st.columns(2)

                with col1:
                    sorted_counts_text = "\n".join([f"{count} {seq}" for seq, count in sorted_counts])
                    st.text(sorted_counts_text)

                    sorted_counts_file = BytesIO()
                    sorted_counts_file.write(sorted_counts_text.encode())
                    sorted_counts_file.seek(0)

                    st.download_button(
                        label="Download Sorted Unique Sequences",
                        data=sorted_counts_file,
                        file_name="sorted_unique_sequences.txt",
                        mime="text/plain"
                    )

                normalized_counts = normalize_counts(sorted_counts)
                st.session_state[session_state_key] = normalized_counts

                with col2:
                    plot_normalized_counts(normalized_counts, label=page_title)
                    plt.xlabel('Rank of Sequences')
                    plt.ylabel('Normalized Recombination')
                    plt.title('Rank vs. Normalized Recombination')
                    plt.grid(True)
                    plt.legend()
                    st.pyplot(plt)
                    
    elif session_state_key in st.session_state and st.session_state[session_state_key] is not None:
        normalized_counts = st.session_state[session_state_key]
        col1, col2 = st.columns(2)

        with col1:
            sorted_counts_text = "\n".join([f"{count:.4f} {seq}" for seq, count in normalized_counts])
            st.text(sorted_counts_text)

            sorted_counts_file = BytesIO()
            sorted_counts_file.write(sorted_counts_text.encode())
            sorted_counts_file.seek(0)

            st.download_button(
                label="Download Sorted Unique Sequences",
                data=sorted_counts_file,
                file_name="sorted_unique_sequences.txt",
                mime="text/plain"
            )

        with col2:
            plot_normalized_counts(normalized_counts, label=page_title)
            plt.xlabel('Rank of Sequences')
            plt.ylabel('Normalized Recombination')
            plt.title('Rank vs. Normalized Recombination')
            plt.grid(True)
            plt.legend()
            st.pyplot(plt)

    return st.session_state[session_state_key]

def average_normalized_page(file1_counts, file2_counts, file3_counts):
    st.title("Average Normalized Read Counts")

    if file1_counts and file2_counts and file3_counts:
        all_sequences = set(seq for seq, _ in file1_counts) | set(seq for seq, _ in file2_counts) | set(seq for seq, _ in file3_counts)
        average_counts = []

        for seq in all_sequences:
            counts = [
                next((count for s, count in file1_counts if s == seq), 0),
                next((count for s, count in file2_counts if s == seq), 0),
                next((count for s, count in file3_counts if s == seq), 0)
            ]
            average_count = np.mean(counts)
            average_counts.append((seq, average_count))

        average_counts.sort(key=lambda item: item[1], reverse=True)

        col1, col2 = st.columns(2)

        with col1:
            st.write("### Sorted Unique Sequences by Average Normalized Frequency")
            sorted_counts_text = "\n".join([f"{count:.4f} {seq}" for seq, count in average_counts])
            st.text(sorted_counts_text)

            sorted_counts_file = BytesIO()
            sorted_counts_file.write(sorted_counts_text.encode())
            sorted_counts_file.seek(0)

            st.download_button(
                label="Download Sorted Unique Sequences by Average Normalized Frequency",
                data=sorted_counts_file,
                file_name="sorted_unique_sequences_by_avg_norm_freq.txt",
                mime="text/plain"
            )

        with col2:
            ranks = list(range(1, len(average_counts) + 1))
            average_values = [count for seq, count in average_counts]

            plt.figure(figsize=(10, 6))
            plt.plot(ranks, average_values, marker='o', linestyle='None', label='Average')

            for file_counts, label in zip([file1_counts, file2_counts, file3_counts], ["File 1", "File 2", "File 3"]):
                counts_dict = {seq: count for seq, count in file_counts}
                file_values = [counts_dict.get(seq, 0) for seq, _ in average_counts]
                plt.plot(ranks, file_values, marker='o', linestyle='None', label=label)

            plt.xlabel('Rank of Sequences')
            plt.ylabel('Normalized Recombination')
            plt.title('Rank vs. Normalized Recombination')
            plt.grid(True)
            plt.legend()
            st.pyplot(plt)

# Initialize session state
if 'file1_counts' not in st.session_state:
    st.session_state.file1_counts = None
if 'file2_counts' not in st.session_state:
    st.session_state.file2_counts = None
if 'file3_counts' not in st.session_state:
    st.session_state.file3_counts = None

# Page navigation
page = st.selectbox("Select a page", ["File 1 Analysis", "File 2 Analysis", "File 3 Analysis", "Average Normalized Read Counts"])

if page == "File 1 Analysis":
    st.session_state.file1_counts = sequence_analysis_page("File 1 Analysis", "file1", "file1_counts")
elif page == "File 2 Analysis":
    st.session_state.file2_counts = sequence_analysis_page("File 2 Analysis", "file2", "file2_counts")
elif page == "File 3 Analysis":
    st.session_state.file3_counts = sequence_analysis_page("File 3 Analysis", "file3", "file3_counts")
elif page == "Average Normalized Read Counts":
    average_normalized_page(st.session_state.file1_counts, st.session_state.file2_counts, st.session_state.file3_counts)

