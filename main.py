import plotly.express as px
import pandas as pd
import subprocess
import streamlit as st
import os
import stat

def my_cache(f):
    @st.cache_data(max_entries=5, ttl=600)
    @functools.wraps(f)
    def inner(*args, **kwargs):
        return f(*args, **kwargs)
    return inner

@st.cache_data
def load_data(r1_fastq, r2_fastq, primes):
    r1_fastq = r1_fastq.getvalue()
    r2_fastq = r2_fastq.getvalue()
    primes = primes.getvalue()
    
    with open("./r1_fastq.fastq.gz", 'wb') as w:
        w.write(r1_fastq)

    with open("./r2_fastq.fastq.gz", 'wb') as w:
        w.write(r2_fastq)

    with open("./primes.csv", 'wb') as w:
        w.write(primes)

    print("hi")
    
    process = subprocess.run("chmod +x ./hello; ./hello server.R", capture_output=True, shell=True)
    print(process)
    result = process.stdout.decode()
    print(result)

def main():
    st.title("NGS Enrichment Website\n")
    st.markdown("---")
    st.markdown("***Website Objective***")
    st.markdown('''NGS data analysis of yeast display repertoires is to analyse the diversity of each selection step,
                and to identify high-affinity binders by quantifying the relative abundance of different protein variants within the population after each selection round.''')
    st.markdown("---")
    st.markdown("***Parameters***")
    st.markdown("R1 and R2 fastq files in .fastq.gz form")
    st.markdown("A .csv file of primers used, i.e:")

    data = [
        {"Primer": "NGS #1", "Fwd Seq": "TCGAAGGCGGAGGGTCGGCTAGC", "Rev Seq": "CTTCGACCTCTTCAGAAATAAGCTTTTGTTCGGATCC"},
        {"Primer": "NGS #2", "Fwd Seq": "ACCTGAGCGGAGGGTCGGCTAGC", "Rev Seq": "TCAGGTCCTCTTCAGAAATAAGCTTTTGTTCGGATCC"},
        {"Primer": "NGS #3", "Fwd Seq": "GCAATCGCGGAGGGTCGGCTAGC", "Rev Seq": "GATTGCCCTCTTCAGAAATAAGCTTTTGTTCGGATCC"},
        {"Primer": "NGS #4", "Fwd Seq": "GAGATTGCGGAGGGTCGGCTAGC", "Rev Seq": "TAATCTCCTCTTCAGAAATAAGCTTTTGTTCGGATCC"},
        {"Primer": "NGS #5", "Fwd Seq": "ATCACGGCGGAGGGTCGGCTAGC", "Rev Seq": "CGATGTCCTCTTCAGAAATAAGCTTTTGTTCGGATCC"},
        {"Primer": "NGS #6", "Fwd Seq": "TTAGGCGCGGAGGGTCGGCTAGC", "Rev Seq": "TGACCACCTCTTCAGAAATAAGCTTTTGTTCGGATCC"},
        {"Primer": "NGS #7", "Fwd Seq": "ACAGTGGCGGAGGGTCGGCTAGC", "Rev Seq": "GCCAATCCTCTTCAGAAATAAGCTTTTGTTCGGATCC"},
        {"Primer": "NGS #8", "Fwd Seq": "CAGATCGCGGAGGGTCGGCTAGC", "Rev Seq": "ACTTGACCTCTTCAGAAATAAGCTTTTGTTCGGATCC"},
        {"Primer": "NGS #9", "Fwd Seq": "ATGCACGCGGAGGGTCGGCTAGC", "Rev Seq": "TACATGCCTCTTCAGAAATAAGCTTTTGTTCGGATCC"}
    ]

    st.dataframe(pd.DataFrame(data))
    st.markdown("---")
    st.markdown("***Github*** [link](%s)" % "https://github.com/codergirl1106/NGS-Enrichment-Website/")
    st.markdown("---")

    
    r1_fastq = st.file_uploader("upload R1 fastq file")
    r2_fastq = st.file_uploader("upload R2 fastq file")
    primes = st.file_uploader("upload a file containing the PCA Primers")

    if r1_fastq != None and r2_fastq != None and primes != None:

        load_data(r1_fastq, r2_fastq, primes)

        dna_seq = pd.read_csv("./dna_sequences.csv", index_col=0)

        amino_seq = pd.read_csv("./amino_sequences.csv", index_col=0)

        dna_sequences_length = list(map(lambda x: len(x), dna_seq["sequence"].tolist()))

        amino_acid_sequences_length = list(map(lambda x: len(x), amino_seq["sequence"].tolist()))

        dna_occurrence_dict = {num: dna_sequences_length.count(num) for num in set(dna_sequences_length)}

        amino_acid_occurrence_dict = {num: amino_acid_sequences_length.count(num) for num in set(amino_acid_sequences_length)}

        dna_distribution = pd.DataFrame(list(dna_occurrence_dict.items()), columns=["Length of Merged DNA Sequences (bp)", "Counts"])

        amino_acid_distribution = pd.DataFrame(list(amino_acid_occurrence_dict.items()), columns=["Length of Amino Acid Sequences", "Counts"])

        dna_dist = px.histogram(dna_distribution, x="Length of Merged DNA Sequences (bp)", y="Counts", nbins=100)

        st.markdown("***DNA Sequences Data***")
        st.markdown("---")
        dna_dist
        
        st.dataframe(dna_seq)

        st.markdown("***Amino Acid Sequences Data***")
        st.markdown("---")
        amino_dist = px.histogram(amino_acid_distribution, x="Length of Amino Acid Sequences", y="Counts", nbins=100)

        amino_dist
        
        st.dataframe(amino_seq)

main()

