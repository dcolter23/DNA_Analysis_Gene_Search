from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.Data import CodonTable
from Bio import SeqIO, File
import streamlit as st
from pathlib import Path
import os

nucleotides = "ATCG"

def validSequence(sequence):
    for nuc in sequence.upper():
        if nuc not in nucleotides:
            return False
    return True

st.title("DNA Statistics")
raw_sequence = st.text_input("Enter a DNA sequence with a length of at least 10 nucleotides.")
space_sequence = raw_sequence.replace(" ", "")
isValid = False

if raw_sequence != "":
    if validSequence(space_sequence) == False or len(space_sequence) < 10:
        st.write("Not a valid sequence!")
        isValid = False
    else:
        st.write("Valid sequence!")
        isValid = True

if isValid:
    final_sequence = Seq(space_sequence)
    st.write("GC%:", 100 * gc_fraction(final_sequence))
    st.write("Reverse Complement:", final_sequence.reverse_complement())
    st.write("Transcript:", final_sequence.transcribe())
    st.text("Nucleotide Counts: \n" + "A: " 
    + str(final_sequence.count("A")) + "\nT: "
    + str(final_sequence.count("T")) + "\nC: "
    + str(final_sequence.count("C")) + "\nG: "
    + str(final_sequence.count("G")) + "\nTotal: "
    + str(len(space_sequence)))

st.header("Create a Text File with Sequence")
new_file = st.text_input("File name")
f = open(new_file + ".txt", "w")
f.write(str(space_sequence))
f.close()
# with File.as_handle(file.name) as fp:
#     for record in SeqIO.parse(file.name, file_split[1][1:]):
#         print(record.id)
#     fp.close()

# print(messenger_rna.translate(table="Vertebrate Mitochondrial"))

# standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
# mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]