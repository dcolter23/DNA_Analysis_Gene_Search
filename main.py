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

def click_button():
        st.session_state.clicked = True

st.title("DNA Statistics")
raw_sequence = st.text_area("Enter a DNA sequence with a length of at least 10 nucleotides.")
formated_seq = raw_sequence.replace(" ", "").replace("\n", "").upper()
isValid = False

if raw_sequence != "":
    if validSequence(formated_seq) == False or len(formated_seq) < 10:
        st.write("Not a valid sequence!")
        isValid = False
    else:
        isValid = True


if isValid:
    final_sequence = Seq(formated_seq)
    st.write("GC%:", 100 * gc_fraction(final_sequence))

    counts = {"A": final_sequence.count("A"), "T": final_sequence.count("T"), "C": final_sequence.count("C"), "G": final_sequence.count("G"), "Total": len(final_sequence)}

    st.text("Nucleotide Counts: \n" + "A: " 
    + str(counts["A"]) + "\nT: "
    + str(counts["T"]) + "\nC: "
    + str(counts["C"]) + "\nG: "
    + str(counts["G"]) + "\nTotal: "
    + str(counts["Total"]))

    st.write("Reverse Complement:", final_sequence.reverse_complement())
    st.write("Transcript:", final_sequence.transcribe())

    st.write("Amino Acid Sequence: ", final_sequence.translate())
    st.text("Using the standard genetic code.")
    st.text(CodonTable.unambiguous_dna_by_id[1])
    
    st.header("Save DNA Statistics")
    new_file = st.text_input("File name")

    st.button("Create file", on_click=click_button)

    if 'clicked' not in st.session_state:
        st.session_state.clicked = False
    
    if new_file != "":
        if st.session_state.clicked:
            f = open(new_file + ".txt", "w")
            f.write("DNA Sequence: \n" + str(final_sequence) + "\n\n")
            
            f.write("Nucleotide Count:\n")
            for nuc, count in counts.items():
                f.write("%s:%s\n" % (nuc, count))
            
            f.write("\nReverse Complement: \n" + str(final_sequence.reverse_complement()) + "\n")
            f.write("\nTranscript: \n" + str(final_sequence.transcribe()) + "\n")
            f.write("\nAmino Acid Sequence: \n" + str(final_sequence.transcribe()) + "\n")
            
            f.close()

            st.write("Created " + new_file + ".txt")
    if new_file == "" and st.session_state.clicked:
        st.write("No file name!")