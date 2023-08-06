from Bio import Entrez, SeqIO
from Bio.Seq import Seq, UndefinedSequenceError
from Bio.Blast import NCBIWWW
from datetime import date
import streamlit as st
import pyperclip

st.markdown("""
<style>
    
    [data-testid=stSidebar] {
        background-color: #800000;
    }
</style>
""", unsafe_allow_html=True)

def click_button():
        st.session_state.clicked = True

current_date = date.today()

st.title("NCBI GenBank Search")
st.text("National Center for Biotechnology Information Database")

term = st.text_input("Enter a GenBank search query")
number = st.number_input("Number of GenBank Identifiers to show", step = 1, min_value = 1, max_value = 100)

if term != "":
    handle = Entrez.esearch(db="nucleotide", term = term, retmax = number, sort = "relevance")
    record = Entrez.read(handle)
    if record["IdList"] != []:
        st.markdown(":red[**Sorted by relevance**]")
        st.write(', '.join(record["IdList"]))
        st.text("Total GenBank Identifiers Found: " + str(record["Count"]))
        
        if len(record["IdList"]) == 1:
            handle = Entrez.efetch(db="nucleotide", id = record["IdList"][0], rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()

            st.divider()
            st.write("**Record ID:**")
            st.text(record.id)

            st.divider()
            st.write("**Record Description:**")
            st.text(record.description)

            if 'clicked' not in st.session_state:
                st.session_state.clicked = False
            
            st.divider()
            try:
                st.write("**Record Sequence:**")
                st.button("Copy", on_click=click_button)
                if st.session_state.clicked:
                    pyperclip.copy(str(record.seq))
                st.markdown(record.seq)
            except UndefinedSequenceError:
                st.markdown(":red[**NOT DEFINED**]")

            st.divider()
            st.write("**Full Record:**")
            st.text(record)
    else:
        st.text("No match.")