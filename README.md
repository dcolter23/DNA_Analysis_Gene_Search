# DNA Analysis and Gene Search Application
#### Video Demo:  https://youtu.be/EZF6yU4zPWs
#### Description:
For my CS50 final project, I designed a simple application in Python that allows you to generate basic DNA analysis and search the NCBI GenBank database.
As a current biotechnology master's student, I'm interested in a variety of biological concepts. Bioinformatics is one of these interests. CS50 has furthered my computational abilities and has allowed me to dive into bioinformatics. When deciding what to do for my final project, I knew I wanted to use my new abilities to gain exposure into the bioinformatics field. After doing some research, I decided to learn more about the Biopython and Streamlit libraries to build my application. The purpose of this application is to provide a simple to use tool to do some preliminary DNA and gene analysis.

This application can be run locally by running:

```console
streamlit run DNA_Analysis.py
```

The Biopython library comes with various DNA string operations that makes DNA analysis much easier. After reading the [biopython documention](https://biopython.org/wiki/Documentation) thoroughly, I discovered several modules that could be used in my application. The Streamlit library is another open-source Python library that was a central part of my project. I used Streamlit to handle the backend and the frontend of my web application. This saved time in these aspects of my project, allowing me to focus on the features of my application. I researched the different [Streamlit components](https://docs.streamlit.io/) and utilized them to quickly create the style and functionality of web application. Streamlit also allows for markdown functionality to insert custom CSS and HTML. However, this project mainly focused on the logic and programming of the application rather than the markup. The Streamlit library also provides a simple way to add additional pages to your application. This application only contains two pages, and one page is the NCBI GenBank search tool which can be easily viewed through the sidebar.

Other packages such as random, pyperclip, and datetime were used to create the various features of the application. Most of my design choices revolved around deciding how to format and save user input. For example, for the DNA text area I had to ensure that the user entered a valid DNA sequence. A valid DNA sequence must contain only letters A, T, C, and G. I decided to set a minimum length of 10 letters to avoid strings that are too short. I also decided to replace new line characters and spaces in the text area with a blank. This decision allowed the DNA sequences to be saved properly in a text file. To generate a random DNA sequence, I created a simple function with the random library to generate a valid DNA string with the user's desired DNA string length. For the GenBank search, I decided to limit the number of GenBank IDs to be shown to 100. This is to prevent the page from becoming too clustered with IDs.

One limitation of this tool is that the translation output in the DNA analysis may not be accurate for all DNA sequences inserted. This is due to different species having different codon usages for the translation of mRNA transcripts. For simplicity, I decided to use the default genetic code table, [The Standard Code (transl_table=1)](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). Therefore, it is encouraged to know how your sequence should be translated to ensure more accurate amino acid sequences.

This was CS50!

