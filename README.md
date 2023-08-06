# DNA Analysis and Gene Search Application
#### Video Demo:  https://youtu.be/EZF6yU4zPWs
#### Description:
For my CS50 final project, I designed a simple application in Python that allows you to generate basic DNA analysis and search the NCBI GenBank database.
As a current biotechnology master's student, I'm interested in a variety of biological concepts. Bioinformatics is one of these interests. CS50 has furthered my computational abilities and has allowed me to dive into bioinformatics. When deciding what to do for my final project, I knew I wanted to use my new abilities to gain insights into the bioinformatics field. After doing some research, I decided to dive into biopython and Streamlit libraries to build my application. The purpose of this application is to provide a simple to use tool to do some preliminary DNA and gene analysis.

This application can be run locally by running:

```console
streamlit run DNA_Analysis.py
```

The biopython library comes with various DNA string operations that makes DNA analysis much easier. After reading the [biopython documention](https://biopython.org/wiki/Documentation) thoroughly, I discovered several modules within biopython that could be used in my application. The Streamlit library is another open-source Python library that was a central part of my project. I used Streamlit to help with the backend and the frontend of my web application. This saved time in these aspects of my project, allowing me to focus on the features of my application. I researched the different [Streamlit components](https://docs.streamlit.io/) and utilized them to quickly create the style and functionality of web application. Streamlit also allows for markdown functionality to insert custom CSS and HTML. However, this project mainly focused on the logic and programming of the application rather than the markup. The Streamlit library also provides a simple way to add additional pages to your application. This application only contains two pages and one page is the NCBI GenBank search tool.

Other packages such as random, pyperclip, and datetime were used to create the various features of the application. Most of my design choices revolved around deciding how to format and save user input. 

