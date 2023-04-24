#!/bin/bash
pdflatex manuscript.tex && bibtex manuscript.aux && pdflatex manuscript.tex && open -a preview && open -a Terminal
pdflatex supplementary.tex && bibtex supplementary.aux && pdflatex supplementary.tex && open -a preview && open -a Terminal
