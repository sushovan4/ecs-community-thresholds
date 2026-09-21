# Two builds from one source.
#   make            -> main.pdf      11pt preprint, for arXiv
#   make esa        -> main_esa.pdf  ESA submission format for Ecological
#                                    Monographs: double-spaced, 12pt,
#                                    continuous line numbers
#   make tables     -> regenerate the application tables from results/
TEX = pdflatex -interaction=nonstopmode

all: main.pdf

main.pdf: main.tex refs.bib tab_program.tex tab_temporal.tex
	$(TEX) main.tex && bibtex main && $(TEX) main.tex && $(TEX) main.tex

esa: main.tex refs.bib
	$(TEX) -jobname=main_esa "\def\ESAMODE{}\input{main}" && \
	bibtex main_esa && \
	$(TEX) -jobname=main_esa "\def\ESAMODE{}\input{main}" && \
	$(TEX) -jobname=main_esa "\def\ESAMODE{}\input{main}"

tables:
	cd ../experiments && python3 32_tables.py

clean:
	rm -f *.aux *.log *.out *.bbl *.blg *.spl

.PHONY: all esa tables clean
