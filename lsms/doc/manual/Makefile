TEX_FILES = lsms_manual.tex \
            constrained_moment.tex instalation.tex units.tex \
            lsms_implementation.tex \
            gpu_details.tex

all: pdf

pdf: $(TEX_FILES)
	pdflatex lsms_manual.tex
	bibtex lsms_manual
	pdflatex lsms_manual.tex
	pdflatex lsms_manual.tex
	pdflatex lsms_implementation.tex
	bibtex lsms_implementation
	pdflatex lsms_implementation.tex
	pdflatex lsms_implementation.tex

clean:
	rm -f lsms_manual.pdf
	rm -f lsms_implementation.pdf
	rm -f *.aux
	rm -f *.log
	rm -f *.toc
	rm -f *.bbl *.blg
