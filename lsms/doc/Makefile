all: pdf_doxygen pdf_manual

pdf_doxygen:
	doxygen
	cd latex; pdflatex refman
	cd latex; pdflatex refman

pdf_manual:
	cd manual && $(MAKE) pdf

clean:
	cd latex && $(MAKE) clean
	cd html && $(MAKE) clean
	cd manual && $(MAKE) clean
