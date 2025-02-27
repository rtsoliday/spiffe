
all:: spiffe.ps spiffe.pdf html

spiffe.ps:
	rm -f spiffe.aux spiffe.dvi spiffe.log spiffe.toc
	latex spiffe
	latex spiffe
	dvips spiffe -o spiffe.ps

spiffe.pdf: spiffe.ps
	ps2pdf13 spiffe.ps spiffe.pdf

html:
	mkdir -p spiffe
	htlatex spiffe "html,3" "" "-dspiffe/"
	rm -f *.html *.png *.css *.4ct *.tmp *.dvi *.idv *.4tc *.lg *.xref *.log *.aux dictionary.tex *.out


clean::
	rm -f spiffe.aux spiffe.dvi spiffe.log spiffe.out spiffe.ps spiffe.pdf
	rm -rf spiffe
