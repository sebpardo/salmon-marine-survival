MS = ms
DRAFT = Pardo-marine-survival
MAIN_BIB = /home/mavzz/Documents/bib/library.bib
TODAY = $(shell date +%Y-%m-%d)

all: pdf draft diff

pdf: 
	latexmk -g $(MS).tex

doc:
	pandoc $(MS).tex --reference-doc=reference.docx --bibliography=subset.bib --csl ices.csl -o $(DRAFT)-$(TODAY).docx

draft: doc	
	cp $(MS).pdf $(DRAFT)-$(TODAY).pdf

rtf: 
	pdflatex $(MS).tex
	pdflatex $(MS).tex
	latex2rtf -E0 $(MS).tex

figs:
	cp -ru ../figures/ ./
#	cp -r ../figures/ ~/Dropbox/2020-ASCF-webinar/

todo:
	pandoc -o tasklist.html tasklist.md 

supp:
	latexmk -g suppmat.tex

bib: subset.bib

subset.bib: $(MS).aux
	bibtex $^
	bibtool -i $(MAIN_BIB) -x $^ -- 'expand.macros = ON' -- 'preserve.key.case = ON' > $@
	sed -i 's/\(annote\)\(\s\+= {in press}\)/year\2/' $@
	sed -i 's/\(annote\)\(\s\+= {in review}\)/year\2/' $@
	cat cosewic2010.bib >> $@

diff:
	latexdiff ms.tex.ices-submission $(MS).tex > ices-revisions-diff.tex
	pdflatex ices-revisions-diff.tex 
	bibtex ices-revisions-diff.aux 
	pdflatex ices-revisions-diff.tex 
	pdflatex ices-revisions-diff.tex 
	

