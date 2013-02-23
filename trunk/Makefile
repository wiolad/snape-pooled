snape-pooled-debug: pooled.ml
	ocamlc -g -o snape-pooled str.cma pooled.ml 
snape-pooled-mac: pooled.ml
	ocamlopt.opt  -o snape-pooled str.cmxa pooled.ml 
snape-pooled: pooled.ml
	ocamlopt.opt -ccopt -static -o snape-pooled str.cmxa pooled.ml
docs: snape-pooled-manual.tex
	latex $<
	dvips -o snape-pooled-manual.ps snape-pooled-manual.dvi
snape-pooled-fast: snape-pooled.c
	gcc -g -o $@ $< -lm
clean:
	rm *.cmx *.cmi *.cmo pooled.o snape-pooled *.log *.ps *.dvi *.aux
