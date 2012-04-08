snape-pooled-mac: pooled.ml
	ocamlopt.opt  -o snape-pooled str.cmxa pooled.ml 
snape-pooled: pooled.ml
	ocamlopt.opt -ccopt -static -o snape-pooled str.cmxa pooled.ml
clean:
	rm *.cmx *.cmi pooled.o snape-pooled
