snape-pooled:pooled.ml
	ocamlopt.opt -o snape-pooled str.cmxa pooled.ml
clean:
	rm *.cmx *.cmi snape-pooled
