# SNAPE-pooled #
## Introduction ##
**SNAPE-pooled** computes the probability distribution for the frequency of the
minor allele in a certain population, at a certain position in the genome.
Obviously if this probability is high enough, then you have a segregating position, a.k.a.
SNP.
## Futher information ##
If you want to use SNAPE-pooled, I guess you should first read the
paper _SNP calling by sequencing pooled samples_
(Raineri E, Ferretti L, Esteve-Codina A, Nevado B, Heath S, Pérez-Enciso M.
, BMC Bioinformatics. 2012 Sep 20;13:239. doi: 10.1186/1471-2105-13-239) which explains the formulae used in it.

There is a Postscript manual in the Downloads section of this site. We are distributing the source code, but
not an executable : hence you'll have to create a binary yourself, by using the enclosed Makefile. If you don't know how to do that, send me an email.

SNAPE is not a difficult software to run. The input format is SAMTOOLS' pileup and a typical incantation goes as follows:
```
./snape-pooled -nchr 10 -theta 0.001 -D 0.1 -fold folded \
-priortype informative   < example.pileup
```

## Source code ##

You can browse the source code from this website,
namely [here](https://code.google.com/p/snape-pooled/source/browse/#svn%2Ftrunk)
and you can download it using SVN ,e.g. by doing:

```
svn checkout http://snape-pooled.googlecode.com/svn/trunk/ snape-pooled-read-only
```

The code is distributed under the GNU GPLv3 License.





If you have doubts on how exactly to download or use the package, feel free to contact me or the other project administrators.

Emanuele