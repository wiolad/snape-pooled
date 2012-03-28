<TeXmacs|1.0.7.3>

<style|tmweb>

<\body>
  <doc-data|<doc-title|SNAPE-pooled>|<\doc-author-data|<author-name|emanuele.raineri@gmail.com>>
    \;
  <|doc-author-data>
    \;
  </doc-author-data>>

  <section|introduction>

  SNAPE-pooled computes the probability distribution for the frequency of the
  minor allele in a certain population, at a certain position in the genome
  of the population. If you decide to use SNAPE-pooled, you should first read
  the accompanying paper which describes the formulae used in it.\ 

  <section|input format>

  The input data must be formatted according to the pileup specifications
  [see <with|color|blue|<with|font-family|tt|http://samtools.sourceforge.net/pileup.shtml>>],
  <em|i.e.> the following fields must be present:

  \;

  <\enumerate-numeric>
    <item>a chromosome field, part of the genomic coordinate

    <item>an integer, specifying the position along the chromosome

    <item>the reference nucleotide, <em|i.e.> the content of the reference
    genome for the population at that position of the given chromosome

    <item>the coverage, that is the number of bases in the pileup

    <item>the pileup, a list of all the nucleotides aligned with the position
    specified in (1) and (2). Each nucleotide comes from a different read,
    each read might (or not) come from a different individual.

    <item>the quality pileup, that is a quality symbol for each of the
    nucleotides in (4).

    \;
  </enumerate-numeric>

  <section|command line options>

  It is also necessary to specify some of the parameters used in the
  calculations, which can be done through a set of command line options.
  These are:

  <block|<tformat|<twith|table-width|45>|<twith|table-hmode|max>|<cwith|6|6|2|2|cell-row-span|34>|<cwith|6|6|2|2|cell-col-span|45>|<table|<row|<cell|<with|font-family|tt|nchr>>|<cell|Number
  of different individuals in the pool>>|<row|<cell|<with|font-family|tt|theta>>|<cell|<with|mode|math|\<theta\>>
  the nucleotide diversity>>|<row|<cell|<with|font-family|tt|D>>|<cell|Prior
  genetic difference between reference genome and
  population>>|<row|<cell|<with|font-family|tt|priortype>>|<cell|Can be
  <with|color|blue|informative> or <with|color|blue|flat>>>|<row|<cell|<with|font-family|tt|fold>>|<cell|<with|color|blue|folded>
  or <with|color|blue|unfolded>>>|<row|<cell|<with|font-family|tt|spectrum>>|<cell|If
  present, print the full pdf for the minor allele frequency.>>>>>

  \;

  if <with|font-family|tt|-spectrum> is not specified, ony summary values
  will be printed, see following section.

  \;

  \;

  \;

  \;

  \;

  \;

  <section|output format>

  the output contains a minimum of <math|10> fields,
  <with|font-family|tt|TAB>-separated, as in the following list:

  <\enumerate>
    <item>chr (1) and (2) are the genomic coordinates

    <item>position along the chromosome

    <item># reference nucleotides

    <item># number of minor (alternative) nucleotides

    <item>average quality of the reference nucleotides

    <item>average quality of the alternative nucleotides

    <item>first and second most frequent nucleotides in the pileup

    <item><with|mode|math|1-p(0)> where <math|p(f)> is the probability
    distribution function for the minor allele freqeuncy

    <item><with|mode|math|p(1)>

    <item><math|E(f)> mean value of <math|f>
  </enumerate>

  In addition, if <with|font-family|tt|-spectrum> is specified on the command
  line, the full pdf for <math|f> is printed after the fields listed above.

  <section|example>

  A typical command line:

  <with|font-base-size|8|.<with|font-family|tt|<with|font-base-size|9|/snape-pooled
  -nchr 9 -theta 0.1 -D 0.1 -priortype flat -fold folded \<less\>
  input_file.pool>><with|font-base-size|9|<with|font-base-size|10|>>>

  \;

  \;
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
    <associate|page-type|a4>
  </collection>
</initial>