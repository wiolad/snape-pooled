#!/bin/bash -x

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=12000M
#SBATCH --partition=general

date
source /opt/Modules/3.2.9/init/Modules4bash.sh
module load bwa-0.6.1
module load samtools-0.1.18-sl61
module load perl-libs-5.10
UTILSPATH=/projects/039-snape/bin

######################  RUN OPTIONS (input/output)  ###############################
#SBATCH --job-name=N100D20R100
DIRDATA=/projects/039-snape/N100D20R100      #  working folder
PREFIX=N100D20R100 # output prefix

REFSEQ=testseq1m_ref.fa  # original (i.e. "ancestral") fasta sequence to use in simulations. length should match next line
SEQLEN=1000000  # nr bps
PATHTOREFSEQ=/projects/039-snape/input # path to sequence


######################  RUN OPTIONS (parameters)  ###############################
NUMIND=100        # number of individuals (haploids) in the pool. This should match the second pop size in the MSCMD above
NUMREPS=100       # number of replicates to run....
AVERAGEDEPTH=20	 # average depth in ART, for each individual
MINDEP=5        # minimum depth for filtering
MAXDEP=40       # maximum depth for filtering

############ YET MORE OPTIONS (probably unchanged for most runs) ################
MSCMD=''$UTILSPATH/'ms-mod_'$NUMIND'_1_-t_500_-r_500_1000000'   	# ms command line to simulate sequences (use underscores instead of spaces)
MAXTRIES=10000  # maximum number of iterations in ms until a "good" dataset is found
                # by good, i mean a dataset without multiple hits. This is to prevent the pipeline 
                #  from getting stuck if the MSCMD is not OK
NMP=10          # number of mismatches permitted in bwa aln
MAPQ=20         # min map quality
NP=4            # no. of processors for alignment step
MINVARFREQ=$(echo "scale=2; 0.5/$NUMIND" | bc) # minimum frequency for alternative allele
STEPSIZE=1000000 # stepsize to use for popoolation
WINSIZE=1000000  # window size to use for popoolation
MINCOVFRACT=0.000001 # another parameter value used for popoolation
SNPQ=20          # minimum quality for popoolation
MINCOUNT=2       # minimum count of the minor allele for popoolation
THETA=0.001      # theta for snape
PD=0.01           # Prior genetic difference between reference genome and population (snape)
PF=0.9           # 1 − p(0) (for snape, value in 9th column)

################################  PIPELINE  #####################################

# Working directory
cd $DIRDATA
# move final results file from previous run, in case it is present
find . -maxdepth 1 -name "summary.$PREFIX.txt" | xargs -i mv {} summary.$PREFIX.txt.old
find . -maxdepth 1 -name "summary2.$PREFIX.txt" | xargs -i mv {} summary2.$PREFIX.txt.old
# clean possible garbage from previous (failed/aborted) runs
rm *.out.ms.txt *.out.all.fas *.reads1.fq *.reads2.fq
# Index sequence (is fast, even in case it has been indexed already in previous runs)
# first copy the sequence into working directory, to avoid permission problems
cp $PATHTOREFSEQ/$REFSEQ $DIRDATA
bwa index -a is $DIRDATA/$REFSEQ


# start loop for replicates
for ((REP = 1; REP <= $NUMREPS; REP++ ))
do
##### RUN MS ##### recursively until a dataset is found without multiple hits
# first create a seedms file with random numbers for seeding ms random generator
$UTILSPATH/seedMS.pl > seedms
# now run ms from perl script
STARTTIME=$(date +%s)
$UTILSPATH/msExplorer-diallelic.pl  $MSCMD $SEQLEN $PREFIX.out.ms.txt $MAXTRIES
ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds to run msExplorer-diallelic.pl"
STARTTIME=$(date +%s)
# ms file written to $prefix.out.ms.txt. now we make a fasta file out of it, using the reference sequence
$UTILSPATH/msToFas.pl $PREFIX.out.ms.txt $PATHTOREFSEQ/$REFSEQ > $PREFIX.out.all.fas
ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds to run msToFas.pl"
STARTTIME=$(date +%s)

# no longer needed:
# extract the first sequence from the fasta file. this is the reference "genome", so we can index it now too
# sed -n '1,2p' $PREFIX.out.all.fas > $PREFIX.replicate$REP.reference.fas
# bwa index -a is $PREFIX.replicate$REP.reference.fas

# start loop for each individual (to run ART)
for ((IND = 1; IND <= $NUMIND; IND++ ))
do
INDNAME=$(((($IND-1)*2)+1)) # writing this in shell is a pain. INDNAME is an index, for the line of $PREFIX.out.all.fas,
INDSEQ=$(((($IND-1)*2)+2))  # where the name of current individual is. INDSEQ is the index for its sequence. 

# extract current individual from $PREFIX.out.all.fas
SEDCMD='sed -n "'${INDNAME},${INDSEQ}'p" '$PREFIX.out.all.fas' > '$PREFIX.currentInd.fas''
eval $SEDCMD

##### RUN ART ##### (paired end, illumina) with the average depth for each individual sequence
$UTILSPATH/art_illumina -i $PREFIX.currentInd.fas -p -l 75 -f $AVERAGEDEPTH -m 500 -s 10 -o $DIRDATA/$PREFIX.replicate$REP.individual$IND.
# finish ART loop for each individual
done
ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds to run the ART step"
STARTTIME=$(date +%s)


##### MERGE FQ FILES (make the pool) #####
# for now, we will sample equally from each individual (this can be changed later on by creating a different PROPORTIONS array)
# so here create an array with number of entries=number of individuals, and values the proportion to sample from each ind
PROPORTIONS[1]=1
for ((i = 1; i <= $NUMIND; i++ ))
do
PROPORTIONS[$i]=$(echo "scale=2; 1 / $NUMIND" | bc)
# again these weird arithmetic operators from shell scripts... but it works it seems, it just divides 1 by the number of inds
done
ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds to prepare PROPORTIONS"
STARTTIME=$(date +%s)


for ((IND = 1; IND <= $NUMIND; IND++ ))
do
# run perl script to subsample reads from each .fq file
$UTILSPATH/subsampleReads.pl ${PROPORTIONS[$IND]} $PREFIX.replicate$REP.individual$IND.1.fq $PREFIX.replicate$REP.individual$IND.2.fq Individual
cat $PREFIX.replicate$REP.individual$IND.1.fq.subsample.fq >> $PREFIX.replicate$REP.reads1.fq
cat $PREFIX.replicate$REP.individual$IND.2.fq.subsample.fq >> $PREFIX.replicate$REP.reads2.fq
done
##### FINISHED MERGING FQ FILES INTO POOL #####
ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds to make the pooled reads"
STARTTIME=$(date +%s)

# align reads
bwa aln -n $NMP -t $NP $PATHTOREFSEQ/$REFSEQ $PREFIX.replicate$REP.reads1.fq > $PREFIX.replicate$REP.1.sai
bwa aln -n $NMP -t $NP $PATHTOREFSEQ/$REFSEQ $PREFIX.replicate$REP.reads2.fq > $PREFIX.replicate$REP.2.sai
# generate alignment in BAM format
bwa sampe -P $PATHTOREFSEQ/$REFSEQ $PREFIX.replicate$REP.1.sai $PREFIX.replicate$REP.2.sai $PREFIX.replicate$REP.reads1.fq $PREFIX.replicate$REP.reads2.fq | samtools view -q $MAPQ -Suh - | samtools sort -m 5000000000 - $PREFIX.replicate$REP

# Extracts coverage per site
samtools view -q $MAPQ -u $PREFIX.replicate$REP.bam  | samtools mpileup -BQ0 -d 100000 -f $PATHTOREFSEQ/$REFSEQ - | awk '{print $2 " " $4}' > $PREFIX.replicate$REP.coverages
ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds to align into bam file and extract coverages"
STARTTIME=$(date +%s)

#################   SNP CALLING ######################
# first do the mpileup file, which is needed for other programs as well. filter mpileup by depth already
samtools mpileup -f $PATHTOREFSEQ/$REFSEQ $PREFIX.replicate$REP.bam | awk '$4>='$MINDEP' && $4<='$MAXDEP'' > $PREFIX.replicate$REP.flt.mpileup
ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds to get flt.mpileup file from samtools"
TIMEMPILEUP=$DIFF
STARTTIME=$(date +%s)

# with Samtools, filtering by minimum and maximum depth. also removing indels and homozygous positions
samtools mpileup -uf $PATHTOREFSEQ/$REFSEQ $PREFIX.replicate$REP.bam | bcftools view -bvcg - > $PREFIX.replicate$REP.bcf 
bcftools view $PREFIX.replicate$REP.bcf | vcfutils.pl varFilter -d $MINDEP -D $MAXDEP | grep -v "INDEL"|grep -v "1/1" - > $PREFIX.replicate$REP.samtools.flt.vcf
# now create coverages file for samtools (exclude homozygous SNPS and indels from .coverages file)
#bcftools view $PREFIX.replicate$REP.bcf | vcfutils.pl varFilter  -d $MINDEP -D $MAXDEP | grep -v "#" | grep "1/1\|INDEL" | awk '{print $2}' >  $PREFIX.replicate$REP.samtools.toRemove
#awk 'FNR==NR {a[$1];next} !($1 in a) '  $PREFIX.replicate$REP.samtools.toRemove $PREFIX.replicate$REP.coverages > $PREFIX.replicate$REP.coverages.samtools

ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds in samtools mpileup SNP calling"
STARTTIME=$(date +%s)

# with old Samtools (pileup).  also removing indels and homozygous positions
/opt/samtools-0.1.16/samtools view -q $MAPQ -u $PREFIX.replicate$REP.bam | /opt/samtools-0.1.16/samtools pileup -vcf $PATHTOREFSEQ/$REFSEQ - > $PREFIX.replicate$REP.oldSam.vcf
samtools.pl varFilter -d $MINDEP -D $MAXDEP -S $SNPQ $PREFIX.replicate$REP.oldSam.vcf | awk '$4~/[WRSYMK]/' > $PREFIX.replicate$REP.oldSamtools.flt.vcf
# now create coverages file for samtoolsold (exclude homozygous SNPS and indels from .coverages file)
#samtools.pl varFilter -d $MINDEP -D $MAXDEP -S $SNPQ $PREFIX.replicate$REP.oldSam.vcf | awk '$4!~/[WRSYMK]/'   | awk '{print $2}' >  $PREFIX.replicate$REP.samtoolsold.toRemove
#awk 'FNR==NR {a[$1];next} !($1 in a) '  $PREFIX.replicate$REP.samtoolsold.toRemove $PREFIX.replicate$REP.coverages > $PREFIX.replicate$REP.coverages.samtoolsold


ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds in samtools pileup SNP calling"
STARTTIME=$(date +%s)

# with VarScan, only filtering by min depth, but mpileup file is already filtered. also removing indels and homozygous positions
java -jar $UTILSPATH/VarScan.v2.2.8.jar mpileup2snp $PREFIX.replicate$REP.flt.mpileup --min-coverage $MINDEP --min-var-freq $MINVARFREQ --p-value 0.1 --output-vcf | grep -v "INDEL"|grep -v "1/1" - > $PREFIX.replicate$REP.varscan.flt.vcf
# now create coverages file for varscan (exclude homozygous SNPS and indels from .coverages file)
#java -jar $UTILSPATH/VarScan.v2.2.8.jar mpileup2snp $PREFIX.replicate$REP.flt.mpileup --min-coverage $MINDEP --min-var-freq $MINVARFREQ --p-value 1 --output-vcf | grep "1/1" | awk '{print $2}'  > $PREFIX.replicate$REP.varscan.toRemove
#awk 'FNR==NR {a[$1];next} !($1 in a) '  $PREFIX.replicate$REP.varscan.toRemove $PREFIX.replicate$REP.coverages > $PREFIX.replicate$REP.coverages.varscan


ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds in varScan SNP calling (+ $TIMEMPILEUP for flt.mpileup file)"
STARTTIME=$(date +%s)

# with popoolation. remove INDELS and HOMO SNPS
#samtools view -h $PREFIX.replicate$REP.bam > $PREFIX.replicate$REP.sam
#perl $UTILSPATH/popoolation_1.2.2/basic-pipeline/mask-sam-indelregions.pl --input $PREFIX.replicate$REP.sam --output $PREFIX.replicate$REP.indelmasked.sam
#samtools view -Sb $PREFIX.replicate$REP.indelmasked.sam > $PREFIX.replicate$REP.indelmasked.bam
#samtools mpileup -f $PATHTOREFSEQ/$REFSEQ $PREFIX.replicate$REP.indelmasked.bam | awk '$4>='$MINDEP' && $4<='$MAXDEP''> $PREFIX.replicate$REP.indelmasked.flt.mpileup

/opt/samtools-0.1.16/samtools pileup -f $PATHTOREFSEQ/$REFSEQ $PREFIX.replicate$REP.bam | awk '$4>='$MINDEP' && $4<='$MAXDEP''> $PREFIX.replicate$REP.flt.pileup
perl $UTILSPATH/popoolation_1.2.2/Variance-sliding.pl  --measure pi --input $PREFIX.replicate$REP.flt.pileup  --min-count $MINCOUNT --min-qual $SNPQ --min-coverage $MINDEP --max-coverage $MAXDEP --pool-size $NUMIND --snp-output $PREFIX.replicate$REP.SNP.output.temp --min-covered-fraction $MINCOVFRACT --window-size $WINSIZE --step-size $STEPSIZE --output  $PREFIX.replicate$REP.window --fastq-type "sanger" 
cat $PREFIX.replicate$REP.SNP.output.temp | awk '$9==0 && $5!=$4 && $6!=$4 && $7!=$4 && $8!=$4' > $PREFIX.replicate$REP.SNP.output
# now create coverages file for popoolation (exclude homozygous SNPS and indels from .coverages file)
#cat $PREFIX.replicate$REP.SNP.indelmasked.output.temp | grep -v "^>" | awk '$9==1 || $5==$4 || $6==$4 || $7==$4 || $8==$4' > $PREFIX.replicate$REP.popoolation.toRemove
#awk 'FNR==NR {a[$1];next} !($1 in a) '  $PREFIX.replicate$REP.popoolation.toRemove $PREFIX.replicate$REP.coverages > $PREFIX.replicate$REP.coverages.popoolation

ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds in popoolation SNP calling"
STARTTIME=$(date +%s)

# with snape informative. removes homozygous positions
cat $PREFIX.replicate$REP.flt.mpileup | $UTILSPATH/snape-pooled.linux.r17 -nchr $NUMIND -theta $THETA -D $PD -fold folded -priortype informative 2> $PREFIX.replicate$REP.snapeinfo.log | awk '$3!="N" && $3!=$8 && $9>'$PF' && ($4+$5)>='$MINDEP' && ($4+$5)<='$MAXDEP' && $8!="A" && $8!="T" && $8!="C" && $8!="G" && $9!="-nan"' > $PREFIX.replicate$REP.snapeinfo.vcf
ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds in snape informative SNP calling (+ $TIMEMPILEUP for flt.mpileup file)"
STARTTIME=$(date +%s)

# with snape flat.  removes homozygous positions
cat $PREFIX.replicate$REP.flt.mpileup | $UTILSPATH/snape-pooled.linux.r17 -nchr $NUMIND -theta $THETA -D $PD -fold folded -priortype flat 2> $PREFIX.replicate$REP.snapeflat.log | awk '$3!="N" && $3!=$8 && $9>'$PF' && ($4+$5)>='$MINDEP' && ($4+$5)<='$MAXDEP' && $8!="A" && $8!="T" && $8!="C" && $8!="G" && $9!="-nan"' > $PREFIX.replicate$REP.snapeflat.vcf
ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds in snape flat SNP calling (+ $TIMEMPILEUP for flt.mpileup file)"
STARTTIME=$(date +%s)





#################   CALCULATE FDR and POWER ###################### 
# snape informative
$UTILSPATH/sum_vcf2.pl $PREFIX.out.ms.txt  $SEQLEN SNAPEINFO.rep$REP $MINDEP $MAXDEP $PREFIX.replicate$REP.coverages $PREFIX.replicate$REP.snapeinfo.vcf >> summary.$PREFIX.txt 
# with snape flat
$UTILSPATH/sum_vcf2.pl $PREFIX.out.ms.txt  $SEQLEN SNAPEFLAT.rep$REP $MINDEP $MAXDEP $PREFIX.replicate$REP.coverages $PREFIX.replicate$REP.snapeflat.vcf >> summary.$PREFIX.txt 
# with varscan
$UTILSPATH/sum_vcf2.pl $PREFIX.out.ms.txt  $SEQLEN VARSCAN.rep$REP $MINDEP $MAXDEP $PREFIX.replicate$REP.coverages $PREFIX.replicate$REP.varscan.flt.vcf >> summary.$PREFIX.txt
#$UTILSPATH/sum_vcf2.pl $PREFIX.out.ms.txt  $SEQLEN VARSCAN.2.rep$REP $MINDEP $MAXDEP $PREFIX.replicate$REP.coverages.varscan $PREFIX.replicate$REP.varscan.flt.vcf >> summary.$PREFIX.txt
# with popoolation
$UTILSPATH/sum_vcf2.pl $PREFIX.out.ms.txt  $SEQLEN POPOOLATION.rep$REP $MINDEP $MAXDEP $PREFIX.replicate$REP.coverages $PREFIX.replicate$REP.SNP.output >> summary.$PREFIX.txt
#$UTILSPATH/sum_vcf2.pl $PREFIX.out.ms.txt  $SEQLEN POPOOLATION.2.rep$REP $MINDEP $MAXDEP $PREFIX.replicate$REP.coverages.popoolation $PREFIX.replicate$REP.SNP.indelmasked.output >> summary.$PREFIX.txt
# with samtools
$UTILSPATH/sum_vcf2.pl $PREFIX.out.ms.txt  $SEQLEN SAMTOOLS.rep$REP $MINDEP $MAXDEP $PREFIX.replicate$REP.coverages $PREFIX.replicate$REP.samtools.flt.vcf >> summary.$PREFIX.txt
#$UTILSPATH/sum_vcf2.pl $PREFIX.out.ms.txt  $SEQLEN SAMTOOLS.2.rep$REP $MINDEP $MAXDEP $PREFIX.replicate$REP.coverages.samtools $PREFIX.replicate$REP.samtools.flt.vcf >> summary.$PREFIX.txt
# with samtools old (pileup)
$UTILSPATH/sum_vcf2.pl $PREFIX.out.ms.txt  $SEQLEN SAMTOOLSold.rep$REP $MINDEP $MAXDEP $PREFIX.replicate$REP.coverages $PREFIX.replicate$REP.oldSamtools.flt.vcf >> summary.$PREFIX.txt
#$UTILSPATH/sum_vcf2.pl $PREFIX.out.ms.txt  $SEQLEN SAMTOOLSold.2.rep$REP $MINDEP $MAXDEP $PREFIX.replicate$REP.coverages.samtoolsold $PREFIX.replicate$REP.oldSamtools.flt.vcf >> summary.$PREFIX.txt
ENDTIME=$(date +%s)
DIFF=$(( $ENDTIME - $STARTTIME ))
echo "Time: It took $DIFF seconds to summarise all SNP callings (sum_vcf2.pl)"
STARTTIME=$(date +%s)


# clean before next replicate
rm $PREFIX.*

# finish loop for replicates
done
date
exit 0

