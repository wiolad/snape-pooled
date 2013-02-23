/*
Copyright (C) 2011 by Emanuele Raineri
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#define MAXLINE 50000
#define MAXFIELDS 15


void splash(){
    fprintf (stderr, "***************************************************************\n"); 
    fprintf (stderr, "snape-pooled : a method for calling SNPs in pooled samples\n");
    fprintf (stderr, "$Date$ $Rev$\n");
    fprintf (stderr, "***************************************************************\n");
}


void parse_pileup(char ref, char* pileup, char* qualities, int* parsed){
	/* parsed contains the following fields:
	 0 number (#) of '.' or ',' characters (reference)
	 1 #A
	 2 #C
	 3 #G
	 4 #T	
	 5 #N
	 6 qA
	 7 qC
	 8 qG
	 9 qT			
	 10 sum of the quality codes for reference symbols
	 11 sum of quality codes for non reference symbols
	 12 the first 8 LSB bits indicate in which strand each character has been observed 
	 MSB TGCA(R)TGCA(F) LSB
	     0000   0000
	 */
    ref=toupper(ref); 
	int le = strlen(pileup);
	if (0==le) {fprintf(stderr,"empty pileup\n"); exit(1);}
	int i=0,pos=0; /* i is the index along the pileup, pos is the index along the qualities*/			
	while (i<le){
		char c=pileup[i];
		switch(c){
			case 'A':
				parsed[1]++;parsed[6]+=qualities[pos++];parsed[12]|=1;break;
			case 'a':
				parsed[1]++;parsed[6]+=qualities[pos++];parsed[12]|=16;break;
			case 'C':
				parsed[2]++;parsed[7]+=qualities[pos++];parsed[12]|=2;break;
			case 'c':
				parsed[2]++;parsed[7]+=qualities[pos++];parsed[12]|=32;break;
			case 'G':
				parsed[3]++;parsed[8]+=qualities[pos++];parsed[12]|=4;break;
			case 'g':
				parsed[3]++;parsed[8]+=qualities[pos++];parsed[12]|=64;break; 
			case 'T':
				parsed[4]++;parsed[9]+=qualities[pos++];parsed[12]|=8;break;
			case 't':
				parsed[4]++;parsed[9]+=qualities[pos++];parsed[12]|=128;break; 
			case 'N':
				parsed[5]++;parsed[11]+=qualities[pos++];break;
			case '.':
				switch(ref){
					case 'A':parsed[1]++;parsed[6]+=qualities[pos];parsed[12]|=1; break;
					case 'C':parsed[2]++;parsed[7]+=qualities[pos];parsed[12]|=2; break;
					case 'G':parsed[3]++;parsed[8]+=qualities[pos];parsed[12]|=4; break;
					case 'T':parsed[4]++;parsed[9]+=qualities[pos];parsed[12]|=8; break;
					case 'N':parsed[5]++;parsed[11]+=qualities[pos]; break;
					default: fprintf(stderr,"unknown ref %c\n",ref);exit(1);
				}
				parsed[0]++;parsed[10]+=qualities[pos++];break;
			case ',':
				switch(ref){
					case 'A':parsed[1]++;parsed[6]+=qualities[pos];parsed[12]|=16; break;
					case 'C':parsed[2]++;parsed[7]+=qualities[pos];parsed[12]|=32; break;
					case 'G':parsed[3]++;parsed[8]+=qualities[pos];parsed[12]|=64; break;
					case 'T':parsed[4]++;parsed[9]+=qualities[pos];parsed[12]|=128; break;
					case 'N':parsed[5]++;parsed[11]+=qualities[pos]; break;
					default: exit(1);
				}
				parsed[0]++; parsed[10]+=qualities[pos++];break;
			case '$':
				i++; continue;
			case '^':
				i+=2; continue;
			case '<':
			case '>': i++; pos++; continue;	
		}
		i++; /* position along the pileup advances anyway */		
	}
}

int split_line(char* line,char* fields[MAXFIELDS]){
    int ll=strlen(line);
    int i;
    int ncols=0,rl=0,rs=0; //running length running spaces
    char c;
    for(i=0;i<ll;i++){
        c = line[i];
        if (isspace(c)){
            if (rs++==0) {fields[ncols++][rl]='\0';rl=0; continue;}
            continue;
        } else {
            fields[ncols][rl++]=c;
            rs=0;
        }
    }
    return ncols;
}


int main(int argc, char* argv[]){
    splash();
    int c;
    int option_index = 0;
    int nchr;
    float theta;
    float bigd;
    char* priortype=malloc(128);
    char* fold=malloc(128);
    int spectrum=0;
    while(1){
        static struct option long_options[] =
        {
        {"nchr",  required_argument, 0, 'a'},
        {"theta",  required_argument, 0, 'b'},
        {"D", required_argument, 0, 'c'},
        {"priortype", required_argument, 0, 'd'},
        {"fold", required_argument, 0, 'e'},
        {0, 0, 0, 0}
        };
         c=getopt_long (argc, argv, "a:b:c:d:e:f:g:", long_options, &option_index);
        if (c == -1)
             break;
        switch (c){
            case 'a':
                nchr=atoi(optarg);
                fprintf(stderr,"nchr:%d\n",nchr);
                break;               
            case 'b':
                theta=atof(optarg);
                fprintf(stderr,"theta: %f\n",theta);
                break;               
            case 'c':
                bigd=atof(optarg);
                fprintf(stderr,"D:%f \n",bigd);
                break;               
            case 'd':
                priortype=strcpy(priortype,optarg);
                fprintf(stderr,"priortype:%s\n",priortype);
                break;               
            case 'e':
                fold=strcpy(fold,optarg);
                fprintf(stderr,"fold:%s\n",fold);
                break;               
       }     
    } /* end of command line parsing */
    /* chr position ref g pileup qualities */
    char chr[10];
    int pos;
    char ref;
    int g;
    char pileup[1000];
    char qualities[1000];
    int* parsed=malloc(sizeof(int)*13);
    int i;
    while(1){
        fscanf(stdin,"%s\t%d\t%c\t%d\t%s\t%s",chr,&pos,&ref,&g,pileup,qualities);
        if (feof(stdin)) break;
        fprintf(stderr,"chr=%s\tpos=%d\tref=%c\tg=%d\tpileup=%s\tq=%s\n",chr,pos,ref,g,pileup,qualities); 
        for(i=0;i<13;i++) parsed[i]=0;
        if (ref=='*' || 
		strchr(pileup,'+') != NULL ||
        strchr(pileup,'-') != NULL ||
        strchr(pileup,'<') != NULL ||
        strchr(pileup,'>') != NULL)  {
            fprintf(stderr,"skipping line %s:%d\n",chr,pos);
            continue;
        }
        parse_pileup(ref, pileup, qualities, parsed);
        for(i=0;i<13;i++) fprintf(stderr,"%d ",parsed[i]);  
        fprintf(stderr,"\n");
        int na=parsed[1]+parsed[2]+parsed[3]+parsed[4]+parsed[5];
        int nr=parsed[0];
    } 
    return 0;
}     



