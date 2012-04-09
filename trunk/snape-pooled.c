/*
Copyright (C) 2011 by Emanuele Raineri
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>

float q_sanger[94]={
1,            /* 33 ! */
0.794328,     /* 34 " */
0.630957,     /* 35 # */
0.501187,     /* 36 $ */
0.398107,     /* 37 % */
0.316228,     /* 38 & */
0.251189,     /* 39 ' */
0.199526,     /* 40 ( */
0.158489,     /* 41 ) */
0.125893,     /* 42 * */
0.1,          /* 43 + */
0.0794328,    /* 44 , */
0.0630957,    /* 45 - */
0.0501187,    /* 46 . */
0.0398107,    /* 47 / */
0.0316228,    /* 48 0 */
0.0251189,    /* 49 1 */
0.0199526,    /* 50 2 */
0.0158489,    /* 51 3 */
0.0125893,    /* 52 4 */
0.01,         /* 53 5 */
0.00794328,   /* 54 6 */
0.00630957,   /* 55 7 */
0.00501187,   /* 56 8 */
0.00398107,   /* 57 9 */
0.00316228,   /* 58 : */
0.00251189,   /* 59 ; */
0.00199526,   /* 60 < */
0.00158489,   /* 61 = */
0.00125893,   /* 62 > */
0.001,        /* 63 ? */
0.000794328,  /* 64 @ */
0.000630957,  /* 65 A */
0.000501187,  /* 66 B */
0.000398107,  /* 67 C */
0.000316228,  /* 68 D */
0.000251189,  /* 69 E */
0.000199526,  /* 70 F */
0.000158489,  /* 71 G */
0.000125893,  /* 72 H */
0.0001,       /* 73 I */
7.94328e-05,  /* 74 J */
6.30957e-05,  /* 75 K */
5.01187e-05,  /* 76 L */
3.98107e-05,  /* 77 M */
3.16228e-05,  /* 78 N */
2.51189e-05,  /* 79 O */
1.99526e-05,  /* 80 P */
1.58489e-05,  /* 81 Q */
1.25893e-05,  /* 82 R */
1e-05,        /* 83 S */
7.94328e-06,  /* 84 T */
6.30957e-06,  /* 85 U */
5.01187e-06,  /* 86 V */
3.98107e-06,  /* 87 W */
3.16228e-06,  /* 88 X */
2.51189e-06,  /* 89 Y */
1.99526e-06,  /* 90 Z */
1.58489e-06,  /* 91 [ */
1.25893e-06,  /* 92 \ */
1e-06,        /* 93 ] */
7.94328e-07,  /* 94 ^ */
6.30957e-07,  /* 95 _ */
5.01187e-07,  /* 96 ` */
3.98107e-07,  /* 97 a */
3.16228e-07,  /* 98 b */
2.51189e-07,  /* 99 c */
1.99526e-07,  /* 100 d */
1.58489e-07,  /* 101 e */
1.25893e-07,  /* 102 f */
1e-07,        /* 103 g */
7.94328e-08,  /* 104 h */
6.30957e-08,  /* 105 i */
5.01187e-08,  /* 106 j */
3.98107e-08,  /* 107 k */
3.16228e-08,  /* 108 l */
2.51189e-08,  /* 109 m */
1.99526e-08,  /* 110 n */
1.58489e-08,  /* 111 o */
1.25893e-08,  /* 112 p */
1e-08,        /* 113 q */
7.94328e-09,  /* 114 r */
6.30957e-09,  /* 115 s */
5.01187e-09,  /* 116 t */
3.98107e-09,  /* 117 u */
3.16228e-09,  /* 118 v */
2.51189e-09,  /* 119 w */
1.99526e-09,  /* 120 x */
1.58489e-09,  /* 121 y */
1.25893e-09,  /* 122 z */
1e-09,        /* 123 { */
7.94328e-10,  /* 124 | */
6.30957e-10,  /* 125 } */
5.01187e-10,  /* 126 ~ */
};

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



