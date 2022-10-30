
/*****************************************************************************/
/***  (seg.c)                                                              ***/
/*** #include precomputed ln(fact(n)) array in "lnfac.h"                   ***/
/*****************************************************************************/

/*--------------------------------------------------------------(includes)---*/

#include "genwin.h"
#include "lnfac.h"
#include "ln.h"

/*---------------------------------------------------------------(defines)---*/

#define LNFACMAX 10000
#define LENCORRLIM 120
#define MIN(a,b)	((a) <= (b) ? (a) : (b))

/*---------------------------------------------------------------(structs)---*/

struct Segment
  {
   int begin;
   int end;
   struct Segment *next;
  };

/*---------------------------------------------------------------(globals)---*/

int window = 12;
int downset, upset;
double locut = 2.2;
double hicut = 2.5;

int hilenmin = 0;
int overlaps = FALSE;
int hionly = FALSE;
int loonly = FALSE;
int entinfo = TRUE;
int singleseq = FALSE;
int prettyseq = FALSE;
int prettytree = TRUE;
int debug = FALSE;
int charline = 60;
int maxtrim = 100;
char* alphabet = "prot";

double *_lnfacw[LNFACMAX];
void lnfacwinit();
void lnfacwfree();
double getprob(), lnperm(), lnass(), lnfacw();

/*------------------------------------------------------------------(main)---*/


main(argc, argv)
  int argc;
  char *argv[];

  {struct Database *db;
   struct Sequence *seq;
   struct Segment *segs;
   int ctime;


      //fprintf(stderr, "main\n");

/* readlencorr(); */                        /* #include lencorr file */

   getparams(argc, argv);
   genwininit(alphabet);
   lnfacwinit(abet->weight);
   entropy_init(window);

   if ((db=opendbase(argv[1]))==NULL)
     {
      fprintf(stderr, "Error opening file %s\n", argv[1]);
      exit(1);
     }

   for (seq=firstseq(db); seq!=NULL; seq=nextseq(db))
     {
      segs = (struct Segment *) NULL;
      segseq(seq, &segs, 0);

      mergesegs(seq, segs);

      if (singleseq) singreport(seq, segs);
      else if (prettyseq) prettyreport(seq, segs);
      else if (prettytree) pretreereport(seq, segs);
      else report(seq, segs);

      freesegs(segs);
      closeseq(seq);
     }

   closedbase(db);
   genwinterm();
   lnfacwfree();
   exit(0);
  }

/*-------------------------------------------------------------(getparams)---*/

getparams(argc, argv)
  int argc;
  char *argv[];

  {int i;
   int nargc;
   char **nargv;
   extern char *optarg;
   extern int optind;
   int c;


      //fprintf(stderr, "getparams\n");

   if (argc<2)
     {
      usage();
      exit(1);
     }

   for (i=2; argc>i && argv[i][0]!='-'; i++)
     {
      if (i==2)
        {
         window = atoi(argv[2]);
        }
      else if (i==3)
        {
         locut = atof(argv[3]);
         hicut = locut + 0.3;
        }
      else if (i==4)
        {
         hicut = atof(argv[4]);
        }
      else if (i==5){
          alphabet = argv[5];
      }
     }



   if (locut>hicut)
     {
      fprintf(stderr, "Warning: trigger complexity greater than extension\n");
      hicut = locut;
     }


   downset = (window+1)/2 - 1;
   upset = window - downset;

   if (i==argc) return;

   nargc = argc-i+1;
   nargv = argv+(i-1);
   while ((c=getopt(nargc, nargv, "dm:olhaxpqc:nt:"))!=-1)
     {
      switch(c)
        {
         case 'd':
            debug = TRUE;
            break;
         case 'm':
            hilenmin = atoi(optarg);
            break;
         case 'o':
            overlaps = TRUE;
            hilenmin = 0;
            break;
         case 'l':
            loonly = TRUE;
            singleseq = FALSE;
            prettyseq = FALSE;
            prettytree = FALSE;
            break;
         case 'h':
            hionly = TRUE;
            singleseq = FALSE;
            prettyseq = FALSE;
            prettytree = FALSE;
            break;
         case 'a':
            hionly = FALSE;
            loonly = FALSE;
            singleseq = FALSE;
            prettyseq = FALSE;
            prettytree = FALSE;
            break;
         case 'x':
            singleseq = TRUE;
            prettyseq = FALSE;
            prettytree = FALSE;
            hionly = TRUE;
            loonly = FALSE;
            break;
         case 'p':
            prettytree = TRUE;
            prettyseq = FALSE;
            singleseq = FALSE;
            hionly = FALSE;
            loonly = FALSE;
            break;
         case 'q':
            prettyseq = TRUE;
            prettytree = FALSE;
            singleseq = FALSE;
            hionly = FALSE;
            loonly = FALSE;
            break;
         case 'c':
            charline = atoi(optarg);
            break;
         case 'n':
            entinfo = FALSE;
            break;
         case 't':
            maxtrim = atoi(optarg);
            break;
         case '?':
            fprintf(stderr, "Unknown option.\n");
            usage();
            exit(1);
            break;
        }
     }   

   if(debug){fprintf(stderr, "Param: %d %f %f %s\n", window, locut, hicut, alphabet);}

   return;
  }

/*---------------------------------------------------------------(segment)---*/

segseq(seq, segs, offset)
  struct Sequence *seq;
  struct Segment **segs;
  int offset;

  {struct Segment *seg, *leftsegs;
   struct Sequence *leftseq;
   int first, last, lowlim;
   int loi, hii, i;
   int leftend, rightend, lend, rend;
   double *H, *seqent();


      if(debug){fprintf(stderr, "segseq(%s)\n",seq->id);}


   H = seqent(seq);


   if (H==NULL) return;

   first = downset;
   last = seq->length - upset;
   lowlim = first;

   if(debug){
       for(i=0;i<last; i++){
           fprintf(stderr,"%0.2f ",H[i]);
       }
       fprintf(stderr,"\n");
   }

   for (i=first; i<=last; i++)
     {
        if(debug) {fprintf(stderr, "check: %d, H (%0.2f)\n",i,H[i]);}
      if (H[i]<=locut && H[i]!=-1)
        {

         loi = findlo(i, lowlim, H);
         hii = findhi(i, last, H);

         leftend = loi - downset;
         rightend = hii + upset - 1;

         trim(openwin(seq, leftend, rightend-leftend+1), &leftend, &rightend);

         if (i+upset-1<leftend)   /* check for trigger window in left trim */
           {
            lend = loi - downset;
            rend = leftend - 1;

            leftseq = openwin(seq, lend, rend-lend+1);
            leftsegs = (struct Segment *) NULL;
            segseq(leftseq, &leftsegs, offset+lend);
            if (leftsegs!=NULL)
              {
               if (*segs==NULL) *segs = leftsegs;
               else appendseg(*segs, leftsegs);
              }
            closewin(leftseq);

/*          trim(openwin(seq, lend, rend-lend+1), &lend, &rend);
            seg = (struct Segment *) malloc(sizeof(struct Segment));
            seg->begin = lend;
            seg->end = rend;
            seg->next = (struct Segment *) NULL;
            if (segs==NULL) segs = seg;
            else appendseg(segs, seg);  */
           }

         seg = (struct Segment *) malloc(sizeof(struct Segment));
         seg->begin = leftend + offset;
         seg->end = rightend + offset;
         seg->next = (struct Segment *) NULL;

         if (*segs==NULL) *segs = seg;
         else appendseg(*segs, seg);

         i = min(hii, rightend+downset);
         lowlim = i + 1;
/*       i = hii;     this ignores the trimmed residues... */
        }
     }

   free(H);
   return;
  }

/*----------------------------------------------------------------(seqent)---*/

double *seqent(seq)
  struct Sequence *seq;

  {struct Sequence *win;
   double *H;
   int i, first, last;
   

      //fprintf(stderr, "segent\n");

   if (window>seq->length)
     {
      return((double *) NULL);
     }

   H = (double *) malloc(seq->length*sizeof(double));

   for (i=0; i<seq->length; i++)
     {
      H[i] = -1.;
     }

   win = openwin(seq, 0, window);
   enton(win);

   //fprintf(stderr,"comp: ");
   //int x;
   //for(x = 0; x < win->alphabet->size; x++){
   //    fprintf(stderr,"%d ",win->composition[x]);
   //}
   //fprintf(stderr,"\n");

   first = downset;
   last = seq->length - upset;

   for (i=first; i<=last; i++)
     {
      if (seq->punctuation && hasdash(win))
        {H[i] = -1;
         shiftwin1(win);
         continue;}

      H[i] = win->entropy;

      shiftwin1(win);

     }

   closewin(win);
   return(H);
  }

/*---------------------------------------------------------------(hasdash)---*/

hasdash(win)
  struct Sequence *win;
{
	register char	*seq, *seqmax;


      //fprintf(stderr, "hasdash\n");

	seq = win->seq;
	seqmax = seq + win->length;

	while (seq < seqmax) {
		if (*seq++ == '-')
			return TRUE;
	}
	return FALSE;
}

/*----------------------------------------------------------------(findlo)---*/

findlo(i, limit, H)
  int i, limit;
  double *H;

  {int j;


      //fprintf(stderr, "findlo\n");

   for (j=i; j>=limit; j--)
     {
      if (H[j]==-1) break;
      if (H[j]>hicut) break;
     }

   return(j+1);
  }

/*----------------------------------------------------------------(findhi)---*/

findhi(i, limit, H)
  int i, limit;
  double *H;

  {int j;


      //fprintf(stderr, "findhi\n");

   for (j=i; j<=limit; j++)
     {
      if (H[j]==-1) break;
      if (H[j]>hicut) break;
     }

   return(j-1);
  }

/*------------------------------------------------------------------(trim)---*/

trim(seq, leftend, rightend)
  struct Sequence *seq;
  int *leftend, *rightend;

  {struct Sequence *win;
   double prob, minprob;
   int shift, len, i;
   int lend, rend;
   int minlen;

   if(debug){
       fprintf(stderr, "(%d-%d) ", *leftend, *rightend);
       for(i = 0; i <= *rightend-*leftend; i++){
           fputc(*(seq->seq + i),stderr);
       }
       fprintf(stderr,"\n");
   }

   lend = 0;
   rend = seq->length - 1;
   minlen = 1;
   len=seq->length;
   if(len > LNFACMAX){
      fprintf(stderr, "[WARNING] (%s) Found expanded LCR in excess of %d bp long: Trimming\n",seq->id,LNFACMAX);

      len = LNFACMAX;
   }
   if ((len-maxtrim)>minlen) minlen = len-maxtrim;
   if(minlen < 1) minlen = 1;

   minprob = 1.;
   for (len; len>minlen; len--)
     {
      win = openwin(seq, 0, len);
      i = 0;

      shift = TRUE;
      while (shift)
        {
         prob = getprob(win->state, len,win->alphabet->size);
         if (prob<minprob)
           {
            minprob = prob;
            lend = i;
            rend = len + i - 1;
           }
         shift = shiftwin1(win);
         i++;
        }
      closewin(win);
     }

/* fprintf(stderr, "%d-%d ", *leftend, *rightend);  */

   *leftend = *leftend + lend;
   *rightend = *rightend - (seq->length - rend - 1);

/* fprintf(stderr, "%d-%d\n", *leftend, *rightend);  */

   closewin(seq);
   return;
  }

/*---------------------------------------------------------------(getprob)---*/

double getprob(sv, total, abetsize)
  int *sv;
  int total;
  int abetsize;

  {double ans, totseq;
      double lnf, lnomega;

   totseq = ((double) total) * ln[abetsize];

   lnf = lnass(sv,abetsize);
   lnomega = lnperm(sv, total);

   ans = lnf + lnomega - totseq;

   if(debug){
       fprintf(stderr, "getprob) total: %d, Ω: %0.2f, F: %0.2f, totseq: %0.2f, prob: %0.2f\n",total,lnomega,lnf,totseq,ans);
       int i;
       for(i = 0; i < abetsize; i++){
           fprintf(stderr,"%d ",sv[i]);
       }
       fprintf(stderr,"\n");
   }

   return(ans);
  }

/*----------------------------------------------------------------(lnperm)---*/

double lnperm(sv, tot)
  int *sv;
   int tot;

  {double ans;
   int i;

      //fprintf(stderr, "lnperm (sv, %d)\n", tot);

   ans = lnfac[tot];

   for (i=0; sv[i]!=0; i++) 
     {
      ans -= lnfacw(sv[i]);
     }

   return(ans);
  }

/*-----------------------------------------------------------------(lnass)---*/

double lnass(sv,abetsize)
	register int	*sv;
        register int    abetsize;
{
	double	ans;
	register int	svi, svim1;
	register int	class, total;
	register int    i;


      //fprintf(stderr, "lnass\n");

	ans = lnfac[abetsize];
	if (sv[0] == 0)
		return ans;

	total = abetsize;
	class = 1;
	svim1 = sv[0];
	for (i=0;; svim1 = svi) {
	        if (++i==abetsize) {
		        ans -= lnfac[class];
                        break;
		      }
		else if ((svi = *++sv) == svim1) {
			class++;
			continue;
		}
		else {
			total -= class;
			ans -= lnfac[class];
			if (svi == 0) {
				ans -= lnfac[total];
				break;
			}
			else {
				class = 1;
				continue;
			}
		}
	}

	return ans;
}

/*-------------------------------------------------------------(mergesegs)---*/

mergesegs(seq, segs)
  struct Sequence *seq;
  struct Segment *segs;

  {struct Segment *seg, *nextseg;
   int len;


      //fprintf(stderr, "mergeseqs\n");

   if (overlaps) return;
   if (segs==NULL) return;

   if (segs->begin<hilenmin) segs->begin = 0;

   seg = segs;
   nextseg = seg->next;

   while (nextseg!=NULL)
     {
      if (seg->end>=nextseg->begin)               /* overlapping segments */
        {
         seg->end = nextseg->end;
         seg->next = nextseg->next;
         free(nextseg);
         nextseg = seg->next;
         continue;
        }
      len = nextseg->begin - seg->end - 1;
      if (len<hilenmin)                            /* short hient segment */
        {
         seg->end = nextseg->end;
         seg->next = nextseg->next;
         free(nextseg);
         nextseg = seg->next;
         continue;
        }
      seg = nextseg;
      nextseg = seg->next;
     }

   len = seq->length - seg->end - 1;
   if (len<hilenmin) seg->end = seq->length - 1;

   return;
  }

/*----------------------------------------------------------------(report)---*/

report(seq, segs)
  struct Sequence *seq;
  struct Segment *segs;

  {struct Sequence *subseq;
   struct Segment *seg, *nextseg;
   static int hi = 1;
   static int lo = 0;


   if (segs==NULL)
     {
      enton(seq);
      seqout(seq, hi, 1, seq->length);
/*    fputc('\n', stdout);   -- for spacing after each sequence */
      return;
     }

   if (segs->begin>0)
     {
      subseq = openwin(seq, 0, segs->begin);
      enton(subseq);
      seqout(subseq, hi, 1, segs->begin);
      closewin(subseq);
     }

   for (seg=segs; seg!=NULL; seg=seg->next)
     {
      subseq = openwin(seq, seg->begin, seg->end-seg->begin+1);
      enton(subseq);
      seqout(subseq, lo, seg->begin+1, seg->end+1);
      closewin(subseq);

      if (seg->next==NULL)
        {
         break;
        }

      nextseg = seg->next;
      
      if (nextseg->begin<=seg->end)
        {
         fprintf(stderr, "Overlapping segments: %s\n", seq->id);
         continue;
        }

      if (nextseg->begin==seg->end+1)
        {
         continue;
        }

      subseq = openwin(seq, seg->end+1, nextseg->begin-seg->end-1);
      enton(subseq);
      seqout(subseq, hi, seg->end+2, nextseg->begin);
      closewin(subseq);
     }

   if (seg->end+1==seq->length)
     {
/*    fputc('\n', stdout);   -- for spacing after each sequence */
      return;
     }

   subseq = openwin(seq, seg->end+1, seq->length-seg->end-1);
   enton(subseq);
   seqout(subseq, hi, seg->end+2, seq->length);
   closewin(subseq);

/* fputc('\n', stdout);   -- for spacing after each sequence */
   return;
  }

/*------------------------------------------------------------(singreport)---*/

singreport(seq, segs)
	struct Sequence	*seq;
	struct Segment	*segs;
{
	char	*proseq, *proseqmax;
	struct Segment	*seg;
	int	begin, end, i, ctr;

	proseq = seq->seq;
	proseqmax = proseq + seq->length;
	upper(proseq, seq->length);

	for (seg=segs; seg!=NULL; seg=seg->next) {
		begin = seg->begin;
		end = seg->end;
		memset(proseq + begin, 'x', end - begin +1);
	}

	fprintf(stdout, "%s\n", seq->header);

	for (i=0, ctr=0; proseq < proseqmax; ++i, ++ctr, ++proseq) {
		if (ctr==charline) {
			putc('\n', stdout);
			ctr = 0;
		}
		putc(*proseq, stdout);
	}

	putc('\n', stdout);
	if (putc('\n', stdout) == EOF) {
		fprintf(stderr, "premature EOF on write\n");
		exit(2);
	}
}

/*----------------------------------------------------------(prettyreport)---*/

prettyreport(seq, segs)
  struct Sequence *seq;
  struct Segment *segs;

{
	char	*proseq, *proseqmax;
	char	format[10];
	struct Segment	*seg;
	int	begin, end, i, ctr;
	int	leftspace;

	leftspace = (int) ceil(log10((double)seq->length));
	sprintf(format, "%%%dd ", leftspace);

	upper(proseq = seq->seq, seq->length);

	for (seg=segs; seg!=NULL; seg=seg->next) {
		begin = seg->begin;
		end = seg->end;
		lower(proseq + begin, end - begin + 1);
	}

	fprintf(stdout, "%s\n", seq->header);

	space(leftspace+1);
	for (i=0, ctr=1; i<charline; ++i, ++ctr) {
		if (ctr==10) {
			putc('.', stdout);
			ctr = 0;
		}
		else
			putc(' ', stdout);
	}
	putc('\n', stdout);
	fprintf(stdout, format, 1);

	proseqmax = proseq + seq->length;
	for (i=0, ctr=0; proseq < proseqmax; ++i, ++ctr, ++proseq) {
		if (ctr==charline) {
			putc('\n', stdout);
			fprintf(stdout, format, i+1);
			ctr = 0;
		}
		putc(*proseq, stdout);
	}

	fprintf(stdout, " %d\n", seq->length);
	if (putc('\n', stdout) == EOF) {
		fprintf(stderr, "premature EOF on write\n");
		exit(2);
	}
}

/*---------------------------------------------------------(pretreereport)---*/

pretreereport(seq, segs)
  struct Sequence *seq;
  struct Segment *segs;

{
	struct Sequence	*win;
	struct Segment	*seg;
	char	buffer[100], leftfmt[20], rightfmt[20];
	char	*curseq;
	int	i, left, right, len;
	int	current, nextloent;
	int	cline;

	cline = charline / 2;
   
	fprintf(stdout, "%s\n\n", seq->header);
	sprintf(leftfmt, "%%%ds", cline);
	sprintf(rightfmt, "%%-%ds", cline);

	current = 0;

	for (seg=segs; ; seg=seg->next) {
		if (seg==NULL)
			nextloent = seq->length;
		else
			nextloent = seg->begin;

		if (current < nextloent) {
			left = current;
			right = nextloent - 1;
			len = right - left + 1;
			win = openwin(seq, left, len);
			upper(curseq = win->seq, win->length);

			space(cline);
			fprintf(stdout, " %4d-%-4d ", left+1, right+1);

			while (len>0) {
				if (len<=cline) {
					fwrite(curseq, len, 1, stdout);
					putc('\n', stdout);
					break;
				}
				else {
					fwrite(curseq, cline, 1, stdout);
					putc('\n', stdout);
					space(cline+11);
					curseq += cline;
					len -= cline;
				}
			}
			closewin(win);
		}

		if (seg==NULL) break;

      left = seg->begin;
      right = seg->end;
      len = right - left + 1;
      win = openwin(seq, left, len);
      lower(curseq = win->seq, win->length);

		i = MIN(cline, len);
		if (i < cline)
			fprintf(stdout, "%*s", cline-i, "");
		fwrite(curseq, i, 1, stdout);
		fprintf(stdout, " %4d-%-4d ", left+1, right+1);
		putc('\n', stdout);

		len -= cline;
		if (len>0)
			curseq += cline;

		while (len>0) {
			i = MIN(cline, len);
			if (i < cline)
				space(cline - i);
			fwrite(curseq, i, 1, stdout);
			putc('\n', stdout);
			len -= i;
			if (len>0)
				curseq += i;
		}

		closewin(win);
		current = right + 1;
	}

	putc('\n', stdout);
}

/*-----------------------------------------------------------------(space)---*/

space(len)
	register int len;
{
	static char	spaces[] =
		{' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
		' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
		' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
		' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
		' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '
		};
	register int i;

	while (len > 0) {
		i = MIN(len, sizeof(spaces)/sizeof(spaces[0]));
		fwrite(spaces, i, 1, stdout);
		len -= i;
	}
}

/*----------------------------------------------------------------(seqout)---*/

seqout(seq, hilo, begin, end)
  struct Sequence *seq;
  int hilo;
  int begin, end;

#define HDRLEN 60

{
	char	*proseq, *proseqmax, *id, *header;
   char outbuf[HDRLEN+1];
   static int hi = 1;
   static int lo = 0;
   int i, ctr, iend;

   if (hionly && hilo==lo) return;
   if (loonly && hilo==hi) return;

   proseq = seq->seq;
   proseqmax = proseq + seq->length;
   id = seq->id;
   if (id==NULL) id = seq->parent->id;
   header = seq->header;
   if (header==NULL) header = seq->parent->header;

   iend = findchar(header, ' ');
   if (iend!=-1) header = header+iend;

   if (entinfo)
     {
      fprintf(stdout, ">%s(%d-%d)", id, begin, end);
/*    if (iend!=-1 && strlen(header)<=HDRLEN) fprintf(stdout, "%s", header);
      else if (iend!=-1) for (i=0; i<HDRLEN; i++) putc(header[i], stdout); */
      fprintf(stdout, " complexity=%4.2f (%d/%4.2f/%4.2f)\n",
           seq->entropy, window, locut, hicut);
     }
   else
     {
      fprintf(stdout, ">%s(%d-%d)", id, begin, end);
      if (iend!=-1)   /* fprintf(stdout, "%s\n", header); */
        {
		 i = MIN(HDRLEN, strlen(header));
		 fwrite(header, i, 1, stdout);
         putc('\n', stdout);
        }
      else putc('\n', stdout);
     }
   
   if (hilo==lo)
     {
      lower(proseq, seq->length);
     }
   else if (hilo==hi && seq->length>=hilenmin)
     {
      upper(proseq, seq->length);
     }
   else
     {
      lower(proseq, seq->length);
     }

   for (; proseq < proseqmax; proseq+=i) {
		i = MIN(charline, proseqmax - proseq);
		fwrite(proseq, i, 1, stdout);
		putc('\n', stdout);
	}

	if (putc('\n', stdout) == EOF) {
		fprintf(stderr, "premature EOF on write\n");
		exit(2);
	}
}

/*-------------------------------------------------------------(appendseg)---*/

appendseg(segs, seg)
  struct Segment *segs, *seg;

  {struct Segment *temp;

   temp = segs;
   while (1)
     {
      if (temp->next==NULL)
        {
         temp->next = seg;
         break;
        }
      else
        {
         temp = temp->next;
        }
     }

   return;
  }

/*--------------------------------------------------------------(freesegs)---*/

freesegs(segs)
  struct Segment *segs;

  {struct Segment *temp;

   while (segs!=NULL)
     {
      temp = segs->next;
      free(segs);
      segs = temp;
     }
  }

/*-----------------------------------------------------------------(usage)---*/

usage()

  {
   fprintf(stderr, "\
Usage: segA <file> <window> <locut> <hicut> <Σ> <options>\n\
         <file>   - filename containing fasta-formatted sequence(s) \n\
         <window> - OPTIONAL window size (default 12) \n\
         <locut>  - OPTIONAL low (trigger) complexity (default 2.2) \n\
         <hicut>  - OPTIONAL high (extension) complexity (default locut + 0.3) \n\
         <Σ>      - OPTIONAL alphabet name; one of [prot]|dna|rna \n\
	 <options> \n\
            -x  each input sequence is represented by a single output \n\
                sequence with low-complexity regions replaced by \n\
                strings of 'x' characters \n\
            -c <chars> number of sequence characters/line (default 60)\n\
            -m <size> minimum length for a high-complexity segment \n\
                (default 0).  Shorter segments are merged with adjacent \n\
                low-complexity segments \n\
            -l  show only low-complexity segments (fasta format) \n\
            -h  show only high-complexity segments (fasta format) \n\
            -a  show all segments (fasta format) \n\
            -n  do not add complexity information to the header line \n\
            -o  show overlapping low-complexity segments (default merge) \n\
            -t <maxtrim> maximum trimming of raw segment (default 100) \n\
            -p  prettyprint each segmented sequence (tree format) \n\
            -q  prettyprint each segmented sequence (block format) \n");
   exit(1);
  }

//Initializes a 2d matrix of natural log of factorial values, allowing for certain
//  non-integer inputs
//  Rows indexes are the integer part of a number
//  Column indexes are the modulus of the number wrt the given weight
void lnfacwinit(int weight){
    int q, m;
    double w = (double) weight;
    for(q = 0; q < LNFACMAX; q++){
        _lnfacw[q] = (double*) malloc(sizeof(double) * weight);
        for(m = 0; m < weight; m++){
            if(m == 0){
                _lnfacw[q][m] = lnfac[q];
            } else {
                _lnfacw[q][m] = lgamma((double) q + (double) m / w + 1.);
            }
        }
    }
}

void lnfacwfree(){
    int i;
    for(i = 0; i < LNFACMAX; i++){
        free(_lnfacw[i]);
    }
}

//Uses the global struct Alphabet *abet to get the weight
//Given an integer value, looks up the log gamma value of ((double)x / (double)weight + 1)
//Will segfault if given an x value > LNFACMAX * (w+1) - 1;
double lnfacw(int x){
    int q, m;
    q = x / abet->weight;
    m = x % abet->weight;
    return _lnfacw[q][m];
}
