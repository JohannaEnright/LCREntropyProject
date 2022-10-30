
/*****************************************************************************/
/***   (genwin.c)                                                          ***/
/*****************************************************************************/

/*--------------------------------------------------------------(includes)---*/

#include "genwin.h"
#include "alphabet.h"

/*---------------------------------------------------------------(defines)---*/

#define STRSIZE 100

/*----------------------------------------------------------------(protos)---*/

struct Sequence *readentry();

/*---------------------------------------------------------------(globals)---*/

char *blastdbs[] =
  {"bba", "bbn", "embl", "gbupdate", "genbank", "genpept", "gpupdate",
   "nr", "nrdb", "nrdb.shuf", "pir", "pseq", "swissprot", "tfdaa"};

int nblastdbs = 14;

#ifndef MIN
#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#endif

char *blastdir = "/net/cruncher/usr/ncbi/db/fasta/";
char *indexdir = "/net/cruncher/usr/ncbi/db/index/";

//int nabets;
//struct Alphabet **abets;
//int ntvecs;
//struct TransVector **tvecs;
//int nsvecs;
//struct ScoreVector **svecs;
//int nsmats;
//struct ScoreMatrix **smats;

struct Alphabet *abet;

struct strlist
  {
   char string[STRSIZE];
   struct strlist *next;
  } *str, *curstr;

/*---------------------------------------------------------------(tmalloc)---*/

#define TESTMAX 1000
void *tmalloc();
int record_ptrs[TESTMAX] = {0,0,0,0};
int rptr = 0;

/*------------------------------------------------------------(genwininit)---*/

genwininit(char* name){
        abet = openabet(name);
	return;
}

/*------------------------------------------------------------(genwinterm)---*/

genwinterm(){
    closeabet(abet);
}

/*-------------------------------------------------------------(loadabet)----*/

extern struct Alphabet *openabet(char* name){
    struct Alphabet *abet;
    int i,j;
    char c,s;
    char buf[100];
    int bufsize = sizeof(buf)/sizeof(buf[0]);

    abet = (struct Alphabet *) malloc(sizeof(struct Alphabet));

    //Initialize structure members
    abet->name = name;
    abet->size = 0;
    abet->weight = 0;
    //Initialize size and indicies for each symbol in the alphabet
    for (j = 0; j < sizeof(abet->symsize)/sizeof(abet->symsize[0]); ++j) {
    	abet->symsize[j] = 0;
        abet->symweight[j] = 0;
        abet->symindicies[j] = NULL;
    }
    //Initialize the indexes for each residue in the alphabet
    for (j = 0; j < sizeof(abet->index)/sizeof(abet->index[0]); ++j) {
        abet->index[j] = -1;
    }

    //Search through the packaged alphabet file for the requested alphabet name
    i=0;
    while(i < alphabet_txt_len && *name != '\0'){
        do{
            c = alphabet_txt[i++];
        } while(i < alphabet_txt_len && c != '>');
        while(i < alphabet_txt_len && alphabet_txt[i] != '\n' && alphabet_txt[i++] == *name){
            name++;
        }
        if(*name != '\0'){
            name = abet->name;
        }
    }
    if(*name != '\0'){
        fprintf(stderr,"[ERROR] Could not recognize the alphabet '%s'\n",abet->name);
        exit (EXIT_FAILURE);
    }
    i++;
    //Process the Alphabet
    while(i < alphabet_txt_len && alphabet_txt[i] != '>'){
        //First character on each line is the symbol for a potentially ambiguous character
        s = toupper(alphabet_txt[i++]);
        //Each subsequent character is a residue in the alphabet
        while(abet->symsize[s] < bufsize && i < alphabet_txt_len &&
                (c = alphabet_txt[i++]) != '\n'){
            buf[abet->symsize[s]++] = toupper(c);
            if(toupper(c) != tolower(c)){
                abet->symsize[tolower(s)]++;
            }
        }
        if(abet->symsize[s] >= bufsize){
            fprintf(stderr,"[WARNING] Encountered an alphabet symbol representing %d or more residues: Might indicate alphabet corruption\n",bufsize);
            abet->symsize[s] = bufsize - 1;
            abet->symsize[tolower(s)] = bufsize - 1;
        }
        abet->symindicies[s] = (char*) malloc(abet->symsize[s] * sizeof(char));
        memcpy(abet->symindicies[s],buf,abet->symsize[s] * sizeof(char));
        abet->symindicies[tolower(s)] = abet->symindicies[s];
        for(j = 0; j < abet->symsize[s]; j++){
            if(abet->index[buf[j]] == -1){
                abet->index[buf[j]] = abet->size;
                abet->index[tolower(buf[j])] = abet->size++;
            }
        }
        abet->weight = lcm(abet->weight,abet->symsize[s]);
    }
    //Now that the LCM of the symsizes has been found, the weight for each symbol can be
    //calculated
    for (j = 0; j < sizeof(abet->symsize)/sizeof(abet->symsize[0]); ++j) {
        if(abet->symsize[j] > 0){
            abet->symweight[j] = abet->weight / abet->symsize[j];
        }
    }
    //Now that the size of the alphabet is known the index of non-alphabet characters can be
    //set
    for (j = 0; j < sizeof(abet->index)/sizeof(abet->index[0]); ++j) {
        if(abet->index[j] == -1){
            abet->index[j] = abet->size;
        }
    }

    /*
        fprintf(stderr,"Name: %s\nSize: %d\nWeight: %d\n",abet->name,abet->size,abet->weight);
        for(i = 0; i < 128; i++){
            if(abet->index[i] < abet->size){ fprintf(stderr,"%c",i);}
        }
        fprintf(stderr,"\n");
        for(i = 0; i < 128; i++){
            if(abet->symsize[i] > 0){ fprintf(stderr,"%c:%d ",i,abet->symsize[i]);}
        }
        fprintf(stderr,"\n");
     /**/

    return abet;
}


/*-------------------------------------------------------------(closeabet)---*/

extern closeabet(struct Alphabet *abet){
    int i;

    for (i = 0; i < sizeof(abet->symsize)/sizeof(abet->symsize[0]); ++i) {
        if(abet->symindicies[i] != NULL) {
            free(abet->symindicies[i]);
            abet->symindicies[tolower(i)] = NULL;
        }
    }
    free(abet);
}
        
/*-------------------------------------------------------------(opendbase)---*/

extern struct Database *opendbase(name)
  char *name;

  {struct Database *dbase;

   dbase = (struct Database *) malloc(sizeof(struct Database));

   if (blastdb(name))
     {
      dbase->filename = (char *) malloc(strlen(blastdir)+strlen(name)+1);
      dbase->indexname = (char *) malloc(strlen(indexdir)+strlen(name)+1);
      strcpy(dbase->filename, blastdir);
      strcat(dbase->filename, name);
      strcpy(dbase->indexname, indexdir);
      strcat(dbase->indexname, name);
     }
   else
     {
      dbase->filename = (char *) malloc(strlen(name)+1);
      dbase->indexname = (char *) malloc(strlen(name)+1);
      strcpy(dbase->filename, name);
      strcpy(dbase->indexname, name);
     }

   if (strcmp(dbase->filename, "-")==0)
     {
      dbase->fp = stdin;
     }
   else if ((dbase->fp=fopen(dbase->filename, "r"))==NULL)
     {
      free(dbase->filename);
      free(dbase->indexname);
      free(dbase);
      return((struct Database *) NULL);
     }

   dbase->filepos = 0L;

   return(dbase);
  }

/*---------------------------------------------------------------(blastdb)---*/

int blastdb(name)
  char *name;

  {int i;

   for (i=0; i<nblastdbs; i++)
     {
      if (strcmp(name, blastdbs[i])==0) {return(TRUE);}
     }

   return(FALSE);
  }

/*------------------------------------------------------------(closedbase)---*/

extern closedbase(dbase)
  struct Database *dbase;

  {
   fclose(dbase->fp);
   free(dbase->filename);
   free(dbase->indexname);
   free(dbase);

   return;
  }

/*--------------------------------------------------------------(firstseq)---*/

extern struct Sequence *firstseq(dbase)
  struct Database *dbase;

  {
   if (dbase->filepos!=0L)
     {
      dbase->filepos = 0L;
      if (fseek(dbase->fp, dbase->filepos, 0)!=0)
        {fprintf(stderr, "Error positioning file %s for firstseq.\n",
                           dbase->filename);
         exit(1);}
     }

   return(readentry(dbase));
  }

/*---------------------------------------------------------------(nextseq)---*/

extern struct Sequence *nextseq(dbase)
  struct Database *dbase;

  {
   return(readentry(dbase));
  }

/*--------------------------------------------------------------(closeseq)---*/

extern closeseq(seq)
  struct Sequence *seq;

  {
   if (seq==NULL) return;

   if (seq->id!=NULL)          free(seq->id);
   if (seq->name!=NULL)        free(seq->name);
   if (seq->organism!=NULL)    free(seq->organism);
   if (seq->header!=NULL)      free(seq->header);
   if (seq->state!=NULL)       free(seq->state);
   if (seq->composition!=NULL) free(seq->composition);
   free(seq->seq);
   free(seq);
   return;
  }

/*---------------------------------------------------------------(openwin)---*/

extern struct Sequence *openwin(parent, start, length)
  struct Sequence *parent;
  int start, length;

  {struct Sequence *win;
   int i;

   //fprintf(stderr, "openwin: %d+%d\n",start,length);

   if (start<0 || length<0 || start+length>parent->length)
     {
      return((struct Sequence *) NULL);
     }

   win = (struct Sequence *) malloc(sizeof(struct Sequence));

/*---                                          ---[set links, up and down]---*/

   win->parent = parent;
   if (parent->root==NULL)
     {win->root = parent;}
   else
     {win->root = parent->root;}
   win->children = (struct Sequence **) NULL;

/* parent->children = ***foo***                   ---[not yet implemented]---*/

   win->id = (char *) NULL;
   if(parent->id){
       win->id = (char*) malloc (sizeof(char) * (strlen(parent->id) + 1));
       memcpy(win->id,parent->id,sizeof(char) * strlen(parent->id) + 1);
   }
   win->name = (char *) NULL;
   win->organism = (char *) NULL;
   win->header = (char *) NULL;

/*---                          ---[install the local copy of the sequence]---*/

   win->start = start;
   win->length = length;
#if 0
   win->seq = (char *) malloc(sizeof(char)*length + 1);
   memcpy(win->seq, (parent->seq)+start, length);
   win->seq[length] = '\0';
#else
	win->seq = parent->seq + start;
#endif

/*---                          ---[setup window implementation parameters]---*/

/*---                                                 ---[set local flags]---*/

	win->rubberwin = FALSE;
	win->floatwin = FALSE;
	win->punctuation = FALSE;

/*---                                   ---[initially unconfiguerd window]---*/

	win->entropy = -2.;
	win->state = (int *) NULL;
	win->composition = (int *) NULL;
	win->classvec = (char *) NULL;
	win->scorevec = (double *) NULL;
        win->alphabet = parent->alphabet;

	stateon(win);

	return win;
}

/*---------------------------------------------------------------(nextwin)---*/

extern struct Sequence *nextwin(win, shift)
  struct Sequence *win;
  int shift;

  {
   if ((win->start+shift)<0 ||
       (win->start+win->length+shift)>win->parent->length)
     {
      return((struct Sequence *) NULL);
     }
   else
     {
      return(openwin(win->parent, win->start+shift, win->length));
     }
  }

/*--------------------------------------------------------------(shiftwin1)---*/
static void	decrementsv(), incrementsv();

extern int shiftwin1(win)
	struct Sequence	*win;
{
	register int	i, j, length, idx;
	register int	*comp;


	length = win->length;
	comp = win->composition;

        int x;

        //fprintf(stderr,"comp: ");
        //for(x = 0; x < abet->size; x++){
        //    fprintf(stderr,"%d ",comp[x]);
        //}
        //fprintf(stderr,"\n");

	if ((++win->start + length) > win->parent->length) {
		--win->start;
		return FALSE;
	}

	if (abet->symsize[j = win->seq[0]] > 0){
            //fprintf(stderr, "left: %c\n",j);
            for(i = 0; i < abet->symsize[j]; i++){
                idx = abet->index[abet->symindicies[j][i]];
                decrementsv(win->state, comp[idx],abet->symweight[j]);
                comp[idx] -= abet->symweight[j];
            }
        }

	j = win->seq[length];
	++win->seq;


	if (abet->symsize[j] > 0)
            //fprintf(stderr, "right: %c\n",j);
            for(i = 0; i < abet->symsize[j]; i++){
                idx = abet->index[abet->symindicies[j][i]];
                incrementsv(win->state, comp[idx], abet->symweight[j]);
                comp[idx] += abet->symweight[j];
            }

	if (win->entropy > -2.)
		win->entropy = entropy(win->state);

	return TRUE;
}

/*--------------------------------------------------------------(closewin)---*/

extern closewin(win)
  struct Sequence *win;

  {
   if (win==NULL) return;

   if (win->state!=NULL)       free(win->state);
   if (win->composition!=NULL) free(win->composition);
   if (win->classvec!=NULL)    free(win->classvec);
   if (win->scorevec!=NULL)    free(win->scorevec);

   free(win);
   return;
  }

/*----------------------------------------------------------------(compon)---*/

extern compon(win)
	struct Sequence	*win;
{
	register int	*comp;
	register int	sym,res;
	register char	*seq, *seqmax;

	win->composition = comp = (int *) calloc(abet->size*sizeof(*comp), 1);
	seq = win->seq;
	seqmax = seq + win->length;

	while (seq < seqmax) {
		sym = *seq++;
		for (res=0; res < abet->symsize[sym];res++){
                    comp[abet->index[abet->symindicies[sym][res]]] += abet->symweight[sym];
                }
	}

	return;
}

/*---------------------------------------------------------------(stateon)---*/

static int state_cmp(s1, s2)
	int	*s1, *s2;
{
	return *s2 - *s1;
}

extern stateon(win)
	struct Sequence	*win;
{
	register int	res, nel, c;

	if (win->composition == NULL)
		compon(win);

	win->state = (int *) malloc((abet->size+1)*sizeof(win->state[0]));

	for (res = nel = 0; res < abet->size; ++res) {
		if ((c = win->composition[res]) == 0)
			continue;
		win->state[nel++] = c;
	}
	for (res = nel; res <= abet->size; ++res)
		win->state[res] = 0;

	qsort(win->state, nel, sizeof(win->state[0]), state_cmp);

	return;
}

/*-----------------------------------------------------------------(enton)---*/

extern enton(win)
  struct Sequence *win;

  {
   if (win->state==NULL) {stateon(win);}

   win->entropy = entropy(win->state);

   return;
  }

/*---------------------------------------------------------------(entropy)---*/
static int		thewindow;
static double	*entray;

#define LN2	0.69314718055994530941723212145818

void
entropy_init(window)
	int	window;
{
	int		i;
	double	x, xw;


	thewindow = window * abet->weight;

	entray = (double *)malloc((thewindow+1) * sizeof(*entray));
        xw = (double)thewindow;
	for (i = 1; i <= thewindow; ++i) {
		x = i / xw;
		entray[i] = -x * log(x) / LN2;
                //if(isnan(entray[i])){
                //    fprintf(stderr,"H: %d, %0.2f/%d[%d*%d] -> (%0.2f, %0.2f) = %0.2f\n",i,xw,thewindow,window,abet->weight,x,log(x),entray[i]);
                //}
	}

}

extern double entropy(sv)
	register int	*sv;
{
	int	*sv0 = sv;
	register double	ent;
	register int	i, total;
	register int	*svmax;
	register double	xtotrecip, xsv;

	for (total = 0; (i = *sv) != 0; ++sv)
		total += i;
	svmax = sv;
	ent = 0.0;
	if (total == thewindow) {
		for (sv = sv0; sv < svmax; sv++) {
			ent += entray[*sv];
		}
                //if(ent == 0){
                //    fprintf(stderr,"h: %0.2f sv: ",ent);
                //    for(i = 0; i < abet->size; i++){
                //        fprintf(stderr,"%d ",sv0[i]);
                //    }
                //    fprintf(stderr,"\n");
                //    exit(EXIT_FAILURE);
                //}
                //fprintf(stderr,"window\n");
		return ent;
	}
	if (total == 0){
                //fprintf(stderr,"empty\n");
		return 0.;
        }

	xtotrecip = 1./(double)total;
	for (sv = sv0; sv < svmax; ) {
		xsv = *sv++;
		ent += xsv * log(xsv * xtotrecip);
	}

        //fprintf(stderr,"sequence\n");
	return -ent * xtotrecip / LN2;
}

/*-----------------------------------------------------------(decrementsv)---*/

static void
decrementsv(sv, class, weight)
	register int	*sv;
	register int	class;
        register int    weight;
{
	register int	svi;

        //int i;
        //fprintf(stderr,"decrementsv: %d, %d\n",class, weight);
        //for(i = 0; i < abet->size; i++){
        //    fprintf(stderr,"%d ",sv[i]);
        //}
        //fprintf(stderr,"\n");

	while ((svi = *sv++) != 0) {
	    if (svi == class && *sv < class) {
                svi -= weight;
                if(svi < 0){ //Negative counts indicate an error
                    fprintf(stderr,"[ERROR] Negative composition weight encountered: Contact the developer\n");
                    exit(EXIT_FAILURE);
                }
	        sv[-1] = svi;
                //Ensure the state vector remains in descending order
                //the state vector always has a zero at the end and will terminate
                while(svi < *sv){ 
                    sv[-1] = *sv;
                    *sv++ = svi;
                }
	        break;
	    }
	}
}

/*-----------------------------------------------------------(incrementsv)---*/

static void
incrementsv(sv, class, weight)
	register int	*sv;
	register int	class;
        register int    weight;
{
    register int svi; 
    register int *sv0 = sv;


        //int i;
        //fprintf(stderr,"incrementsv: %d, %d\n",class, weight);
        //for(i = 0; i < abet->size; i++){
        //    fprintf(stderr,"%d ",sv[i]);
        //}
        //fprintf(stderr,"\n");

	for (;;) {
	    if (*sv++ == class) {
                svi = sv[-1] + weight;
	        sv[-1] = svi;
                //Maintain sv order
                while(--sv > sv0 && svi > sv[-1]){
                    *sv = sv[-1];
                    sv[-1] = svi;
                }
	        break;
	    }
	}
}

/*-------------------------------------------------------------(readentry)---*/

struct Sequence *readentry(dbase)
  struct Database *dbase;

  {struct Sequence *seq;
   int	c;

   seq = (struct Sequence *) malloc(sizeof(struct Sequence));

   seq->db = dbase;

/*---                                    ---[backpointers null at the top]---*/

   seq->parent = (struct Sequence *) NULL;
   seq->root = (struct Sequence *) NULL;
   seq->children = (struct Sequence **) NULL;

/*---                                                       ---[set flags]---*/

   seq->rubberwin = FALSE;
   seq->floatwin = FALSE;

/*---                                                  ---[read from file]---*/

   if (!readhdr(seq))
     {
      return((struct Sequence *) NULL);
     }
   while (1)  /*---[skip multiple headers]---*/
     {
      c = getc(dbase->fp);
	  if (c == EOF)
		break;
      if (c != '>') {
         ungetc(c, dbase->fp);
         break;
		}
      while ((c=getc(dbase->fp)) != EOF && c !='\n')
		;
		if (c == EOF)
			break;
     }
   readseq(seq);

/*---                                   ---[set implementation parameters]---*/

/*---                                          ---[initially unconfigured]---*/

   seq->entropy = -2.;
   seq->state = (int *) NULL;
   seq->composition = (int *) NULL;
   seq->classvec = (char *) NULL;
   seq->scorevec = (double *) NULL;
   seq->alphabet = abet;

   return(seq);
  }

/*---------------------------------------------------------------(readhdr)---*/

readhdr(seq)
  struct Sequence *seq;

  {FILE *fp;
   char *bptr, *curpos;
   int	c, i, itotal;
   int idend, namend, orgend;

   fp = seq->db->fp;

   if ((c=getc(fp)) == EOF)
     {
      free(seq);
      return(FALSE);
     }
   
   while (c != EOF && isspace(c))
     {
      c = getc(fp);
     }

   if (c!='>')
     {fprintf(stderr, "Error reading fasta format - '>' not found.\n");
      exit(1);}
   ungetc(c, fp);
/*                                               ---[read the header line]---*/
   str = (struct strlist *) malloc (sizeof(struct strlist));
   str->next = NULL;
   curstr = str;

   for (i=0,itotal=0,c=getc(fp); c != EOF; c=getc(fp))
     {
      if (c=='\n') break;

      if (i==STRSIZE-1)
        {curstr->string[i] = '\0';
         curstr->next = (struct strlist *) malloc (sizeof(struct strlist));
         curstr = curstr->next;
         curstr->next = NULL;
         i = 0;}

      curstr->string[i] = c;
      itotal++;
      i++;
     }

   curstr->string[i] = '\0';
   seq->header = (char *) malloc (itotal+2);
   seq->header[0] = '\0';

   for (curstr=str, curpos=seq->header; curstr!=NULL;)
     {
      if (curstr->next==NULL)
        {memccpy(curpos, curstr->string, '\0', STRSIZE);}
      else
        {memccpy(curpos, curstr->string, '\0', STRSIZE-1);}

      str = curstr;
      curstr = curstr->next;
      free (str);

      if (curstr!=NULL) {curpos = curpos+STRSIZE-1;}
     }

   bptr = (seq->header)+1;
   seq->name = (char *) NULL;
   seq->organism = (char *) NULL;
/*                                                   ---[parse out the id]---*/
   idend = findchar(bptr, ' ');
   if (idend==-1) {idend = findchar(bptr, '\n');}
   if (idend==-1) {idend = findchar(bptr, '\0');}
   if (idend==-1)
     {fprintf(stderr, "Error parsing header line - id.\n");
      fputs(seq->header, fp);
      exit(1);}   

   seq->id = (char *) malloc((idend+1)*sizeof(char));
   memcpy(seq->id, bptr, idend);
   seq->id[idend] = '\0';

   if (bptr[idend]=='\n' || bptr[idend]=='\0') {return(TRUE);}

/*                                         ---[parse out the protein name]---*/
   bptr = bptr + idend + 1;
   while (bptr[0]==' ') {bptr++;}

   namend = findchar(bptr, '-');
   if (namend==-1) {namend = findchar(bptr, '\n');}
   if (namend==-1) {namend = findchar(bptr, '\0');}
   if (namend==-1)
     {fprintf(stderr, "Error parsing header line - name.\n");
      fputs(seq->header, fp);
      return(TRUE);}

   seq->name = (char *) malloc((namend+1)*sizeof(char));
   memcpy(seq->name, bptr, namend);
   seq->name[namend] = '\0';

   if (bptr[namend]=='\n' || bptr[namend]=='\0') {return(TRUE);}

/*                                                 ---[parse out organism]---*/
   bptr = bptr + namend + 1;
   while (bptr[0]==' ') {bptr++;}

   orgend = findchar(bptr, '|');
   if (orgend==-1) {orgend = findchar(bptr, '#');}
   if (orgend==-1) {orgend = findchar(bptr, '\n');}
   if (orgend==-1) {orgend = findchar(bptr, '\0');}
   if (orgend==-1)
     {fprintf(stderr, "Error parsing header line - organism.\n");
      fputs(seq->header, fp);
      return(TRUE);}

   seq->organism = (char *) malloc((orgend+1)*sizeof(char));
   memcpy(seq->organism, bptr, orgend);
   seq->organism[orgend] = '\0';

/*                                    ---[skip over multiple header lines]---*/
   while (TRUE)
     {
      c = getc(fp);
	  if (c == EOF)
		return(TRUE);
      if (c=='>')
        {
         skipline(fp);
        }
      else
        {
         ungetc(c,fp);
         break;
        }
     }

   return(TRUE);
  }

/*--------------------------------------------------------------(skipline)---*/

skipline(fp)
  FILE *fp;

  {int	c;

   while ((c=getc(fp))!='\n' && c!=EOF)
     ;

   return;
  }

/*--------------------------------------------------------------(findchar)---*/

extern int findchar(str, chr)
  char *str;
  char chr;

  {int i;

   for (i=0; ; i++)
     {
      if (str[i]==chr)
        {
         return(i);
        }
      if (str[i]=='\0')
        {
         return(-1);
        }
     }
   }

/*---------------------------------------------------------------(readseq)---*/

readseq(seq)
  struct Sequence *seq;

{FILE *fp;
   int i, itotal;
   int	c;
   char *curpos;

   fp = seq->db->fp;

   seq->punctuation = FALSE;   

   str = (struct strlist *) malloc (sizeof(struct strlist));
   str->next = NULL;
   curstr = str;

	for (i = 0, itotal = 0, c = getc(fp); c != EOF; c = getc(fp)) {
		if (abet->symsize[c] > 0) {
Keep:
			if (i < STRSIZE-1) {
				curstr->string[i++] = c;
				continue;
			}
			itotal += STRSIZE-1;
			curstr->string[STRSIZE-1] = '\0';
			curstr->next = (struct strlist *) malloc(sizeof(*curstr));
			curstr = curstr->next;
			curstr->next = NULL;
			curstr->string[0] = c;
			i = 1;
			continue;
                }

		switch (c) {
		case '>':
			ungetc(c, fp);
			goto EndLoop;
		case '*': case '-':
			seq->punctuation = TRUE;
			goto Keep;
		default:
			continue;
		}
	}
EndLoop:
	itotal += i;

	curstr->string[i] = '\0';
	seq->seq = (char *) malloc (itotal+2);
	seq->seq[0] = '\0';

	for (curstr = str, curpos = seq->seq; curstr != NULL;) {
		if (curstr->next == NULL)
			memccpy(curpos, curstr->string, '\0', STRSIZE);
		else
        	memccpy(curpos, curstr->string, '\0', STRSIZE-1);

		str = curstr;
		curstr = curstr->next;
		free(str);

		if (curstr != NULL)
			curpos = curpos+STRSIZE-1;
	}

	seq->length = strlen(seq->seq);

	return;
}
/*-----------------------------------------------------------------(upper)---*/

extern upper(string, len)
	register char	*string;
	size_t	len;
{
	register char	*stringmax, c;

	for (stringmax = string + len; string < stringmax; ++string)
		if (islower(c = *string))
			*string = toupper(c);
}

/*-----------------------------------------------------------------(lower)---*/

extern lower(string, len)
	char	*string;
	size_t	len;
{
	register char	*stringmax, c;

	for (stringmax = string + len; string < stringmax; ++string)
		if (isupper(c = *string))
			*string = tolower(c);
}

/*-------------------------------------------------------------------(min)---*/

int min(a, b)
  int a, b;

  {
   if (a<b) {return(a);}
   else {return(b);}
  }

/*-------------------------------------------------------------------(max)---*/

int max(a, b)
  int a, b;

  {
   if (a<b) {return(b);}
   else {return(a);}
  }

/*-------------------------------------------------------------------(gcd)---*/

extern int gcd(int a, int b){
    int i = (a > b) ? 0 : 1;
    int *num[2];
    if(a == b){
        return a;
    }
    num[0] = &a;
    num[1] = &b;
    do{
        *num[i] = *num[i] % *num[!i];
        i = !i;
    } while(*num[!i] != 0 && *num[i] > 1);
    return *num[i];
}

/*-------------------------------------------------------------------(lcd)---*/

extern int lcm(int a, int b){
    if(min(a,b) == 0 || max(a,b) % min(a,b) == 0){
        return max(a,b);
    }
    return (a*b)/gcd(a,b);
}

/*---------------------------------------------------------------------------*/
