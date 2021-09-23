/******************************************************************* 
*  Preprocessor for fortran programs                               *
*                                                                  *
*  Here is how it works:                                           *
*                                                                  *
*          preproc -Ddef1=val1 -Ddef2=val2 ... code.f tmp.f        *
*                                                                  *
*  The output is in tmp.f. To compile:                             *
*                                                                  *
*          $(F_CMPLR) $(FFLAGS) tmp.f                              *
*                                                                  *
*  Written in Oct. 1996 by Yang Wang at PSC.                       *
*                                                                  *
*  VERSION 2.0                                                     *
*          This new release also takes into account of #ifdef and  *
*          #ifndef                                                 *
*                                                                  *
*          preproc -Darg1 -Darg2 ... code.f                        *
********************************************************************/
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define LCOL 80    /* last col used                               */
#define MAXC 80    /* longest line allowed                        */
#define MAXD 100   /* The maximum number of definitions           */
 
struct str_array {
       char *def;
       char *val;
       char *def_eq_val;
       int  ind;
};

struct str_array cl_define[MAXD];
static int num_cl_define;
static int line_num;

char   lcase(char);
int    compare(char *, char *);
void   replace(char *, char *, char *);
char  *ncs_strstr(char *, char *);
void   string_to_lcase(char *);
 
int main(int argc, char *argv[])
{
     FILE *ifp;
     void procfile(FILE *);
     char *p, *t;
     int  i,j;
     int equ;
 
     if (argc < 2) {
        fprintf(stderr,"In correct number of arguments\n");
        exit(1);
     }
 
     num_cl_define = 0;
     argv++;

     for (i=0; i < argc-2; i++) {
        equ=0;
        if (compare("-D",*argv)) {
           p=*argv+2;
           t=p;
           while (*t != '\0') {
              if(*++t == '=') {
                 equ=1;
                 break;
              }
           }
           j=t-p;
        }
        else {
           fprintf(stderr, "The argument should be: -Dparam[=value]");
           exit(1);
        }
        if (j>0) {
           cl_define[num_cl_define].def=(char *)malloc((j+1)*sizeof(char));
           strncpy(cl_define[num_cl_define].def,p,j);
           *(cl_define[num_cl_define].def+j)='\0';
           if (equ==1) {
              p=t+1;
              cl_define[num_cl_define].def_eq_val=
                 (char *)malloc((j+strlen(p)+3)*sizeof(char));
              strcpy(cl_define[num_cl_define].def_eq_val,
                     cl_define[num_cl_define].def);
              strcat(cl_define[num_cl_define].def_eq_val,"==");
              j=strlen(p);
              if (j < 1) {
                 fprintf(stderr, "The argument should be: -Dparam[=value]");
                 exit(1);   
              }
              cl_define[num_cl_define].val=(char *)malloc((j+1)*sizeof(char));
              strcpy(cl_define[num_cl_define].val,p);
              strcat(cl_define[num_cl_define].def_eq_val,p);
           }
           else {
              cl_define[num_cl_define].def_eq_val=(char *)NULL;
              cl_define[num_cl_define].val=(char *)NULL;
           }
           cl_define[num_cl_define].ind=num_cl_define;
           num_cl_define++;
        }
        argv++;
     }
 
     if ((ifp = fopen(*argv, "r")) == NULL) {
         fprintf(stderr, "format: can't open %s\n", *argv);
         exit(1);   
     }

     procfile(ifp);
     fclose(ifp);    
     return 0;
}
 

/* read and process the input file 
   write to the output file                                     */
void procfile (FILE *ifp)
{
     char line[MAXC+1];
     char *p, *t;
     char cif[]          = "#if";
     char cifdef[]       = "#ifdef";
     char cifndef[]      = "#ifndef";
     char celif[]        = "#elif";
     char celse[]        = "#else";
     char cendif[]       = "#endif";
     char cdefine[]      = "#define";
     struct str_array sc_define[MAXD];
     int  i;
     int  num_def= -1;
     int  num_delete;
     int  if_level;
     int  ifeq_on;

     num_delete = 0;
     line_num=0;
     if_level=0;
     ifeq_on= 0;
     while (fgets(line,MAXC+1,ifp) != NULL) {
        line_num++;
        if (line[0] == '#') {
           if (compare(cif,line))
              if_level++;

           if (compare(cif,line) || compare(celif,line)) {
              p=line;
              while (p < line+LCOL && *p != ' ')
                 p++;

              num_delete = 1;

              if (compare(cifndef,line)) {
                 num_delete = 0;
                 ifeq_on=1;
                 for (i=0; i<num_cl_define; i++) {
                    if (!cl_define[i].def_eq_val) {
                       if (compare(cl_define[i].def,p)) {
                          num_delete = 1;
                          ifeq_on=0;
                          break;
                       }
                    }
                 }
              }
              else {
                 for (i=0; i<num_cl_define; i++) {
                    if (cl_define[i].def_eq_val) {
                       if (compare(cl_define[i].def_eq_val,p) ||
                           compare(cl_define[i].def,p)) {
                          num_delete = 0;
                          ifeq_on=1;
                          break;
                       }
                    }
                    else {
                       if (compare(cl_define[i].def,p)) {
                          num_delete = 0;
                          ifeq_on=1;
                          break;
                       }
                    }
                 }
              }
           }
           else if(compare(celse,line)) {
              num_delete = 1-num_delete+ifeq_on;
           }
           else if(compare(cendif,line)) {
              num_delete = 0;
              ifeq_on= 0;
              if_level--;
           }
           else if(compare(cdefine,line) && num_delete == 0) {
              p=line+1;
              t=strtok(p," \t");
              p=(char *)NULL;
              if ((t = strtok(p," \t")) != (char *)NULL) {
                 num_def++;
                 if (num_def > MAXD) exit(1);
/*               printf("token = %s\n",t);                                   */
                 sc_define[num_def].def=
                    (char *)malloc((strlen(t)+1)*sizeof(char));
                 sc_define[num_def].ind=num_def;
                 string_to_lcase(t);
                 strcpy(sc_define[num_def].def,t);
                 p=(char *)NULL;
                 if ((t = strtok(p," \t")) != (char *)NULL){
/*                  printf("token = %s\n",t)                                 */
                    sc_define[num_def].val=
                       (char *)malloc((strlen(t)+1)*sizeof(char));
                    strcpy(sc_define[num_def].val,t);
                    sc_define[num_def].ind=num_def;
                 }
                 else {
/*                  printf("token = %s\n",t);                                */
                    sc_define[num_def].val=(char *)malloc(sizeof(char)+1);
                    strcpy(sc_define[num_def].val," ");
                 }
              }
           }
           else if(num_delete == 0) {
              printf("PREPROC. SYNTAX ERROR at line no.: %d\n",line_num);
              exit(1);
           }
        }
 
        if (num_delete == 0 && line[0] != '#') {
           /* truncate beyond column 72 and strip trailing blanks  */
           if (line[0] != 'c') {
              p=line;     
              while (p < line+LCOL && *p != '\0' && *p != '\n') p++;
              p--;
              while (*p == ' ') p--;
              *++p = '\0';
              p=line;
              for(i=0;i<=num_def;i++)
                 replace(sc_define[i].def,sc_define[i].val,p);
           }
           puts(line);
        }
     }
     if (!feof(ifp))
        fprintf(stderr,"\nError reading file.\n");
}
 

/* compare string s to t, not case sensitive
   string s is usually shorter than t                           */
int compare (char *s, char *t)  
{
    while (*s != '\0') {
	if (*t == '\0') 
           return(0);
        while (*s == ' ')
           s++;
        while (*t == ' ')
           t++;
        if (lcase(*s++) != lcase(*t++)) 
           return (0);
    }
    return(1);
}
 

/* replace string d with v in s, not case sensitive */
void replace (char *d, char *v, char *s)  
{
    char  *p;
    char  *t;
    char  newl[MAXC+1];
    short fl;
    int   i,nd,nv;

    nd=strlen(d);
    nv=strlen(v);
    p=s;
    while ((p = ncs_strstr(p,d)) != (char *)NULL) {
       fl=1;
       if ((t = p-1) != (char *)NULL) 
          if(('a' <= *t && *t <= 'z') || *t == '_' || ('0' <= *t && *t <= '9'))
             fl=0;
       if ((t = p+nd) != (char *)NULL) 
          if(('a' <= *t && *t <= 'z') || *t == '_' || ('0' <= *t && *t <= '9'))
             fl=0;
       if (fl == 1) {
          t=p+nd;
          strcpy(newl,t);
          t=v;
          for (i=0; i < nv; i++) {
             *p++ = *t++;
          }
          p--;
          i=p-s;
          if (strlen(newl)+i <= LCOL)
             strcpy(p,newl);
          else {
             fprintf(stderr,"Line %d: too long after the subst.\n",line_num);
             exit(1);
          }
       }  
    }
}


/* non-case sensitive strstr */
char *ncs_strstr(char *s1, char *s2)
{
      int  i;
      char *s0;
      char *p;

      s0=(char *)malloc((strlen(s1)+1)*sizeof(char));
      strcpy(s0,s1);
      string_to_lcase(s0);
      p=strstr(s0,s2);
      if(p != (char *)NULL) {
         i=p-s0;
         free(s0);
         return(s1+i);
      }
      else {
         free(s0);
         return((char *)NULL);
      }
}


/* convert the string to low case */
void string_to_lcase(char *s)
{
     char *p;

     p=s;
     while (*p != '\0') {
        *p = lcase(*p);
        p++;
     }
}


/* convert to lower case                            */
char lcase(char c)         
{
    if (c < 'A') return(c);
    else if (c > 'Z') return(c);
    else return (c+'a'-'A');
}

