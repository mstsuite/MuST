#include <iotk.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#ifndef NoXDR_format
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif

struct FileListStruct {
    FILE *fp;
    char *fname;
    char *fstatus;
    char *fform;
    int index;
#ifndef NoXDR_format
    XDR xdrs;
#endif
    struct FileListStruct *next;
    struct FileListStruct *prev;
};

int NumFiles = 0;
struct FileListStruct *flh = 0;
struct FileListStruct *pfl = 0;

#define XDR_STRING_SIZE(x) (x + (4 - x % 4) % 4 + 4)

#ifdef Integer2long
#define INTEGER long
#define xdr_INTEGER xdr_long
#else
#define INTEGER int
#define xdr_INTEGER xdr_int
#endif

char lcase(char);
void string_to_lcase(char *);
void searchFile(const INTEGER *);

void c_dtsize(INTEGER *intger4_size, INTEGER *real8_size) {
     *intger4_size=sizeof(INTEGER);
     *real8_size=sizeof(double);
}

void c_gopen(INTEGER *findex, char *fname, INTEGER *fname_len, 
              char *fstatus, INTEGER *fstatus_len,
              char *fform, INTEGER *fform_len,
              INTEGER lfname, INTEGER lfstatus) {
    int i, jx;
    int n = *fname_len;
    int m = *fstatus_len;
    int k = *fform_len;

    if (NumFiles == 0) {
        flh = (struct FileListStruct *)malloc(sizeof(struct FileListStruct));
        flh->prev = 0;
        pfl = flh;
        jx = 0;
    }
    else {
        pfl = flh;
        for (i=0; i<NumFiles; i++) {
            pfl = pfl->next;
        }
        jx = pfl->index;
        pfl->next =
            (struct FileListStruct *)malloc(sizeof(struct FileListStruct));
        (pfl->next)->prev = pfl;
        pfl = pfl->next;
    }
    pfl->next = 0;

    pfl->fname = (char *)malloc((n+1)*sizeof(char));
    strncpy(pfl->fname,fname,n);
    pfl->fname[n] = '\0';

    pfl->fstatus = (char *)malloc((m+1)*sizeof(char));
    strncpy(pfl->fstatus,fstatus,m);
    pfl->fstatus[m] = '\0';
    string_to_lcase(pfl->fstatus);

    pfl->fform = (char *)malloc((k+1)*sizeof(char));
    strncpy(pfl->fform,fform,k);
    pfl->fform[k] = '\0';
    string_to_lcase(pfl->fform);

    if (strcmp(pfl->fstatus,"old") == 0 || strcmp(pfl->fstatus,"read") == 0
                                        || strcmp(pfl->fstatus,"r") == 0) {
        pfl->fp = fopen(pfl->fname, "r");
        if (pfl->fp == 0) {
            printf("Error in openning file: %s\n",pfl->fname);
            perror("Error in openning file");
            exit(1);
        }
#ifndef NoXDR_format
        else if (strcmp(pfl->fform,"xdr") == 0) {
            xdrstdio_create(&(pfl->xdrs), pfl->fp, XDR_DECODE);
        }
#endif
    }
    else if (   strcmp(pfl->fstatus,"unknown") == 0
             || strcmp(pfl->fstatus,"write") == 0
             || strcmp(pfl->fstatus,"w") == 0) {
        pfl->fp = fopen(pfl->fname, "w");
        if (pfl->fp == 0) {
            printf("Error in openning file: %s\n",pfl->fname);
            perror("Error in openning file");
            exit(1);
        }
#ifndef NoXDR_format
        else if (strcmp(pfl->fform,"xdr") == 0) {
            xdrstdio_create(&(pfl->xdrs), pfl->fp, XDR_ENCODE);
        }
#endif
    }
    else {
        printf("Unknown file status: %s\n",pfl->fstatus);
        perror("Error in openning file");
        exit(1);
    }

    NumFiles++;
    pfl->index = jx + 1;
    *findex = pfl->index;
}

void c_close(const INTEGER *findex) {
    struct FileListStruct *tfl;
 
    searchFile(findex);

    if (fclose(pfl->fp) == 0) {
        free(pfl->fname);
        free(pfl->fstatus);
#ifndef NoXDR_format
        if (strcmp(pfl->fform,"xdr") == 0) {
            xdr_destroy(&(pfl->xdrs));
        }
#endif
        free(pfl->fform);
        if (NumFiles > 1) {
            if (pfl->index == 1) {
                flh = pfl->next;
                flh->prev = 0;
            }
            else {
                (pfl->prev)->next = pfl->next;
                (pfl->next)->prev = pfl->prev;
                pfl->prev = 0;
            }
            tfl = pfl;
            pfl = pfl->next;
            tfl->next = 0;
            free(tfl);
        }
        else {
            pfl = 0;
            free(flh);
        }
        NumFiles--;
        return;
    }
    else {
        printf("Error in closing file: %s\n",pfl->fname);
    }
    perror("Error in closing file");
    exit(1);
}

/* convert the string to low case */
void string_to_lcase(char *s) {
    char *p;
    p=s;
    while (*p != '\0') {
        *p = lcase(*p);
        p++;
    }
}

/* convert to lower case                            */
char lcase(char c) {
    if (c < 'A') {
        return(c);
    }
    else if (c > 'Z') {
        return(c);
    }
    else {
        return (c+'a'-'A');
    }
}

void c_string_padsize(INTEGER *findex, INTEGER *n, INTEGER *m) {
    searchFile(findex);

#ifndef NoXDR_format
    if (strcmp(pfl->fform,"xdr") == 0) {
        *m = XDR_STRING_SIZE(*n) - *n;
    }
    else
#endif
        *m = 0;
}

void c_fseek(INTEGER *findex, INTEGER *offset, INTEGER *whence) {
    searchFile(findex);

    if (*whence == 0) {
        if (!fseek(pfl->fp, (long)(*offset-1), SEEK_SET)) return;
    }
    else {
        if (!fseek(pfl->fp, (long)(*offset-1), SEEK_CUR)) return;
    }
    printf("Error in fseek file pointer for: %s\n",pfl->fname);
    perror("Error in calling fseek");
    exit(1);
}

void c_write_double(const INTEGER *findex, double *data, const INTEGER *n) {
    int j, success;
 
    searchFile(findex);

    success = 1;
#ifndef NoXDR_format
    if (strcmp(pfl->fform,"xdr") == 0) {
        for (j=0; j<*n; j++) {
            if ( !xdr_double(&(pfl->xdrs), &data[j]) ) {
                success = 0;
                break;
            }
        }
    }
    else 
#endif
    if (fwrite( (void *)data, (size_t)(sizeof(double)), (size_t)*n, 
                 pfl->fp ) != *n) {
        success = 0;
    }

    if (success) {
        return;
    }
    else {
        printf("Error in writing double to file: %s\n",pfl->fname);
    }
    perror("Error in writing data");
    exit(1);
}

void c_write_integer(const INTEGER *findex, INTEGER *data, const INTEGER *n) {
    int j, success;
 
    searchFile(findex);
 
    success = 1;
#ifndef NoXDR_format
    if (strcmp(pfl->fform,"xdr") == 0) {
        for (j=0; j<*n; j++) {
            if (!xdr_INTEGER(&(pfl->xdrs), &data[j])) {
                success = 0;
                break;
            }
        }
    }
    else
#endif
    if (fwrite( (void *)data, (size_t)(sizeof(INTEGER)), (size_t)*n, 
                 pfl->fp ) != *n) {
        success = 0;
    }

    if (success) {
        return;
    }
    else {
        printf("Error in writing integer to file: %s\n",pfl->fname);
    }
    perror("Error in writing data");
    exit(1);
}

void c_write_string(const INTEGER *findex, char *data, const INTEGER *n,
                     INTEGER *m, INTEGER data_length) {
    int j;
    char *ctmp;

    searchFile(findex);

#ifndef NoXDR_format
    if (strcmp(pfl->fform,"xdr") == 0) {
        *m = XDR_STRING_SIZE(*n);
        ctmp = (char *)malloc((*n+1)*sizeof(char));
        strncpy(ctmp,data,*n);
        ctmp[*n] = '\0';
        if (xdr_string(&(pfl->xdrs), &ctmp, (u_int)*n)) {
            free(ctmp);
            return;
        }
        else {
            printf("Error in writing string to file: %s\n",pfl->fname);
        }
    }
    else
#endif
    if (fwrite( (void *)data, (size_t)(*n*sizeof(char)), 1, pfl->fp ) == 1) {
        *m = *n;
        return;
    }
    else {
        printf("Error in writing string to file: %s\n",pfl->fname);
    }
    perror("Error in writing string");
    exit(1);
}

void c_read_double(const INTEGER *findex, double *data, const INTEGER *n) {
    int j, success;
 
    searchFile(findex);

    success = 1;
#ifndef NoXDR_format
    if (strcmp(pfl->fform,"xdr") == 0) {
        for (j=0; j<*n; j++) {
            if (!xdr_double(&(pfl->xdrs), &data[j])) {
                success = 0;
                break;
            }
        }
    }
    else
#endif
    if (fread( (void *)data, (size_t)(sizeof(double)), (size_t)*n, 
                pfl->fp ) != *n) {
        success = 0;
    }

    if (success) {
        return;
    }
    else {
        printf("Error in reading double to file: %s\n",pfl->fname);
    }
    perror("Error in reading data");
    exit(1);
}

void c_read_integer(const INTEGER *findex, INTEGER *data, const INTEGER *n) {
    int j, success;
 
    searchFile(findex);

    success = 1;
#ifndef NoXDR_format
    if (strcmp(pfl->fform,"xdr") == 0) {
        for (j=0; j<*n; j++) {
            if (!xdr_INTEGER(&(pfl->xdrs), &data[j])) {
                success = 0;
                break;
            }
        }
    }
    else
#endif
    if (fread( (void *)data, (size_t)(sizeof(INTEGER)), (size_t)*n, 
                pfl->fp ) != *n) {
        success = 0;
    }

    if (success) {
        return;
    }
    else {
        printf("Error in reading integer to file: %s\n",pfl->fname);
    }
    perror("Error in reading data");
    exit(1);
}

void c_read_string(const INTEGER *findex, char *data, const INTEGER *n,
                    INTEGER *m, INTEGER data_length) {
    int j;
    char *ctmp;
 
    searchFile(findex);

#ifndef NoXDR_format
    if (strcmp(pfl->fform,"xdr") == 0) {
        *m = XDR_STRING_SIZE(*n);
        ctmp = (char *)malloc((*n+1)*sizeof(char));

        if (xdr_string(&(pfl->xdrs), &ctmp, (u_int)*n)) {
            strncpy(data,ctmp,*n);
            free(ctmp);
            for (j=*n; j<data_length; j++) data[j] = ' ';
            return;
        }
        else {
            printf("Error in reading string to file: %s\n",pfl->fname);
        }
    }
    else
#endif
    if (fread( (void *)data, (size_t)(*n*sizeof(char)), 1, pfl->fp ) == 1) {
        *m = *n;
        for (j=*n; j<data_length; j++) data[j] = ' ';
        return;
    }
    else {
        printf("Error in reading string to file: %s\n",pfl->fname);
    }
    perror("Error in reading data");
    exit(1);
}

void searchFile(const INTEGER *findex) {
    int i;
 
    if (pfl->index == *findex) {
        return;
    }
    else {
        pfl = flh;
        for (i=0; i<NumFiles; i++) {
            if (pfl->index == *findex) {
                return;
            }
            else {
                pfl = pfl->next;
            }
        }
    }
    printf("Invalid file index: %d\n",*findex);
    perror("Error in searchFile");
    exit(1);
}
