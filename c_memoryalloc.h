#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

long *vector_long(long nl,long nh){
        long *v;
        v=(long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
        return v-nl+NR_END;
}

void free_vector_long(long *v,long nl,long nh){
       free((FREE_ARG) (v+nl-NR_END));
}

long **matrix_long(long nrl,long nrh,long ncl,long nch){
        long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        long **m;

        m=(long **)malloc((size_t)((nrow+NR_END)*sizeof(long*)));
        m+=NR_END;
        m-=nrl;

        m[nrl]=(long *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(long)));
        m[nrl]+=NR_END;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
        return m;
}

void free_matrix_long(long **m,long nrl,long nrh,long ncl,long nch){
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

double *vector_double(long nl,long nh){
        double *v;
        v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
        return v-nl+NR_END;
}

void free_vector_double(double *v,long nl,long nh){
       free((FREE_ARG) (v+nl-NR_END));
}

double **matrix_double(long nrl,long nrh,long ncl,long nch){
        long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        double **m;

        m=(double **)malloc((size_t)((nrow+NR_END)*sizeof(double*)));
        m+=NR_END;
        m-=nrl;

        m[nrl]=(double *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
        m[nrl]+=NR_END;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
        return m;
}

void free_matrix_double(double **m,long nrl,long nrh,long ncl,long nch){
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}
