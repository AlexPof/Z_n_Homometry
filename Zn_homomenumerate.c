#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
* EnsZN is the structure definition for subsets of Z_n
*/
typedef struct EnsZN {
  int* A0;
  char* hash;
  int N;
} EnsZN;

void enumerate_homometric(int N, int P, FILE* output_file);
void initialize_EnsZN(EnsZN* ensemble, int* sequence, int N);
void get_IV(EnsZN* X,EnsZN* IV);
int is_equal(EnsZN* X,EnsZN* Y);

void niceprint_EnsZN(EnsZN* ensemble, char* str_rep);
int is_trivially_related(EnsZN* X,EnsZN* Y);
int mod(int a, int b);
void next_kbit_seq(int* seq,int* nextseq,int N);
void binaryprint(unsigned int a,int N);

/*
* main takes 3 arguments
*
* N: order of the Z_N dihedral group
* P: cardinality of the subsets to be considered
* output_file: the output for the enumeration
*/

int main(int argc, char *argv[]) {
  int N,P;
  FILE* output_file=NULL;

  if (argc<3) {
    printf("Not enough arguments - Need at least 4\n");
    exit(1);
  }

  N = atoi(argv[1]);
  P = atoi(argv[2]);
  output_file=fopen(argv[3],"w");
  if (output_file==NULL) {
    printf("Error creating the output file...\n");
    exit(1);
  }

  enumerate_homometric(N, P, output_file);

  fclose(output_file);
  exit(0);
}

void enumerate_homometric(int N, int P, FILE* output_file) {
    int total_ensembles=0;
    int i,j,flag,c;
    EnsZN *ensembles,*ensembles_left_IV;
    EnsZN X,X_left_IV;
    int sequence[N],next_sequence[N];
    char str_rep[200]="";
    int p,sign,endflag=0;

    printf("Building the collection of Z_%d subsets of cardinality %d...\n",N,P);
    for(i=0;i<N;i++) {
      if(i<P)
        sequence[i]=1;
      else
        sequence[i]=0;
    }

    while(!endflag) {
      initialize_EnsZN(&X, sequence, N);
      get_IV(&X,&X_left_IV);
        // If this is the first subset found, we initialize the ensemble table
        if (total_ensembles==0) {
          total_ensembles++;
          ensembles = (EnsZN*)malloc(sizeof(EnsZN));
          initialize_EnsZN(ensembles, sequence, N);

          ensembles_left_IV = (EnsZN*)malloc(sizeof(EnsZN));
          get_IV(ensembles,ensembles_left_IV);
        } else {
          // Otherwise, we check if the IV already exists in our enumeration
          flag=0;
          for (i=0;i<total_ensembles;i++) {
            if(is_equal(&X_left_IV,&ensembles_left_IV[i])) {
              // Most commonly, the interval vector already exists because it is
              // a translate or an inversion of a previously enumerated ensemble,
              // in which case we flag it so we now we don't need it.
              if (is_trivially_related(&X,&ensembles[i])) {
                flag=1;
              }
            }
          }

          if (flag==0) {
            // If we have not flagged it before, it means that we have an ensemble
            // with a common IV, but which is not trivially related to the previous
            // ones so we add it to the enumeration.
            total_ensembles++;
            ensembles = (EnsZN*)realloc(ensembles,total_ensembles*sizeof(EnsZN));
            ensembles_left_IV = (EnsZN*)realloc(ensembles_left_IV,total_ensembles*sizeof(EnsZN));

            initialize_EnsZN(&ensembles[total_ensembles-1], sequence, N);
            get_IV(&ensembles[total_ensembles-1],&ensembles_left_IV[total_ensembles-1]);
          }
        }

      endflag=1;
      for(i=0;i<P;i++)
        endflag &= sequence[N-1-i];
      next_kbit_seq(sequence,next_sequence,N);
      for(i=0;i<N;i++)
        sequence[i]=next_sequence[i];
    }

    c=0;
    for (i=0;i<total_ensembles;i++) {
      for (j=i+1;j<total_ensembles;j++) {
        if(is_equal(&ensembles_left_IV[i],&ensembles_left_IV[j]) && !is_equal(&ensembles[i],&ensembles[j])) {
          fprintf(output_file,"===== %d ======\n",c);
          fprintf(output_file,(&ensembles[i])->hash);
          fprintf(output_file," - ");
          fprintf(output_file,(&ensembles[j])->hash);
          fprintf(output_file,"\n");
          niceprint_EnsZN(&ensembles[i], str_rep);
          fprintf(output_file,str_rep);
          fprintf(output_file,"\n");
          niceprint_EnsZN(&ensembles[j], str_rep);
          fprintf(output_file,str_rep);
          fprintf(output_file,"\n");
          c+=1;
        }
      }
    }
}

void get_IV(EnsZN* X,EnsZN* IV) {
  /*
  * Function:  get_IV
  * --------------------
  * Calculates the Interval Vector (IV) of a subset of Z_n
  *
  * X: the Z_n subset, whose left IV is calculated
  *
  *  returns: the IV as a subset 'IV' of Z_n
  */
  int i,j;

  IV->N = X->N;
  IV->A0 = calloc(X->N,sizeof(int));

  for(i=0;i<X->N;i++) {
    IV->A0[i] = 0;
  }

  for(i=0;i<X->N;i++) {
    for(j=0;j<X->N;j++) {
      IV->A0[mod((i-j),X->N)] += X->A0[i]*X->A0[j];
    }
  }
}


int is_equal(EnsZN* X,EnsZN* Y){
  int i;

  for(i=0;i<X->N;i++) {
    if (X->A0[i]!=Y->A0[i])
      return 0;
  }
  return 1;
}

///////////////////////////////////////////////////////

int is_trivially_related(EnsZN* X,EnsZN* Y){
  /*
  * Function:  is_trivially_related
  * --------------------
  * Checks if two Z_n subsets X and Y are trivially related, i.e. either
  * translated or inversionally related.
  *
  * X,Y: the D_2n subsets to be checked
  *
  *  returns: True if the subsets X and Y are trivially related, False otherwise
  */
  int p,i,N,c,d;

  N=X->N;
  for(p=0;p<N;p++) {
    c=0;
    d=0;
    for(i=0;i<N;i++) {
      if (X->A0[i] == Y->A0[ mod((i+p),N)])
        c++;
      if (X->A0[i] == Y->A0[ mod((p-i),N)])
        d++;
    }
    if (c==N || d==N)
      return 1;
  }

  return 0;
}


void initialize_EnsZN(EnsZN* ensemble, int* sequence, int N) {
  /*
  * Function:  initialize_EnsZN
  * --------------------
  * Initialize a Z_n subset with the given values
  *
  * ensemble: the D_2n subset to be initialized
  * a: the integer representation of the D_2n subset
  * N: the order n of D_2n
  *
  *  returns: None
  */
  int i,j;
  int countX;
  int numbytes=(N/4)+2;

  ensemble->N = N;
  ensemble->A0 = calloc(N,sizeof(int));
  ensemble->hash = calloc(numbytes,sizeof(char));

  for(i=0;i<numbytes;i++) {
    ensemble->hash[i]=65;
  }
  ensemble->hash[numbytes-1]=0;
  for(i=0;i<N;i++) {
    ensemble->hash[i/4] = ensemble->hash[i/4] | (sequence[i]<<((i%4)+1));
  }

  for(i=0;i<N;i++) {
    ensemble->A0[i] = sequence[i];
  }
}

void niceprint_EnsZN(EnsZN* ensemble, char* str_rep) {
  /*
  * Function:  niceprint_EnsZN
  * --------------------
  * Print an interpretable version of a D_2n subset
  *
  * ensemble: the D_2n subset to be printed
  * str_rep: the destination string. Prints an element (g,h) of the subset
  *          as g+ if h=1, g- otherwise.
  *
  *
  *  returns: None.
  *
  */
  int i;
  char* temp_str;

  strcpy(str_rep,"");
  // 200 is the max size of the string
  // Not the better way to code it, though
  snprintf(str_rep, 200, "%s%s", str_rep, "{");
  for(i=0;i<ensemble->N;i++) {
    if (ensemble->A0[i])
      snprintf(str_rep, 200, "%s%d,", str_rep, i);
  }
  snprintf(str_rep, 200, "%s%s", str_rep, "}");
}

void next_kbit_seq(int* seq,int* nextseq,int N) {
  int i,smallest,next_smallest,c,nc;
  int ripples[N],ones[N],final[N];

  for(i=0;i<N;i++) {
    ripples[i]=seq[i];
    ones[i]=0;
    final[i]=0;
  }

  smallest=0;
  while(seq[smallest]==0)
    smallest+=1;

  ripples[smallest]=0;
  c=1;
  for(i=smallest+1;i<N;i++) {
    nc = c & ripples[i];
    ripples[i] = (ripples[i]+c)%2;
    c=nc;
  }

  next_smallest=0;
  while(ripples[next_smallest]==0)
    next_smallest+=1;

  for(i=0;i<(next_smallest-smallest-1);i++)
    ones[i]=1;
  for(i=0;i<N;i++)
    nextseq[i] = ripples[i] | ones[i];
}


int mod(int a, int b) {
  /*
  * Function:  mod
  * --------------------
  * Arithmetic modulo of a by b
  *
  *  returns: the value of mod(a,b).
  */
  return (a >= 0 ? a % b : b - (-a) % b);
}

void binaryprint(unsigned int a,int N) {
  /*
  * Function:  binaryprint
  * --------------------
  * Print the binary representation of an integer a considered as a word
  * of length N. For debug purposes.
  */
  int i;
  for(i=0;i<N;i++)
  printf("%d",(a>>i)&1);
  printf("\n");
}
