#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265


/*...................... Define Complex Number ...........................*/
typedef struct complex
{
   double realpart, imaginary;
}Complex;

/*................... checking if N is power of 2 ........................*/
int CheckPowerOfTwo(int N){
    double c=N;
    while (c>1)
        c=c/2;
    if(c<1||N==1){
        printf("N must be power of 2:\n");
        return 0;
    }
    else
        return 1;
    }

/*........................ Evaluate log2(N) .............................*/
int EvaluateLogN(int N){
    int count=0;
    while(N!=1){
        N=N/2;
        count++;
    }
    return count;
}

/*........... Evaluate the new index for the input x[N] ................*/
void EvaluateNewInputIndex(int N,int log2N,int index[]){
    int a[N];//a[N] holding the binary of a number 0->N
    int i,j; //for loops
    for(i=0;i<N;i++){
        index[i]=0;
    }
    for(i=0;i<N;i++){
            int k=i;     //k the numbers from 0->N-1
            for(j=0;j<log2N;j++){
                a[j]=k%2;    //store the binary of k in the array a[]
                k=k/2;
            }
            j=log2N-1;
            k=0;
            for(j;j>=0;j--){
                index[i]+=pow(2,k)*a[j];  //reverse the binary number and transfer it back to decimal
                k++;
            }
    }
}


/*................. Evaluate the output of FFT ......................*/
void EvaluateFFT(int N,int log2N,int index[],Complex x[],Complex X[]){
    int k=2;  //number of operation in a butterfly in the first stage
    int m=1;  //half of k
    int i,j,n; //for loops
    double W;
    for(i=0;i<log2N;i++){
        W=((2*180)/k)*(PI/180);   //W=2*pi/N   degree
        int c=0;   //index for the output
        int f=0;   //the start of a butterfly
        for(n=0;n<N/k;n++){
            //first equation X[k]=f1(k)+ WnK * f2(k)
            for(j=f;j<m+f;j++){
                X[c].realpart=x[index[j]].realpart+x[index[j+m]].realpart*cos(W*j)+x[index[j+m]].imaginary*sin(W*j);
                X[c].imaginary=x[index[j]].imaginary-x[index[j+m]].realpart*sin(W*j)+x[index[j+m]].imaginary*cos(W*j);
                c++;
            }
            //second equation X[k+N/2]=f1(k)- WnK * f2(k)
            for(j=f;j<m+f;j++){
                X[c].realpart=x[index[j]].realpart-x[index[j+m]].realpart*cos(W*j)-x[index[j+m]].imaginary*sin(W*j);
                X[c].imaginary=x[index[j]].imaginary+x[index[j+m]].realpart*sin(W*j)-x[index[j+m]].imaginary*cos(W*j);
                c++;
            }
            f=c;  //to start the next butterfly in the same stage
        }
        for(j=0;j<N;j++){
            x[j].realpart=X[j].realpart;      //the real part of input for the next stage
            x[j].imaginary=X[j].imaginary;     //the imaginary part of input for the next stage
            index[j]=j;
        }
        k*=2;            //for the next stage of butterfly
        m*=2;
    }
}

/*...................... Printing Discrete Fourier Transform X[K] ...........................*/
void PrintDFT(int N, Complex X[])
{   int i;

    /* Printing X[K] & Handling Real & Img if they equals Zeros */
    for(i = 0; i < N; i++)
    {
        if(i == 0)
        printf("X[K} { ");

        if(i == N-1)
        {
            if ( (int)(X[i].imaginary * 100) == 0 )
                printf("%0.3lf, ",X[i].realpart);
            else if( (int)(X[i].realpart * 100) == 0 )
            {
                if(X[i].imaginary > 0)
                printf("%0.3lfi ",X[i].imaginary);
            else
                printf("%0.3lfi ",X[i].imaginary);
            }
            else if(X[i].imaginary > 0)
                printf("%0.3lf + %0.3lfi ",X[i].realpart,X[i].imaginary);
            else
                printf("%0.3lf %0.3lfi ",X[i].realpart,X[i].imaginary);
                printf("}");
        }
        else{
            if ( (int)(X[i].imaginary * 100) == 0 )
                printf("%0.3lf, ",X[i].realpart);
            else if( (int)(X[i].realpart * 100) == 0 )
            {
                if(X[i].imaginary > 0)
                printf("%0.3lfi, ",X[i].imaginary);
                else
                printf("%0.3lfi, ",X[i].imaginary);
            }
            else if(X[i].imaginary > 0)
                printf(" %0.3lf + %0.3lfi, ",X[i].realpart,X[i].imaginary);
            else
                printf(" %0.3lf %0.3lfi, ",X[i].realpart,X[i].imaginary);
        }

    }
    printf("\n");
}

/*...................... Printing Magnitude of X[K] ...........................*/
void PrintMagDFT(int N,double MagDFT[]){

        /* Printing Mag[X[K]] */
    int i;
    for(i = 0; i < N; i++)
    {
        if(i == 0)
            printf("Mag[X[K]] { ");
        if(i == N-1)
        {
            printf("%0.3lf ",MagDFT[i]);
            printf("}");
        }
        else{
            printf("%0.3lf ,",MagDFT[i]);
        }
    }
    printf("\n");
}

/*...................... Printing Phase of X[K] ...........................*/
void PrintPhaseDFT(int N,double PhaseDFT[])
{
    /* Printing Phase[X[K]] */
    int i;
    for(i = 0; i < N; i++)
    {
        if(i == 0)
            printf("Phase[X[K]] { ");
        if(i == N-1)
        {
            printf("%0.3lf ",PhaseDFT[i]);
            printf("}");
        }
        else{
            printf("%0.3lf, ",PhaseDFT[i]);
        }

    }
}

/*...................... Evaluate Magnitude of X[K] ...........................*/
void EvaluateMagintudeDFT(int N,Complex X[],double MagDFT[]){
    int i;

    for(i = 0; i < N; i++)
    {
        MagDFT[i] = sqrt( ( X[i].realpart * X[i].realpart ) + ( X[i].imaginary * X[i].imaginary ) );
    }
}

/*...................... Evaluate Phase of X[K] ...........................*/
void EvaluatePhaseDFT(int N,Complex X[],double PhaseDFT[]){
        int i;

    for(i = 0; i < N; i++)
    {
        if ( (int)(X[i].imaginary * 100) == 0)
        {
            if(X[i].realpart > 0)
                PhaseDFT[i] = 0;
            else
                PhaseDFT[i] = 180;
        }
        else if ( (int)(X[i].realpart * 100) == 0)
            {
                if(X[i].imaginary > 0)
                    PhaseDFT[i] = 90;
                else
                    PhaseDFT[i] = -90;
            }
        else if(X[i].imaginary > 0 && X[i].realpart > 0)
            PhaseDFT[i] = atan( (X[i].imaginary) / (X[i].realpart) ) * 180 / PI;
        else if (X[i].imaginary > 0 && X[i].realpart < 0)
            PhaseDFT[i] = 180 - ( atan( (X[i].imaginary) / (-1 * (X[i].realpart) ) ) * 180 / PI );
        else if (X[i].imaginary < 0 && X[i].realpart < 0)
            PhaseDFT[i] = -180 + ( atan ( (-1 * (X[i].imaginary) ) / (-1 * (X[i].realpart) ) ) * 180 / PI  );
        else if (X[i].imaginary< 0 && X[i].realpart> 0)
            PhaseDFT[i] = - ( atan( (-1 * (X[i].imaginary) ) / (X[i].realpart) )  * 180 / PI );
    }
}
int main(int argc, char const *argv[])
{
    int N,log2N;
    int i,j;
    printf("How many N-Points do you want to calculate FFT? \n");
    while(1){
        scanf("%d",&N);                 // Reading Number of Points
        if(CheckPowerOfTwo(N))          //check if N is power of 2
            break;
    }
    int index[N];                       //index[] holding the new indexes for the input signal
    Complex X[N],x[N];                  //x[N]the input signal ,X[N] the output DFT
    double MagDFT[N],PhaseDFT[N];             //M[N] the magnitude P[N] the phase


    // Reading Input x[n] in Sequence Form
    printf("Enter the sequence form of x(N): Ex: 1 2 3 4 \n");
    for(i=0;i<N;i++){
        scanf("%lf",&x[i].realpart);
        x[i].imaginary=0;
    }

    log2N=EvaluateLogN(N);                                      //Evaluate log2(N)
    EvaluateNewInputIndex(N,log2N,index);                       //calculate the index of x[N]
    EvaluateFFT(N,log2N,index,x,X);                             //calculate the output X[N] with radix-2
    EvaluateMagintudeDFT(N, X, MagDFT);                         // Evaluate Magnitude
    EvaluatePhaseDFT(N,X,PhaseDFT);                             // Evaluate Phase
    PrintDFT(N, X);                                             // Printing X[K]
    PrintMagDFT(N, MagDFT);                                     // Printing Magnitude
    PrintPhaseDFT(N, PhaseDFT);                                 // Printing Phase


    return(0);
}
