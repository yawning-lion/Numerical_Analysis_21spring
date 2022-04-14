#include<stdio.h>
#include<math.h>
#include<stdio.h>
#define MAXNUM 64
int LDLTDecompose(double matA[][MAXNUM],double matL[][MAXNUM],double *matD,int n)
{
    double matT[n];//temp matrix
    matD[0]=matA[0][0];
    for(int j=1;j<n;j++)
    {
        for(int k=0;k<j;k++)
        {
            matT[k]=matA[k][j];
            for(int i=0;i<k;i++)matT[k]-=matL[k][i]*matT[i];
            matL[j][k]=matT[k]/matD[k];
        }
        matD[j]=matA[j][j];
        for(int i=0;i<j;i++)matD[j]-=matL[j][i]*matT[i];
    }//execute LDLT algorithm
    return 0;
}
int SolveB(double matL[][MAXNUM],double *matD,double matB[][MAXNUM],int n,int m)
{
    double yTemp[n];
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)//compute y
        {
            yTemp[j]=matB[j][i];
            for(int l=0;l<j;l++)yTemp[j]-=matL[j][l]*yTemp[l];
        }
        for(int l=0;l<n;l++)yTemp[l]=yTemp[l]/matD[l]; 
        for(int j=n-1;j>-1;j--)
        {
            matB[j][i]=yTemp[j];
            for(int l=j+1;l<n;l++)matB[j][i]-=matL[l][j]*matB[l][i];//for reduce storage, store the solution in matB straightly
        }
    }
    return 0;
}
int OutputMat(double matB[][MAXNUM],int n,int m,FILE *pOutput)//display all the element in a file
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++) fprintf(pOutput,"%g ",matB[i][j]);
        fprintf(pOutput,"\n");
    }
    return 0;
}
int main()
{
    int n;//matrix's size
    FILE *pInput;
    pInput=fopen("Input.txt","r");
    fscanf(pInput,"%d",&n);// to determine the row of matrix
    if(n>64)
    {
        printf("overflow,need a smaller n");
        return 0;
    }
    double matA[MAXNUM][MAXNUM],matL[MAXNUM][MAXNUM];//define unsolved matrix A, Lower triangular matrix L
    double matD[MAXNUM];//define Diagonal matrix D
    for(int i=0;i<n;i++)//initialize all the matrix
    {
        matD[i]=0.0;
        for(int j=0;j<n;j++)matA[i][j]=matL[i][j]=0.0;
        matL[i][i]=1.0;
    }
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)fscanf(pInput,"%lf",&matA[i][j]);
    LDLTDecompose(matA,matL,matD,n);
    FILE *pOutput;
    pOutput=fopen("Output.txt","w");
    fprintf(pOutput,"The  element of the matrixL are as follows\n");
    OutputMat(matL,n,n,pOutput);
    fprintf(pOutput,"\n");
    fprintf(pOutput,"The diagonal element of the matrixD are as follows\n");
    for(int i=0;i<n;i++)fprintf(pOutput,"%g ",matD[i]);
    fprintf(pOutput,"\n");
    fprintf(pOutput,"\n");
    int m=0;
    fscanf(pInput,"%d",&m);//determine the column of B
    if(m>64)
    {
        printf("overflow,need a smaller m");
        return 0;
    }
    double matB[MAXNUM][MAXNUM];
    for(int i=0;i<n;i++)for(int j=0;j<m;j++)fscanf(pInput,"%lf",&matB[i][j]);//input the element of matrix B
    SolveB(matL,matD,matB,n,m);
    fprintf(pOutput,"The element of the solution are as follows\n");
    OutputMat(matB,n,m,pOutput);
    printf("please check the output in file output.txt in the folder code");
}