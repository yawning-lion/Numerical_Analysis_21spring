#include<stdio.h>
#include<math.h>
#include<stdio.h>
# define EPS 0.0000001
int RenewJ(double* jocA,double* arrX,int n,int m)
{
    n=m=3;
    for(int i=0;i<m;i++) *(jocA+i)=2*arrX[i];
    *(jocA+1*n+0)=4*arrX[0];
    *(jocA+1*n+1)=2*arrX[1];
    *(jocA+1*n+2)=-4;
    *(jocA+2*n+0)=6*arrX[0];
    *(jocA+2*n+1)=-4;
    *(jocA+2*n+2)=2*arrX[2];
    return 0;
}// make jocobi matrix of A from x
int RenewB(double* arrX,double* arrB,int n)
{
    double x1=arrX[0];
    double x2=arrX[1];
    double x3=arrX[2];
    arrB[0]=-(x1*x1+x2*x2+x3*x3-1);
    arrB[1]=-(2*x1*x1+x2*x2-4*x3);
    arrB[2]=-(3*x1*x1-4*x2+x3*x3);
    return 0;
}// make the right side of the equation
int RenewX(double* jocA,double* arrB,double* arrX,int n,int m)
{
    double w=0;
    for(int i=0;i<n-1;i++)
    {
        for(int j=i+1;j<n;j++)
        {
            w=(*(jocA+n*j+i))/(*(jocA+n*i+i));
            for(int k=i;k<m;k++)
            {
                *(jocA+n*j+k)=*(jocA+n*j+k)-*(jocA+n*i+k)*w;
            }
            arrB[j]=arrB[j]-arrB[i]*w;
        }
    }
    double deltX[n];
    for(int i=n-1;i>-1;i--)
    {
        double temp=0;
        for(int j=i+1;j<n;j++)
        {
            temp+=*(jocA+n*i+j)*deltX[j];
        }
        deltX[i]=(arrB[i]-temp)/ *(jocA+n*i+i);
    }// Calculate delta x by Gaussian elimination method
    for(int i=0;i<n;i++) arrX[i]+=deltX[i];// renew x
    return 0;
}
int ErrorCheck(int n,double* lastX,double*arrX )
{
    double error=0,dif;
    for(int i=0;i<n;i++)
    {
        dif=lastX[i]-arrX[i];
        dif=(dif>-dif)?(dif):(-dif);
        error=(dif>error)?dif:error;
    }
    return (error>EPS)?0:1;;//shu ru liang ge xiang liang ,ji suan wu qiong fan shu yi yi xia de wu cha, bing yu rong xu wu cha bi jiao
}
int main()
{
    int n=3,m=3;
    int k=0;
    double jocA[n][m];
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) jocA[i][j]=0;
    double arrX[n]={1.0,1.0,1,0};
    double arrB[n]={0,0,0};
    double lastX[n]={0,0,0};
    while(!ErrorCheck(n,lastX,arrX))
    {
        for(int i=0;i<n;i++) lastX[i]=arrX[i];
        RenewB(arrX,arrB,n);
        RenewJ(*jocA,arrX,n,m);
        RenewX(*jocA,arrB,arrX,n,m);
        k++;
    }
    printf("solution=\n");
    for(int i=0;i<n;i++) printf("%.7lf ",arrX[i]);
    printf("\nThe number of iterations is %d",k);
}