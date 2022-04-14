#include<stdio.h>
#include<math.h>
#include<stdio.h>
#define MAXNUM 10000000//Maximum number of iterations
#define TOL 0.0001//allowable error
int ErrorCheck(int n,double* iniX,double* pX)
{
    double error=0,dif;
    for(int i=0;i<n;i++)
    {
        dif=iniX[i]-pX[i];
        dif=(dif>-dif)?(dif):(-dif);
        error=(dif>error)?dif:error;
    }
    return (error>TOL)?0:1;;//shu ru liang ge xiang liang ,ji suan wu qiong fan shu yi yi xia de wu cha, bing yu rong xu wu cha bi jiao
}
int OrIte(double *pX,int n,double* matA,double *b,double omega)
{
    int iteNum=0;
    int state=0;
    double iniX[n]={25,-300,1050,-1400,630};//yu zhen shi jie bi jiao ji suan wu cha
    for(iteNum;iteNum<=MAXNUM;iteNum++)
    {
        for(int i=0;i<n;i++)
        {
            double temp=0;
            pX[i]*=(1-omega);
            for(int j=0;j<=i-1;j++)temp+=(*(matA+i*n+j))*pX[j];
            for(int j=i+1;j<n;j++)temp+=(*(matA+i*n+j))*pX[j];
            pX[i]+=omega*(b[i]-temp)/(*(matA+i*n+i));
        }
        if(ErrorCheck(n,iniX,pX))
        {
            state=1;
            break;
        }
    }
    if(state)return iteNum;
    else return 0;//ru guo yin wei wu cha man zu yao qiu tui chu, state bei xiu gai wei 1 ,fan hui die dai ci shu
}
int main()
{
    int n=5;//a specified size
    double matA[n][n];//define the matA
    for(int i=0;i<n;i++)//intializing with specified element
    {
        for(int j=0;j<n;j++)
        {
            matA[i][j]=1.0/(i+j+1);
        }
    }
    double b[n]={1,0,0,0,0};
    int m=21;
    int state;//store the number of iterations
    double omega;
    for(int i=0;i<m;i++)
    {
        double pX[n]={0,0,0,0,0};//solution,initial value
        omega=i/10.0;
        state=OrIte(pX,n,*matA,b,omega);
        printf("when omega=%g\n",omega);
        if(state)
        {
            printf("the number of iteration is %d\n",state);
            printf("the approximate solution is\n");
            for(int j=0;j<n;j++) printf("%.6lf ",pX[j]);
        }else printf("exceed the set maximum number of iteration,treated as un-convergent");
        printf("\n\n");
    }
}