#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define MAXNUM 64
const double PIE=3.14159265358979323846264338328;
const double A=0;
const double B=PIE/2;

typedef struct Oval
{
    double majorAxis;
    double minorAxis;
}COval;//define struction to restore all information of a oval
typedef COval* pCOval;

double OvalFx(double Theta,pCOval oval)//Give an expression for the integrand
{
    double yFx;
    double ySin=sin(Theta);//The value of sin theta is calculated in advance， to reduce repeated calculation
    double maj=oval->majorAxis;
    double min=oval->minorAxis;
    double a=min*min;
    double b=maj*maj-a;
    yFx=4*sqrt(a+b*ySin*ySin);
    return yFx;
}
pCOval IniOva(double majorAxis,double minorAxis)
{
    pCOval oval;
    oval=(pCOval)malloc(sizeof(COval));
    oval->majorAxis=majorAxis;
    oval->minorAxis=minorAxis;//initialize a oval
    return oval;
}
int Lengthen_T0(int k,double h,int n,double integralT[][MAXNUM],pCOval oval)//calculate the new T when Add sub-divisions to n
{
    double tempF=0.0;//used for calculating T[0][k];
    for(int i=1;i<n+1;i++) tempF+=OvalFx(A+(2*i-1)*h,oval);
    integralT[0][k]=integralT[0][k-1]/2+h*tempF;
    return 0;
}

double CalErr(double a,double b)//calculate error
{
    double c=a-b;
    return (c>0)?(c):(-c);
}

double Lengthen_Tm(int k,double integralT[][MAXNUM],double EPSILON)//apply romberg
{
    double error;
    for(int m=1;m<k+1;m++)
    {
        double temp=pow(4,m);
        integralT[m][k-m]=(temp*integralT[m-1][k-m+1]-integralT[m-1][k-m])/(temp-1);
        error=CalErr(integralT[m][0],integralT[m-1][0]);
        if(error<EPSILON) return integralT[m][0];//If the error meets the requirements, exit and return the integral value
    }
    return -1;//If the error  does't meets the requirements at last,then return -1 as sign
}
int main()
{
    double minorAxis;
    double majorAxis;//define two axis of the oval
    double EPSILON;//difine the allowable error
    printf("please enter the value of the majorAxis:\n");
    scanf("%lf",&majorAxis);
    printf("please enter the value of the minorAxis:\n");
    scanf("%lf",&minorAxis);
    printf("please enter the value of the allowable error:\n");
    scanf("%lf",&EPSILON);
    //define important parameter;
    double h=(B-A)/2;//define the length the step
    int k=1;
    int n=1;//control the upper bound of a cycle
    pCOval oval;
    oval=IniOva(majorAxis,minorAxis);
    double integralT[MAXNUM][MAXNUM];//define a matrix to restore all the value of integral when running Romberg
    integralT[0][0]=h*(OvalFx(A,oval)+OvalFx(B,oval));
    double INTEGRAL=0;
    double state;
    do
    {
        Lengthen_T0(k,h,n,integralT,oval);
        state=Lengthen_Tm(k,integralT,EPSILON);
        if(state<0)
        {
            h=h/2;
            n=n*2;
            k++;//if error is too large, apply romberg once again
        }else
        {
            INTEGRAL=state;
        }
    }while(state<0);
    FILE *pOutput;
    pOutput=fopen("output.txt","a");
    for(int i=0;i<k+1;i++) fprintf(pOutput,"当外推次数为m=%d时，积分结果为%.12lf\n",i,integralT[i][0]);
    fprintf(pOutput,"当外推次数为m=%d，积分值符合精度要求，保留12位小数的结果为%.12lf",k,INTEGRAL);
    free(oval);
}

