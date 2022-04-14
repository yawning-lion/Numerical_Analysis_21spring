#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
//下面函数中输入中都有计算好的数组，是为了避免重复计算部分数据比如节点值，使用时利用区间编号取出需要的值
double Omega(int num,double x,double *pX_K)//定义重心坐标加权中的ω函数，输入等分的份数num和这个ω函数的所有零点信息，输出函数值
{
    double y_Omega=(double)1;//为最后要输出的函数值
    for(int k=0;k<num+1;k++)
    {
        y_Omega*=(x-pX_K[k]);//利用循环累乘来计算ω函数
    }
    return y_Omega;
}
//数组pX_K存储节点处所有变量值，数组pX_M储存需要计算函数值的变量，E_Y存储节点处函数值
double Lagrange(int num,double *pX_K,double x)//定义拉格朗日插值函数，使用重心坐标加权计算，输入等分的份数num和变量值x，输出函数值,变量值利用x_k从pX_M中取出
{
    double yLag=(double)0;//最后要输出的函数值
    double pX_I[num+1];//为方便重心加权公式计算的方便，重心加权计算中的分母具有ω函数的形式
    memcpy(pX_I,pX_K,(num+1)*sizeof(double));
    for(int i=0;i<num+1;i++)
    {
        pX_I[i]=x;//将原pX_K中第i+1个值替换为x即为公式分母omega函数的所有零点信息
        yLag-=((double)exp((float)-2*pX_K[i]))*Omega(num,x,pX_K)/Omega(num,pX_K[i],pX_I);//根据定义，Omega(num,pX_K[i],pX_I)即为计算公式中的分母部分
        pX_I[i]=pX_K[i];//将pX_I[i]还原，用于下一步计算
    }
    free(pX_I);
    return yLag;
}

double PieceLinear(int x_k,double *pX_M,double *pX_K)//x_k用于确定区间
{
    double yLin;
    double x=pX_M[x_k];//取出变量值
    yLin=(double)exp(-2*pX_K[x_k])+(exp(-2*pX_K[x_k])-exp(-2*pX_K[x_k+1]))*(-x/pX_K[1]+x_k);//为减小误差，直接用步长做计算，pX_K[1]大小即为步长
    return yLin;
}

int Make_M(int num,double *pX_K,double *E_Y,double *M_DIF)//M_DIF数组用来储存节点一阶导数值
{
    double ALPHA=1.0/2.0;//α为定常0.5
    double BETA[num];
    BETA[0]=(-2)*E_Y[0];
    double A[num];
    double B[num];
    A[0]=0;
    B[0]=BETA[0];//利用固支边界条件
    for(int i=1;i<num;i++)
    {
        BETA[i]=3*ALPHA*(E_Y[i+1]-E_Y[i-1])/pX_K[1];//简化后的计算公式,顺次计算各β值
        A[i]=-ALPHA/(2+ALPHA*A[i-1]);
        B[i]=(BETA[i]-ALPHA*B[i-1])/(2+ALPHA*A[i-1]);//逐次计算递推公式中的系数值
    }
    M_DIF[num]=-2*E_Y[num];//利用边界条件得到M_DIF最后一个值
    for(int i=num-1;i>-1;i--)
    {
        M_DIF[i]=A[i]*M_DIF[i+1]+B[i];//利用递推公式计算M
    }
    free(A);
    free(B);
    free(BETA);
    return 0;
}
double Spline(int num,int x_k,double *pX_M,double *pX_K,double *M_DIF,double *E_Y)//输入各节点一阶导数值，区间编号和包含各数据的数组，
{
    double ySpl=0;//为最后返回的函数值
    double x=pX_M[x_k];//变量值
    double X0=pX_K[x_k],X1=pX_K[x_k+1],H=pX_K[1];//H为步长，X0为第i个节点值，X1为第i+1个
    double Y0=E_Y[x_k],Y1=E_Y[x_k+1];
    ySpl=(1+2*(x-X0)/H)*(x-X1)*(x-X1)*Y0/(H*H)+(1-2*(x-X1)/H)*(x-X0)*(x-X0)*Y1/(H*H)+(x-X0)*M_DIF[x_k]*(x-X1)*(x-X1)/(H*H)+(x-X1)*M_DIF[x_k+1]*(x-X0)*(x-X0)/(H*H);//利用等分区间的特殊性简化计算
    return ySpl;
}

int main()
{
    int num;
    printf("please enter the value of n\n");
    scanf("%d",&num);//读入区间等分的份数
    if(num<0||num>700)
    {
        printf("illegal input,please run again");
        exit(0);
    };
    double *pX_K;
    double *pX_M;//储存用于插值的节点变量值和需要计算函数值的变量值
    double E_Y[num+1];//存储节点函数值
    pX_K=(double*)malloc((num+1)*sizeof(double));
    pX_M=(double*)malloc((num+1)*sizeof(double));//有n+1格，实际上只用到前n格
    for(int k=0;k<num+1;k++)//生成节点数组、变量、函数值数组
    {
        pX_K[k]=(double)(6*k)/(double)num;
        pX_M[k]=(double)(6*k+3)/(double)num;
        E_Y[k]=exp(-2*pX_K[k]);
    }
    double M_DIF[num+1];//存储样条插值各节点一阶导数值
    Make_M(num,pX_K,E_Y,M_DIF);//生成样条插值各节点一阶导数值
    double M_error[3]={0,0,0};//定义最大误差
    double error[3]={0,0,0};
    FILE *pFile;
    pFile=fopen("output.txt","a");//本程序一次只能输入一个区间数并输出，每次运行都在原有文档内继续写入数据而不覆盖
    fprintf(pFile,"当n的值为%d时\n",num);
    fprintf(pFile,"|  变量值  |  拉格朗日插值 |   分段线性插值 |  三次样条插值  |    原函数值    |\n");//输出格式为按行输出，每行依次输出某一个变量下各插值函数之值
    for(int i=0;i<num;i++)
    {
        double yLag=Lagrange(num,pX_K,pX_M[i]);
        double yLin=PieceLinear(i,pX_M,pX_K);
        double ySpl=Spline(num,i,pX_M,pX_K,M_DIF,E_Y);
        error[0]=(double)abs(yLag-exp(-2*pX_M[i]));
        error[1]=(double)abs(yLin-exp(-2*pX_M[i]));
        error[2]=(double)abs(ySpl-exp(-2*pX_M[i]));
        M_error[0]=(M_error[0]>error[0])?M_error[0]:error[0];
        M_error[1]=(M_error[1]>error[1])?M_error[1]:error[1];
        M_error[2]=(M_error[2]>error[2])?M_error[2]:error[2];//依次对比，较大的替换为M_error的值
        fprintf(pFile,"| x=%.3f | %12g |\t%12g |\t%12g |\t%12g |\n",pX_M[i],yLag,yLin,ySpl,exp(-2*(pX_M[i])));//在变量值取定下，依次输出各插值函数值
    }
    fprintf(pFile,"三种插值函数的最大绝对误差分别为\n");
    fprintf(pFile,"%12g\t%12g\t%12g\n\n",M_error[0],M_error[1],M_error[2]);//依次输出最大误差值
    free(E_Y);
    free(M_DIF);
    free(pX_K);
    free(pX_M);
    printf("\nEND");
    return 0;
}