//���޹p�g4�ඥ�t�ΡA�ɶ��b����
#include <iostream>
#include <cstdlib>
#include <cctype>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;
#include "dislin.h"

double ReRabi(double &x);//�߽ĥ]���u���(�곡)
double ImRabi(double &x);//�߽ĥ]���u���(�곡)
void fun(double ,double ,double ,int );//�һݸѪ��p�ߤ�{��
void solve(double*,double &,int &,double &);//�t��k

double dx1,dx2,gamma_1,gamma_13,gamma_14,gamma_2,gamma_23,gamma_24,gamma_p,gamma_c;
//dx1:���߽Įɪ��B�e�Cdx2:�S���߽Įɪ��B�e�Cgamma_i�G�ඥi���I�ܲv�Cgamma_ij�G�ඥi�I�ܨ�ඥj���I�ܲv�C
//gamma_p�Mgamma_c����Arelaxation rate�Agamma_p�v�Tpopulation�Agamma_c�v�Tcoherence
double omega_12=1.264428211,omega_13=1.264428211,omega_14=59.02343708,omega_23=0,omega_24=57.75900887,omega_34=57.75900887;
//omega_ij���ඥi�M�ඥj�����W�v�t�C
double ReH13,ReH14,ReH23,ReH24;//�U�ඥ�������Ԥ��W�v(�곡)
double ImH13,ImH14,ImH23,ImH24;//�U�ඥ�������Ԥ��W�v(�곡)
double frequency,peroid,FWHM,peak;//frequency:���i���W�v�Cperoid�G�߽ĩP���CFWHM�G�߽ĥb���e�Cpeak�G�Ԥ��W�v�̤j��
const int n=16;//n:�p�ߤ�{���ƥ�
const double pi=3.141592654;//��P�v
const double peroid0=10.87827757;//�w���ۦ쪺�߽ĩP��
int interval,interval2;//interval:���߽İϰ쪺���μ�(��@�P��)�Cinterval2:�S���߽İϰ쪺���μ�(��@�P��)
int numFWHM=5;//interval���Ѽ�
int totalterm;//�i�}���ƥ�
double phase;//�ۦ�


double ReRabi(double &x)//�߽ĥ]���u���(�곡)�A�������*Re[e^{-i*phase}]
{
  double value=0,time=0,factor=0;
  int i=0;
  if(x<0.5*peroid)
    value=exp(-pow(x/FWHM,2));
  else
  {
  i=(x-0.5*peroid)/peroid+1;
  value=exp(-pow((x-i*peroid)/FWHM,2))*cos(-i*phase);
  }
  return peak*value;
}

double ImRabi(double &x)//�߽ĥ]���u���(�곡)�A�������*Im[e^{-i*phase}]
{
  double value=0,time=0,factor=0;
  int i=0;
  if(x<0.5*peroid)
    value=exp(-pow(x/FWHM,2));
  else
  {
  i=(x-0.5*peroid)/peroid+1;
  value=exp(-pow((x-i*peroid)/FWHM,2))*sin(-i*phase);
  }
  return peak*value;
}

void fun(long double k[],long double y[],double x,int i)
{
//�̧Ǭ�rho_11,rho_22,rho_33,rho_44,Im[rho_21],Re[rho_21],Im[rho_31],Re[rho_31],Im[rho_41],Re[rho_41],Im[rho_32],Re[rho_32],
//Im[rho_42],Re[rho_42],Im[rho_43],Re[rho_43]����{��
  k[i*n]=(-ReH13*y[6]-ReH14*y[8]-ImH13*y[7]-ImH14*y[9]-gamma_1*y[0])/i;
  k[i*n+1]=(-ReH23*y[10]-ReH24*y[12]-ImH23*y[11]-ImH24*y[13]-gamma_2*y[1])/i;
  k[i*n+2]=(ReH13*y[6]+ReH23*y[10]+ImH13*y[7]+ImH23*y[11]+gamma_13*y[0]+gamma_23*y[1]-0.5*gamma_p*(y[3]-y[4]))/i;
  k[i*n+3]=(ReH14*y[8]+ReH24*y[12]+ImH14*y[9]+ImH24*y[13]+gamma_14*y[0]+gamma_24*y[1]-0.5*gamma_p*(y[4]-y[3]))/i;
  k[i*n+4]=(-0.5*(gamma_1+gamma_2)*y[4]+omega_12*y[5]+0.5*(ReH23*y[7]-ReH13*y[11]+ReH24*y[9]-ReH14*y[9]-ImH23*y[6]+ImH13*y[10]-ImH24*y[8]+ImH14*y[12]))/i;
  k[i*n+5]=(-0.5*(gamma_1+gamma_2)*y[5]-omega_12*y[4]+0.5*(-ReH23*y[6]-ReH13*y[10]-ReH24*y[8]+ReH14*y[8]-ImH23*y[7]-ImH13*y[11]-ImH24*y[9]-ImH14*y[13]))/i;
  k[i*n+6]=(-0.5*gamma_1*y[6]-(frequency-omega_13)*y[7]+0.5*(ReH13*(y[0]-y[2])+ReH23*y[5]-ReH14*y[15]+ImH23*y[4]+ImH14*y[14]))/i;
  k[i*n+7]=(-0.5*gamma_1*y[7]+(frequency-omega_13)*y[6]+0.5*(ImH13*(y[0]-y[2])-ReH23*y[4]-ReH14*y[14]+ImH23*y[5]-ImH14*y[15]))/i;
  k[i*n+8]=(-0.5*gamma_1*y[8]-(frequency-omega_14)*y[9]+0.5*(ReH14*(y[0]-y[3])+ReH24*y[5]-ReH13*y[15]+ImH24*y[4]-ImH13*y[14]))/i;
  k[i*n+9]=(-0.5*gamma_1*y[9]+(frequency-omega_14)*y[8]+0.5*(ImH14*(y[0]-y[3])-ReH24*y[4]+ReH13*y[14]+ImH24*y[5]-ImH13*y[15]))/i;
  k[i*n+10]=(-0.5*gamma_2*y[10]-(frequency-omega_23)*y[11]+0.5*(ReH23*(y[1]-y[2])+ReH13*y[5]-ReH24*y[15]-ImH13*y[4]+ImH24*y[14]))/i;
  k[i*n+11]=(-0.5*gamma_2*y[11]+(frequency-omega_23)*y[10]+0.5*(ImH23*(y[1]-y[2])+ReH13*y[4]-ReH24*y[14]+ImH13*y[5]-ImH24*y[15]))/i;
  k[i*n+12]=(-0.5*gamma_2*y[12]-(frequency-omega_24)*y[13]+0.5*(ReH24*(y[1]-y[3])+ReH14*y[5]-ReH23*y[15]-ImH14*y[4]-ImH23*y[14]))/i;
  k[i*n+13]=(-0.5*gamma_2*y[13]+(frequency-omega_24)*y[12]+0.5*(ImH24*(y[1]-y[3])+ReH14*y[4]+ReH23*y[14]+ImH14*y[5]-ImH23*y[15]))/i;
  k[i*n+14]=(-0.5*gamma_c*y[14]+(omega_34)*y[15]+0.5*(ReH14*y[7]-ReH13*y[9]+ReH24*y[11]-ReH23*y[13]-ImH14*y[6]+ImH13*y[8]-ImH24*y[10]+ImH23*y[12]))/i;
  k[i*n+15]=(-0.5*gamma_c*y[15]-(omega_34)*y[14]+0.5*(ReH14*y[6]+ReH13*y[8]+ReH24*y[10]+ReH23*y[12]+ImH14*y[7]+ImH13*y[9]+ImH24*y[11]+ImH23*y[13]))/i;
}

void solve(double *result,double &x,int &t,double &dx)
{
  long double factor[totalterm*n],mid[n];
  double next[n]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for(int i=0;i<n;i++)//�x�st_n
    factor[i]=*(result+t*n+i);

  for(int i=1;i<totalterm;i++)//���o�Y��
  {
    for(int j=0;j<n;j++)
      mid[j]=factor[(i-1)*n+j];
    fun(factor,mid,x,i);
  }
  for(int j=0;j<n;j++)//�Dt_n+1
    for(int k=0;k<totalterm;k++)
      next[j]+=factor[k*n+j]*pow(dx,k);

  for(int j=0;j<n;j++)
    *(result+(t+1)*n+j)=next[j];
} 

int main()
{
  fstream file1,file2;//file1�G�����p�⪺�ѼơCfile2:�����ƾ�
  file1.open("input.txt", ios::out | ios::trunc);
  file2.open("data.txt", ios::out | ios::trunc);
  file1.precision(10);  
  file2.precision(10); 

  double endvalue,x0=0,y0[n]={0,0,0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0};//x0:�ɶ��Cy0[n]:�K�ׯx�}�_�l����
  int factor,maxpower,numpluse;
//facter�Ginterval2���ѼơCinterval_peroid:�e�Ϫ��I�ơCmaxpower:�̤j���i�}������ơCnumofpulse:�򥻯߽ļ�

//�U���i�H���ܭp�⪺�Ѽ�  
  frequency=0*2*pi;
  FWHM=0.00085;
  peak=37.6834589;
  gamma_1=0.0052227*2*pi;
  gamma_2=0.0052227*2*pi;
  gamma_13=gamma_14=gamma_1/2.0;
  gamma_23=gamma_24=gamma_2/2.0;
  gamma_p=0.000001*2*pi;
  gamma_c=0.0000005*2*pi;
  maxpower=10;
  numofpulse=1000;
  interval=200;
  factor=50;
  const int nplot=5;//�g�Jfile2������`�ƥءA
  int line[4]={0,1,2,3};//��ƪ���m�Cex:0�Grho_11�A1:rho_22  
  int const cutofnumber=1050;//�פ�߽ļ�
  double phase0=0*pi/180.0;//�ۦ�
//���ܭp�⪺�ѼƵ���

  int maxinterval;
  interval2=peroid*factor;//
  maxinterval=interval2+interval;

  double *presult = new double [(maxinterval+1)*n];
    if(!presult)
      exit(1);

  totalterm=maxpower+1;

  for (int i=0;i<n;i++)
    *(presult+i)=y0[i];
  dx1=numFWHM*FWHM/interval;
  dx2=(peroid-2*numFWHM*FWHM)/interval2;


  double value_0,value_1;

  int k=0;

  for(int j=0;j<interval;j++)
  {
    ReH13=2.61*ReRabi(x0);ReH14=1.302*ReRabi(x0);ReH23=2.269*ReRabi(x0);ReH24=1.906*ReRabi(x0);
    ImH13=2.61*ImRabi(x0);ImH14=1.302*ImRabi(x0);ImH23=2.269*ImRabi(x0);ImH24=1.906*ImRabi(x0);
    solve(presult,x0,j,dx1);
    file2<<setiosflags(ios::left)<<setw(20)<<x0;
    for(int i=0;i<nplot;i++)
      file2<<setiosflags(ios::left)<<setw(20)<<*(presult+j*n+line[i]); 
    file2<<endl;
    x0+=dx1;
    k=k+1;
  }

  dx1=2*numFWHM*FWHM/interval;
  for (int i=0;i<n;i++)
    *(presult+i)=*(presult+interval*n+i);

  while(x0<numpulse*peroid)//���^��p��ܰ򥻯߽ļ�
  {
    for(int j=0;j<interval2;j++)
    {
      ReH13=2.61*ReRabi(x0);ReH14=1.302*ReRabi(x0);ReH23=2.269*ReRabi(x0);ReH24=1.906*ReRabi(x0);
      ImH13=2.61*ImRabi(x0);ImH14=1.302*ImRabi(x0);ImH23=2.269*ImRabi(x0);ImH24=1.906*ImRabi(x0);
      solve(presult,x0,j,dx2);
      file2<<setiosflags(ios::left)<<setw(20)<<x0;
      for(int i=0;i<nplot;i++)
	file2<<setiosflags(ios::left)<<setw(20)<<*(presult+j*n+line[i]); 
      file2<<endl;
      x0+=dx2;
      k=k+1;
    }

    for(int j=interval2;j<maxinterval;j++)
    {
      ReH13=2.61*ReRabi(x0);ReH14=1.302*ReRabi(x0);ReH23=2.269*ReRabi(x0);ReH24=1.906*ReRabi(x0);
      ImH13=2.61*ImRabi(x0);ImH14=1.302*ImRabi(x0);ImH23=2.269*ImRabi(x0);ImH24=1.906*ImRabi(x0);
      solve(presult,x0,j,dx1);
      file2<<setiosflags(ios::left)<<setw(20)<<x0;
      for(int i=0;i<nplot;i++)
	file2<<setiosflags(ios::left)<<setw(20)<<*(presult+j*n+line[i]); 
      file2<<endl;
      x0+=dx1;
      k=k+1;
    }

    for (int i=0;i<n;i++)
      *(presult+i)=*(presult+maxinterval*n+i);
  }    


  int FLAG=1;
  while(FLAG==1)//�p��ܲŦX�������
  {
    value_0=value_1;
    value_1=1;
    for(int count=0;count<10;count++)
    {
    for(int j=0;j<interval2;j++)
    {
      ReH13=2.61*ReRabi(x0);ReH14=1.302*ReRabi(x0);ReH23=2.269*ReRabi(x0);ReH24=1.906*ReRabi(x0);
      ImH13=2.61*ImRabi(x0);ImH14=1.302*ImRabi(x0);ImH23=2.269*ImRabi(x0);ImH24=1.906*ImRabi(x0);
      solve(presult,x0,j,dx2);
      if(value_1>*(presult+j*n))
	value_1=*(presult+j*n);
      file2<<setiosflags(ios::left)<<setw(20)<<x0;
      for(int i=0;i<nplot;i++)
	file2<<setiosflags(ios::left)<<setw(20)<<*(presult+j*n+line[i]); 
      file2<<endl;
      x0+=dx2;
      k=k+1;
    }

    for(int j=interval2;j<maxinterval;j++)
    {
      ReH13=2.61*ReRabi(x0);ReH14=1.302*ReRabi(x0);ReH23=2.269*ReRabi(x0);ReH24=1.906*ReRabi(x0);
      ImH13=2.61*ImRabi(x0);ImH14=1.302*ImRabi(x0);ImH23=2.269*ImRabi(x0);ImH24=1.906*ImRabi(x0);
      solve(presult,x0,j,dx1);
      if(value_1>*(presult+j*n))
	value_1=*(presult+j*n);

      file2<<setiosflags(ios::left)<<setw(20)<<x0;
      for(int i=0;i<nplot;i++)
	file2<<setiosflags(ios::left)<<setw(20)<<*(presult+j*n+line[i]); 
      file2<<endl;
      k+=1;
      x0+=dx1;

    }
    for (int i=0;i<n;i++)
      *(presult+i)=*(presult+maxinterval*n+i);
    }
    if((fabs(value_0-value_1)/value_1)<0.00001)
      FLAG=0;
    else if(x0/peroid>cutofnumber)
    {
      cout<<123<<endl;
      FLAG=0;
    }

  }    
  endvalue=x0;

  file1 <<setiosflags(ios::left)<<setw(15)<<nplot<<setw(15)<<0<<setw(15)<<endvalue<<setw(15)
	<<k<<endl;
  delete [] presult;

  return 0;
} 
