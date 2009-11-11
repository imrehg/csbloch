//CW�p�g�T�ඥ�t��
#include <iostream>
#include <cstdlib>
#include <cctype>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;


//�ඥi�P�ඥj�������Ԥ��W�v
void fun(double ,double ,double ,int );//�p�ߤ�{��
void solve(double*,double &,int &);//�t��k

//�ඥa=�ඥ1�A�ඥb=�ඥ2�A�ඥc=�ඥ3
double gamma_a,gamma_1,gamma_2,gamma_ab,gamma_ac,rabi_ab,rabi_ac,dx,detuning_ab,detuning_ac;
//gamma_i:�ඥ1���I�ܲv gamma_ij:�ඥi�I�ܨ�ඥj���I�ܲv�Crabi_ij:�Ԥ��W�v(���Ҽ{��ӹp�g������t)�Cdetuning_ij�G���նq
double ReH12,ReH13,ImH12,ImH13;//�Ԥ��W�v(�����ҵ��A�Ҽ{��ӹp�g������t)
const int n=9;//�p�ߤ�{���ƥ�
int nplot;//�Ʊ��X����{���ƥ�
int interval;//�ɶ����Ϊ��϶���
int totalterm;//�i�}�����ƥ�


void fun(long double k[],long double y[],double x,int i)
{
  k[i*n]=(-gamma_a*y[0]-2*(ReH12*y[3]+ReH13*y[5]+ImH12*y[4]+ImH13*y[6]))/i;
  k[i*n+1]=(gamma_ab*y[0]-0.5*gamma_1*(y[1]-y[2])+2*(ReH12*y[3]+ImH12*y[4]))/i;
  k[i*n+2]=(gamma_ac*y[0]-0.5*gamma_1*(y[2]-y[1])+2*(ReH13*y[5]+ImH13*y[6]))/i;
  k[i*n+3]=(-0.5*gamma_a*y[3]-detuning_ac*y[4]+ReH12*(y[0]-y[1])-ReH13*y[8]+ImH13*y[7])/i;
  k[i*n+4]=(-0.5*gamma_a*y[4]+detuning_ac*y[3]+ImH12*(y[0]-y[1])-ReH13*y[7]-ImH13*y[8])/i;
  k[i*n+5]=(-0.5*gamma_a*y[5]-detuning_ab*y[6]+ReH13*(y[0]-y[2])-ReH12*y[8]-ImH12*y[7])/i;
  k[i*n+6]=(-0.5*gamma_a*y[6]+detuning_ab*y[5]+ImH13*(y[0]-y[2])+ReH12*y[7]-ImH12*y[8])/i;
  k[i*n+7]=(-0.5*gamma_2*y[7]-(detuning_ac-detuning_ab)*y[8]+ReH13*y[4]-ReH12*y[6]-ImH13*y[3]+ImH12*y[5])/i;
  k[i*n+8]=(-0.5*gamma_2*y[8]+(detuning_ac-detuning_ab)*y[7]+ReH13*y[3]+ReH12*y[5]+ImH13*y[4]+ImH12*y[6])/i;
   
}

void solve(double *result,double &x,int &t)
{

  long double factor[totalterm*n],mid[n];
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
      *(result+(t+1)*n+j)+=factor[k*n+j]*pow(dx,k);

} 

int main()
{
  int maxpower;
  double phase;
  fstream file1,file2;//file1:������J���ѼơCfile2://�����p�⵲�G
  file1.open("input.txt", ios::out | ios::trunc);
  file2.open("data.txt", ios::out | ios::trunc);
  double endvalue,x0=0,y0[n]={0,0.5,0.5,0,0,0,0,0,0};//endcalue:�p��I��ɶ��Cx0:�ɶ��Cy0[n]�G�K�ׯx�}�_�l����
  
  cout<<"input rabifrequency_ab and rabifrequency_ac"<<endl;
  cin>>rabi_ab>>rabi_ac;
  cout<<"input detuning_ab and detuning_ac"<<endl;
  cin>>detuning_ab>>detuning_ac;
  cout<<"input difference phase of the two laser"<<endl;
  cin>>phase;
  cout<<"input gamma_a"<<endl;
  cin>>gamma_a;
  gamma_ab=gamma_ac=0.5*gamma_a;
  cout<<"input gamma_1 and gamma_2"<<endl;
  cin>>gamma_1>>gamma_2;
  cout<<"input maximum power of ADM"<<endl;
  cin>>maxpower;
  cout<<"Input the end value of x"<<endl;
  cin>>endvalue;
  cout<<"Input the initial number of intervals:"<<endl;
  cin>>interval;
  cout<<"input nplot"<<endl;
  cin>>nplot;

  int line[nplot];

  cout<<"input plut line"<<endl;
  for(int i=0;i<nplot;i++)
    cin>>line[i];
  cout<<"!!!  start  !!!  "<<endl;

  file1 <<setiosflags(ios::left)<<setw(10)<<nplot<<setw(10)<<0<<setw(10)<<endvalue<<setw(10)
	<<interval<<endl;
  file1 <<"rabi_ab="<<rabi_ab<<"  rabi_ac="<<rabi_ac<<"  phase="<<phase<<"  gamma_a="<<gamma_a<<endl;
	
  totalterm=maxpower+1;
  double *presult = new double [(interval+1)*n];//�Ҧ��ɶ��I���ƭȦs�󦹫���
    if(!presult)
      exit(1);

  for (int i=0;i<n;i++)
    *(presult+i)=y0[i];
  dx=(endvalue)/interval;
  
  ReH12=rabi_ab;  
  ReH13=rabi_ac*cos(phase);
  ImH12=0;
  ImH13=rabi_ac*sin(phase);    
  for(int k=0;k<interval;k++)
  { 
    solve(presult,x0,k);
    file2<<setiosflags(ios::left)<<setw(15)<<x0;
    for(int i=0;i<nplot;i++)
      file2<<setiosflags(ios::left)<<setw(15)<<*(presult+k*n+line[i]); 
    file2<<endl;
    x0+=dx;
  }

  delete [] presult;

  return 0;
} 
