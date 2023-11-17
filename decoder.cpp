#include <iostream>
#include <iterator>
#include <random>
#include <cmath>
#include <ctime>

using namespace std;
const double E=1.75;
const int N=40;
const int STATE=8;
const double INF=1e5;

//добавляет белый шум
void addNoise(double *r, const int n);
//отображение {0,1}->{-E,+E}
void mapping(const int* v, double *r, const int n);
//кодирование по схеме
void coder(const int *u, int *v);

void MaxLogMap(const double *r, int *u);
double get_v(const int j, const int m);
int get_backward_state(const int j, const int m);
int get_forward_state(const int j, const int m);
//Выводить вектор с целыми числами
void print(const int *arrays,const int n,const char tap=' ');
//Выводить вектор с вещественными числами
void print(const double *arrays, const int n,const char tap=' ');
//Преобразует строку в массив с размером N
int strToArrays(const string x_str, int*x,const int n);
//Количество неправильно декадированных битов
const int error(const int *x, const int *u);
const int error(const double *r, const int *u);
int sgn(const double a);
int main(){
    //иницилизация
    int x[N],y[N+3],v[N+3], mistake=0,u[N],x_len, iter;
    double r[N+3];
    unsigned int time;//время работы некоторого блока кода
    
    //считывание строки
    string x_str;
    cout<<endl<<"  Исходное сообщение X"<<endl;
    cin>>x_str;
    
    //преобразование исходной строки в массив
    x_len=x_str.length();
    for(int i=0; i<N-x_len;++i)
        x_str='0'+x_str;
    if(strToArrays(x_str,x,N)==0) return 0;
    
    //
    cout<<endl<<"  Переданное сообщение X"<<endl;
    print(x,N,'\0');
    //кодирование по схеме
    coder(x,y);
    cout<<endl<<"  Закодированное сообщение y"<<endl;
    print(y,N+3,'\0');
    //отображание {0,1} -> {-E,E}, где E - положительное число
    mapping(y,r,N+3);
    //добавление шума АБГШ
    addNoise(r,N+3);
    cout<<endl<<"  Закодированное сообщение r c шумом"<<endl;
    print(r,N+3);

    //декодирование с мягким решением
    time=clock();
    MaxLogMap(r,u);
    time=clock()-time;
    cout<<endl<<"  Декодированное сообщение u"<<endl;
    print(u,N,'\0');
    cout<<endl<<"  Закодируем декодированное сообщение u и обозначим v"<<endl;
    coder(u,v);
    print(v,N+3,'\0');
    
    cout<<endl<<"  Количество ошибок y от r "<<error(r,y)<<"/43"<<endl;
    cout<<endl<<"  Количество ошибок v от r "<<error(r,v)<<"/43"<<endl;

    cout<<endl<<"  Количество ошибок u от x "<<error(x,u)<<"/40     Время работы "<<double(time)/CLOCKS_PER_SEC <<" с";
    cout<<endl;
    
    return 0;
}


void print(const int *arrays,const int n,const char tap){
    for(int i=0;i<n;++i)
        cout<<arrays[i]<<tap;
    cout<<endl;
}
void print(const double *arrays, const int n,const char tap){
    for(int i=0;i<n;++i)
        cout<<arrays[i]<<tap;
    cout<<endl;
}
void addNoise(double *r, const int n){
    const double mean = 0.0;
    const double stddev = 1.;
    default_random_engine generator;
    normal_distribution<double> dist(mean, stddev);
    for(int i=0;i<n;++i)
        r[i]=r[i]+dist(generator);
}
void mapping(const int* v, double *r, const int n){
    for(int i=0; i<n;++i)
        if(v[i]==0) r[i]=-E;
        else r[i]=E;
}
void coder(const int *u, int *v){
    bool d1=0,d2=0,d3=0,t=0;
    for(int i=0;i<N; ++i){
        v[i]=(u[i]+d1+d2)%2;
        t=d3; d3=d2; d2=d1;
        d1=(t+d3+u[i])%2;
    }
    for(int i=0; i<3; ++i){
        v[N+i]=(d1+d3)%2;
        d3=d2;d2=d1;d1=0;
    }
    
}

double get_v(const int j, const int m){
    int d1=m>>2, d2=(m>>1)&1, d3=m&1;
    return (((j+d1+d2)%2)*2-1)*E;
}

int get_backward_state(const int j, const int m){
    //k<N
    int d1=m>>2, d2=(m>>1)&1, d3=m&1;
    return d2*4+d3*2+(d1+d3+j)%2;
}
int get_forward_state(const int j, const int m){
    int d1=m>>2, d2=(m>>1)&1, d3=m&1;
    //int d1=m/4, d2=(m%4)/2, d3=m%2;
    int t=(d2+d3+j)%2;
    d3=d2;d2=d1; d1=t;
    return d1*4+d2*2+d3;
}
void MaxLogMap(const double *r, int *u){
    //log(gamma_0)    log(gamma_1)
    double gamma_0[STATE*N], gamma_1[STATE*N];
    //log(alpha)       log(beta)
    double alpha[STATE*N], beta[STATE*N];
    
    //Инициализация alpha
    //состояние начинается с 000
    alpha[0]=0;// т.к. alpha[0]=1, то log(alpha[0])=1
    for(int m=1; m<STATE;++m) alpha[m]=-INF; // т.к. alpha[m]=0, то log(alpha[0])=-INF
    

    //Инициализация beta
    //состояние после N шага может быть любым
    if(r[N]==0 or r[N+1]==0 or r[N+2]==0)
        for(int m=0; m<STATE;++m)
            beta[m+(N-1)*STATE]=log(1./STATE);
    else{
        for(int m=0; m<STATE;++m)
            beta[m+(N-1)*STATE]=-INF;
        int d2=(r[N+1]>0?1:0), d1=(r[N+2]>0?1:0), d3=((r[N]>0?1:0)+d1)%2;
        beta[d1*4+d2*2+d3+(N-1)*STATE]=0;
    }

    
    //Расчет log(gamma)
    for(int k=0; k<N;++k){
        for(int m=0; m<STATE; ++m){
            gamma_0[m+k*STATE]=-pow(r[k]-get_v(0,m),2)/2;//по формуле, где стандартное отклонение=1, а свободные коэфф сокращаются
            gamma_1[m+k*STATE]=-pow(r[k]-get_v(1,m),2)/2;
            //cout<<get_v(0,m)<<" ";
        }
    }

    //Расчет log(alpha)
    for(int k=1;k<N;++k){
        for(int m=0;m<STATE;++m){
            int b0=get_backward_state(0,m),b1=get_backward_state(1,m);
            double alp0=alpha[b0+(k-1)*STATE]+gamma_0[b0+(k-1)*STATE];
            double alp1=alpha[b1+(k-1)*STATE]+gamma_1[b1+(k-1)*STATE];
            alpha[m+k*STATE]=max(alp0,alp1);
        }
    }
    cout<<endl<<"alpha "<<endl;
    for(int k=0;k<N;++k){
        for(int m=0;m<STATE;++m)
            cout<<alpha[m+k*STATE]<<' ';
        cout<<endl;
    }

    
    //Расчет log(beta)
    for(int k=N-2;k>=0;--k){
        for(int m=0;m<STATE;++m){
            int f0=get_forward_state(0,m),f1=get_forward_state(1,m);
            double beta0=gamma_0[m+(k+1)*STATE]+beta[f0+(k+1)*STATE];
            double beta1=gamma_1[m+(k+1)*STATE]+beta[f1+(k+1)*STATE];
            beta[m+k*STATE]=max(beta0,beta1);
        }
    }
    cout<<endl<<"beta "<<endl;
    for(int k=0;k<N;++k){
        for(int m=0;m<STATE;++m)
            cout<<beta[m+k*STATE]<<' ';
        cout<<endl;
    }
    
    for(int k=0;k<N;++k){
        double max_1=-INF,max_0=-INF;
        for(int m=0; m<STATE;++m){
            int f0=get_forward_state(0,m),f1=get_forward_state(1,m);
            double value1=alpha[m+k*STATE]+gamma_1[m+k*STATE]+beta[f1+k*STATE];
            double value0=alpha[m+k*STATE]+gamma_0[m+k*STATE]+beta[f0+k*STATE];
            max_1=max(max_1,value1);
            max_0=max(max_0,value0);
        }
        
        if(max_1>max_0) u[k]=1;
        else u[k]=0;
    }
}


int strToArrays(const string x_str, int*x,const int n){
    for(int i=0;i<n;++i){
        if(x_str[i]=='0')
            x[i]=0;
        else if(x_str[i]=='1')
            x[i]=1;
        else{
            cout<<"Некорректный ввод последовательности"<<endl;
            return 0;
        }
    }
    return 1;
}
const int error(const int *x, const int *u){
    int loss=0;
    for(int i=0; i<N;++i)
        if(x[i]!=u[i]) loss++;
    return loss;
}
const int error(const double *r, const int *u){
    int loss=0;
    for(int i=0; i<N;++i)
        
        if(sgn(r[i])!=(2*u[i]-1)) loss++;
    return loss;
}

int sgn(const double a){
    if(a>0) return 1;
    if(a<0) return -1;
    return 0;
}
