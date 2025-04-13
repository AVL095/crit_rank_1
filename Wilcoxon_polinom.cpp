#include <string.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

void pnorm(double *x, int idimx, double eps);
void pmpy(double * z, int& idimz,double * x, int idimx,double * y, int idimy);
void pdiv(double * p, int& idimp, double * x, int& idimx, double * y, int &idimy, double eps, int& ier);
void wilcoxon(int n1,int n2,vector <double> &pw,vector <double> &wrange);


/////////////////////////////////////////////////////////////////////////////////////////////

void pnorm(double *x, int idimx, double eps) {
    double z;
    if (idimx <= 0)  return;
    if (idimx > 0) {
        z = abs(x[idimx]);
        z -= eps;
        if (z <= 0) {
            idimx = idimx - 1;
            pnorm(x, idimx, eps);
        }
    }
}

/////////////////Умножение полиномов/////////////////////////////////////////////////////

void pmpy(double * z, int& idimz, double * x, int idimx, double * y, int idimy) {
    double z1, z2;
    int i, j, k;
    z1 = double(idimx)*double(idimy);
    z2 = double(idimx)*double(idimx);
    if (z1 <= 0) {
        idimz = 0;
        return;
    }
    if(z2<=0) {
      idimz=0;
      return;
    }
   idimz = idimx + idimy - 1;
    for(i=0;i<=idimz;i++) z[i]=0; 
    for (i = 1; i <= idimx; i++) {
        for (j = 1; j <= idimy; j++) {
            k = i + j - 1;
            z[k] = x[i] * y[j] + z[k];
        }
    }
}

////////////////////Деление полиномов//////////////////////////////////////////////////////////////

void pdiv(double * p, int& idimp, double * x, int& idimx, double * y, int& idimy, double tol, int& ier) {

    int i, j, ii, k;

    pnorm(y, idimy, tol);
    if (idimy <= 0) {
        ier = 1;return;
    }
        idimp = idimx - idimy + 1;
        if (idimp == 0) {
            ier = 0;return;
        }
        if (idimp < 0) {
            idimp = 0; ier=0;return;
        }
        idimx=idimy-1;i=idimp;
        while(i>0) { 
           ii = i + idimx;
           p[i]= x[ii]/y[idimy];
           for (k=1; k<=idimx; k++) {
             j=k-1+i;
             x[j]=x[j]-p[i]*y[k];
           }
           i--;
        }
        pnorm(x, idimx, tol);ier=0;
}

///////Точное распределение критерия Уилкоксона с помощью производящей функции моментов (умножение и деление полиномов)////////////////////

void wilcoxon(int n1,int n2,vector <double> &pw,vector <double> &wrange) {

    double *x,*y,*z,*p,s;
    int  n, nm, k, i, j, kb,idimx, idimy, idimp, idimz, ier;
    long int iw;
   
    n = n1 + n2;
    nm = (n1 < n2) ? n1 : n2;
    iw=1000;
    x=new double[iw];
    y=new double[iw];
    z=new double[iw];
    p=new double[iw];
  
    k=n+1; x[1]=-1;y[1]=-1;x[k]=1;idimx = k;
 
    for (i = 1; i <= nm - 1; i++) {
        k = k - 1; y[k] = 1;
        pmpy(z, idimz, x, idimx, y, k);  //умножение полиномов
        for (j = 1; j <= idimz; j++)   x[j] = z[j];
        idimx = idimz;
        for (j = 2; j <= n; j++)  y[j] = 0;
    }
    y[1] = -1; kb = nm + 1;
    for (i = 1; i <= nm; i++) {
        y[kb] = 1;
        pdiv(p, idimp, x, idimx, y, kb, 0, ier); //деление полиномов
        kb = kb - 1;
        for (j = 1; j <= idimp; j++) {
            x[j] = p[j];
        }
        idimx = idimp;
        for (j = 2; j <= n; j++)  y[j] = 0;
    }

     ier=1;
     if (idimx==0) ier=0;
     s=0;
     for(i=1;i<=idimp;i++) {
        s+=p[i];
        y[i]=(i-1)+nm*(nm+1)/2; //Wilcoxon
        //y[i]=double(i-1); //Mann-Whitney
        p[i]=s;
    }
    for (i =0; i < idimp; i++) {
          pw.push_back(p[i+1]/s);
          wrange.push_back(y[i+1]);
    }
    delete[] x,y,z,p;
 }
    
 ////////////////////////////////////////////////////////////////

int main() {

  int i,k,*m;
  long int j,nn;
  string st,ff;
  vector <double> pw;
  vector <double> wrange;
 
  ff="Wilcoxon";
  ifstream inp("Inp/" + ff + ".inp");
  ofstream out("Out/" + ff + ".out");

  inp>>st;
  inp>>k;
  m=new int[k];
  inp>>st;
  for(i=0;i<k;i++) inp>>m[i];
  inp.close();

  out<<"Criterion:"<<ff<<endl;
  for(i=0;i<k;i++) out<<m[i]<<";";
  out <<endl;

   wilcoxon(m[0],m[1],pw,wrange);  

  nn=pw.size();
  out << "Size=" << nn << endl;
  for (j=0;j<pw.size();j++) out<<(j+1)<<". "<<setprecision(0)<<fixed<<wrange[j]<<"  "<<setprecision(20)<<fixed<<pw[j]<<endl;
  out <<endl;
  out.close();

  pw.clear();
  wrange.clear(); 
  delete [] m;
  return 0;
}
