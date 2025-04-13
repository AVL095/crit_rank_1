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

/////////////////////////////////////////////////////////////////

void ansari(int n1,int n2,double *p,int &idimp) {

    double *x,*y,*z;
    int  n, nm, k, i, j, kb,idimx, idimy,idimz, ier;
    long int iw;
   
    if(n1<=0) {p[1]=1;idimp=1;return;}
    if(n1>=n2) {p[1]=1;idimp=1;return;}
    n=n1+n2;
    nm = (n1 < n2) ? n1 : n2;
    iw=1000;
    x=new double[iw];
    y=new double[iw];
    z=new double[iw];
    for (i=0;i<=iw;i++) {x[i]=0;y[i]=0;z[i]=0;}  


    k=n2+1;x[1]=1;y[1]=1;x[k]=-1;idimx=k;
    for (i=1;i<=n1-1;i++) {
        k--; y[k]=-1;
        pmpy(z,idimz,x,idimx,y,k);  //умножение полиномов
        for (j=1;j<=idimz;j++) x[j]=z[j];
        idimx=idimz;
        for (j=1;j<=n;j++) y[j]=0;
        y[1]=1;
    }
    for (i=0;i<=iw;i++) y[i]=0;
    y[1]=1;kb=nm+1;
    for (i=1;i<=nm; i++) {
        y[kb]=-1;
        pdiv(p,idimp,x,idimx,y,kb,0,ier); //деление полиномов
        kb--;
        for (j= 1;j<=idimp;j++) {
            x[j]=p[j];
        }
        idimx=idimp;
        for (j=1;j<=n;j++) y[j] = 0;
        y[1]=1;
    }

    delete[] x,y,z;
 }
    
    
 ////////////////////////////////////////////////////////////////

int main() {

  int i,jj,m,n,N,k,ier,a0;
  long int j,iw;
  string st,ff;
  double s,s1;
  vector <double> pw;
  vector <int> wrange;
 
  ff="Ansari";
  ifstream inp("Inp/" + ff + ".inp");
  ofstream out("Out/" + ff + ".out");

  inp>>st;
  inp>>m; //Sample 1
  inp>>st;
  inp>>n;  //Sample 2
  inp.close();

  out<<"Criterion:"<<ff<<endl;
  out<<"m="<<m<<endl;
  out<<"n="<<n<<endl;
  N=m+n;

  iw=1000;
  int idimx1,idimx2,idimx3,idimx4,idimx5,idimx6,idimx;
  double *x,*x1,*x2,*x3,*x4,*x5,*x6;
  x=new double[iw];
  x1=new double[iw];
  x2=new double[iw];
  x3=new double[iw];
  x4=new double[iw];
  x5=new double[iw];
  x6=new double[iw];
  
/////////////////////////////////////////////////////////

for(i=0;i<=iw;i++) {x1[i]=0;x2[i]=0;x3[i]=0;x4[i];x5[i]=0;x6[i]=0;}
  
if(N%2==0) {

for(i=0;i<=m;i++) {
    for(j=0;j<=iw;j++) x[j]=0;
    idimx=i*(i+1)/2+(m-i)*(m-i+1)/2+1; 
    x[idimx]=1;
    ansari(i,N/2,x1,idimx1);
    ansari(m-i,N/2,x2,idimx2);
    pmpy(x3,idimx3,x1,idimx1,x2,idimx2);
    pmpy(x4,idimx4,x3,idimx3,x,idimx);
    for(j=1;j<=idimx4;j++) x5[j]=x5[j]+x4[j];
} //end for
   k=0;
   for(j=1;j<=2*idimx4;j++)  {
        if(x5[j] !=0) {
         k++;
         x6[k]=x5[j];
        }
   } 
 }
 else {
  for(jj=0;jj<2;jj++) {
    for(i=0;i<=iw;i++) x1[i]=0;
   for(i=0;i<m-jj+1;i++) {
      for(j=0;j<=iw;j++) x[j]=0;
      idimx=i*(i+1)/2+(m-jj-i)*(m-jj-i+1)/2+1; 
      x[idimx]=1;
    ansari(i,(N-1)/2,x1,idimx1);
    ansari(m-jj-i,(N-1)/2,x2,idimx2);
    pmpy(x3,idimx3,x1,idimx1,x2,idimx2);
    pmpy(x4,idimx4,x3,idimx3,x,idimx);
    for(j=1;j<=idimx4;j++) x5[j]=x5[j]+x4[j];
  } //end for i
    k=0;
    for(j=1;j<=2*idimx4;j++)  {
        if(x5[j] !=0) {
         k++;x6[k]=x5[j];
        }
    } 
  } //end for jj 
 } //end if

  s=0;
  for(i=1;i<=k;i++) {x5[i]=x6[k-i+1];s+=x5[i];}
  s1=0;
  a0=int((m+1)/2)*(1+int(m/2));
  s=tgamma(N+1)/(tgamma(m+1)*tgamma(n+1));
   for(i=1;i<=k;i++) {
     s1+=x5[i]/s;
     x6[i]=s1;
    }
    for (i=0;i<k;i++)  {
        pw.push_back(x6[i+1]);
        wrange.push_back(a0+i);
    }
    out << "Var="<<s<<endl;
    out << "Size="<<k<<endl;
    for (i=0;i<k;i++) out<<i+1<<". "<<setprecision(0)<<fixed<<x5[i+1]<<"  "<<setprecision(0)<<fixed<<wrange[i]<<"  "<<setprecision(15)<<fixed<<pw[i]<<endl;
    out <<endl;

  pw.clear();
  wrange.clear(); 
  out.close();
  delete [] x,x1,x2,x3,x4,x5,x6;
  return 0;
}
