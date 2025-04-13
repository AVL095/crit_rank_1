function pnorm(x,idimx,eps) {
    let z;

    if (idimx <= 0)  return;
    if (idimx > 0) {
        z = Math.abs(x[idimx]);
        z-=eps;
        if (z <= 0) {
            idimx = idimx - 1;
            pnorm(x, idimx, eps);
        }
    }
}

/////////////////Умножение полиномов/////////////////////////////////////////////////////

function pmpy(z,x,idimx,y,idimy) {
    let z1, z2,idimz;
    let i, j, k;

    z1 =idimx*idimy;
    z2 = idimx*idimx;
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
  return idimz;
}

////////////////////Деление полиномов//////////////////////////////////////////////////////////////

function pdiv(p,x,idimx,y,idimy,tol,ier) {

    let i, j, ii, k,idimp;

    pnorm(y, idimy, tol);
    if (idimy <= 0) {
        ier = 1;return;
    }
        idimp = idimx - idimy + 1;
        if (idimp === 0) {
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
        pnorm(x, idimx, tol);
        ier=0;
        return idimp;
}

///////Точное распределение критерия Уилкоксона с помощью производящей функции моментов (умножение и деление полиномов)////////////////////

function wilcoxon(n1,n2,pw,wrange) {

    let x=[],y=[],z=[],p=[],s;
    let  n, nm, k, i, j, kb,idimx, idimy, idimp, idimz, ier;
   
    n = n1 + n2;
    nm=(n1 < n2) ? n1 : n2;
  
     for(i=0;i<=n+1;i++) {
       x[i]=0;y[i]=0;p[i]=0;z[i]=0;
    }


    k=n+1; x[1]=-1;y[1]=-1;x[k]=1;idimx = k;
 
    for (i = 1; i <= nm - 1; i++) {
        k = k - 1; y[k] = 1;
        idimz=pmpy(z,x,idimx, y, k);  //умножение полиномов
        for (j = 1; j <= idimz; j++)   x[j] = z[j];
        idimx = idimz;
        for (j = 2; j <= n; j++)  y[j] = 0;
    }
    y[1] = -1; kb = nm + 1;
    for (i = 1; i <= nm; i++) {
        y[kb] = 1;
        idimp=pdiv(p,x, idimx, y, kb, 0, ier); //деление полиномов
        kb = kb - 1;
        for (j = 1; j <= idimp; j++) {
            x[j] = p[j];
        }
        idimx = idimp;
        for (j = 2; j <= n; j++)  y[j] = 0;
    }

     ier=1;
     if (idimx===0) ier=0;
     s=0;
     for(i=1;i<=idimp;i++) {
        s+=p[i];
        y[i]=parseFloat(i-1+nm*(nm+1)*0.5); //Wilcoxon
        //y[i]=double(i-1); //Mann-Whitney
        p[i]=parseFloat(s);
    }

    for (i=0;i<idimp; i++) {
          pw[i]=p[i+1]/s;
          wrange[i]=y[i+1];
    }
 }
    
 ////////////////////////////////////////////////////////////////

function Test() {

  let i,k,nn,ff;
  let pw=[],wrange=[];
 
  ff="Wilcoxon";
  m=4;
  n=5; 

  document.write("Criterion:",ff,"<br>");
  document.write("m=",m," n=",n,"<br>");

  wilcoxon(m,n,pw,wrange);  

  nn=pw.length;
  document.write("Size=",nn,"<br>");

  for (i=0;i<nn;i++) document.write(wrange[i],"   ",pw[i],"<br>");

}
