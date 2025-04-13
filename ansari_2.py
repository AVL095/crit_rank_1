from sympy import symbols, expand, simplify, binomial, product,Rational
from sympy import degree

N,m,q,k = symbols('N m q k ')

def c(i):
    return q**(Rational(1, 2)*i*(i + 1)) 

def frs(s,r):
    if s==0: return 1
    t1=product(1-q**k,(k,s+1,r))
    t2=product(1-q**k,(k,1,r-s))
    return(expand(simplify(t1/t2)))

def ansari(N, m):
    s = 0
    binomial_coefficient=binomial(N, m)
    if (N%2==0):
         q_term =q**(-Rational(m, 2))
         for i in range(m+1):s=s+c(i)*c(m - i)*frs(i,N//2)*frs(m-i,N//2)
         result =expand(simplify(s*q_term))
    else:
         for j in range(2):
              for i in range(m-j+1):  s=s+c(i)*c(m -j- i)*frs(i,(N-1)//2)*frs(m-j-i,(N-1)//2)
         result =expand(simplify(s))
    return result/binomial_coefficient

#########################################################

m=4
n=4
N=m+n
z =ansari(N, m)

fout=open("ansari_2.out",'w')
z=z.as_ordered_terms() 
k=len(z)
wrange=[]
pw=[]
w=[]
for i in range(k): wrange.append(degree(z[k-i-1])+4)
for i in range(k):w.append((z[i].subs([(q,1)])))

sum=0
for i in range(k): sum=sum+w[i]
s=0
for i in range(k):
    s=s+float(w[i])
    pw.append(s*1.0/sum)

print("Ansari-Bradley",file=fout)
print("Size",";",len(w),file=fout)
print("m",";",m,";","n",";",n,file=fout)
print("Freq",";","Rang",";","Prob",";",file=fout)

for i in range(k):print(w[i],";",wrange[i],";",pw[i],file=fout)

fout.close()
