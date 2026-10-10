"""Independent rational check of the frozen twelve-fibre profile.

Only Python's standard library is used. No floating-point value controls
acceptance. The candidate coefficients came from the companion BVP scout;
the trig enclosures, polynomial integration, and residual check are separate.
This proves a method obstruction, NOT a zero-proportion lower bound.
"""
from fractions import Fraction as F
from math import factorial

if not __debug__:
    raise RuntimeError('Proof checks require assertions; do not run with -O or -OO.')

DEN = 10**14
HALF = [
    [81452785069661,7292373590683,-495291728588,-8832261778,468777856,3579876,-172336],
    [94237762244332,5325472480437,-540048842503,-5392125542,498965786,1817587,-181576],
    [102829051521759,3211889787491,-561973129360,-1752444969,508473653,-36353,-183294],
    [106611512235993,1220224746606,-555948929227,1548835301,491383198,-1705511,-175034],
    [105872959612158,3026144941748,-520652451189,-4834131202,443772398,2367437,-154759],
    [110067398521096,987168788008,-542216720202,-1573291099,463005280,777705,-161751],
]

def add(a,b):
    return [(a[k] if k<len(a) else F(0))+(b[k] if k<len(b) else F(0))
            for k in range(max(len(a),len(b)))]

def scale(a,t):
    return [t*x for x in a]

def product(a,b):
    r=[F(0)]*(len(a)+len(b)-1)
    for j,x in enumerate(a):
        for k,y in enumerate(b): r[j+k]+=x*y
    return r

def primitive(a):
    return [F(0)]+[x/F(k+1) for k,x in enumerate(a)]

def at(a,x):
    v=F(0)
    for t in reversed(a): v=v*x+t
    return v

def integral(a):
    # Each physical fibre has length 1/12 and y runs from -1 to 1.
    return sum((2*x/F(k+1) for k,x in enumerate(a) if k%2==0),F(0))/24

def atan_bounds(x,n=32):
    a=sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(n)),F(0))
    b=a+(-1)**n*x**(2*n+1)/F(2*n+1)
    return min(a,b),max(a,b)

def outward(interval, digits=26):
    s=10**digits
    lo,hi=interval
    return F((lo*s).__floor__(),s),F((hi*s).__ceil__(),s)

a=atan_bounds(F(1,5)); b=atan_bounds(F(1,239))
# Machin's identity pi=16 atan(1/5)-4 atan(1/239).
PI=outward((16*a[0]-4*b[1],16*a[1]-4*b[0]),36)
assert F(314,100)<PI[0]<PI[1]<F(22,7)

def cos_pi(r):
    # Cosine period and parity put the rational multiplier in [0,1].
    r=r%2
    if r>1: r=2-r
    x=r*(PI[0]+PI[1])/2
    value=sum(((-1)**k*x**(2*k)/F(factorial(2*k)) for k in range(31)),F(0))
    # Degree-61 Taylor formula (odd coefficient zero) plus |cos'|<=1.
    error=x**62/F(factorial(62))+r*(PI[1]-PI[0])/2
    return outward((value-error,value+error))

def iadd(a,b): return a[0]+b[0],a[1]+b[1]
def imul(a,b):
    v=[x*y for x in a for y in b]
    return min(v),max(v)
def iscale(a,t): return imul(a,(t,t))
def point(x): return x,x
def isub(a,b): return a[0]-b[1],a[1]-b[0]

def matrix(q):
    lam=9-6/q
    alpha=F(5189,5000); gamma=F(10289,10000)
    c0=(1-lam)*F(2,3)+lam*F(3,4)
    shifts={}
    for j in (1,2):
        shifts[4*j]=iscale(iadd(point(F(2)),iscale(cos_pi(2*j*alpha/3),2)),(1-lam)/3)
    for j in (1,2,3):
        terms=iadd(point(F(3)),iadd(iscale(cos_pi(j*gamma/2),4),iscale(cos_pi(j*gamma),2)))
        shifts[3*j]=iscale(terms,lam/6)
    return [[point(1-q*c0) if j==k else iscale(shifts[abs(j-k)],-q/2)
             if abs(j-k) in shifts else point(F(0)) for k in range(12)] for j in range(12)]

p=[[F(x,DEN) for x in row] for row in HALF]
p += [[(-1)**k*x for k,x in enumerate(row)] for row in reversed(p)]
normalizer=sum((integral(row) for row in p),F(0))
assert normalizer>0
p=[scale(row,1/normalizer) for row in p]
assert sum((integral(row) for row in p),F(0))==1
assert all(p[11-j]==[(-1)**k*x for k,x in enumerate(p[j])] for j in range(12))
positivity=min(row[0]-sum(map(abs,row[1:])) for row in p)
assert positivity>F(7,10)

# Exact integral potential: g(-1/2)=1/2 by evenness and normalization,
# Fcdf'=f, g'=2Fcdf-1. Integrate in y with dx=dy/24.
potentials=[]; cdfs=[]; mass=F(0); g_start=F(1,2)
for row in p:
    cdf=scale(primitive(row),F(1,24))
    cdf[0]+=mass-at(cdf,-1)
    derivative=scale(cdf,F(1,12)); derivative[0]-=F(1,24)
    g=primitive(derivative); g[0]+=g_start-at(g,-1)
    potentials.append(g); cdfs.append(cdf)
    mass=at(cdf,1); g_start=at(g,1)
assert mass==1 and g_start==F(1,2)
J=sum((integral(product(row,g)) for row,g in zip(p,potentials)),F(0))
J_check=2*sum((integral(add(cdf,scale(product(cdf,cdf),-1))) for cdf in cdfs),F(0))
assert J==J_check
inner=[[integral(product(a,b)) for b in p] for a in p]

def audit(q):
    mat=matrix(q)
    diagonal_lower=min(mat[j][j][0]-sum(max(abs(mat[j][k][0]),abs(mat[j][k][1]))
                       for k in range(12) if k!=j) for j in range(12))
    kappa=F(1,10)
    assert diagonal_lower-2/PI[0]**2>kappa
    energy=point(J)
    for j in range(12):
        for k in range(12): energy=iadd(energy,iscale(mat[j][k],inner[j][k]))
    midpoint=[[(a+b)/2 for a,b in row] for row in mat]
    c=(energy[0]+energy[1])/2
    r2=F(0); eps2=F(0)
    for j in range(12):
        r=list(potentials[j])
        for k in range(12): r=add(r,scale(p[k],midpoint[j][k]))
        r[0]-=c
        r2+=integral(product(r,r))
        eps=sum(((mat[j][k][1]-mat[j][k][0])/2*sum(map(abs,p[k]))
                 for k in range(12)),F(0))
        eps2+=eps**2/12
    assert r2>=0
    residual_upper=2*r2+2*eps2
    minimum_lower=energy[0]-residual_upper/kappa
    target=2-2*q
    return {'q':q,'matrix_lower':diagonal_lower,'energy':energy,
            'minimum_lower':minimum_lower,'target':target,'residual_upper':residual_upper}

old=audit(F(337017,500000))
new=audit(F(16851,25000))
assert old['energy'][1]<old['target']-F(5,10**6)
assert new['minimum_lower']>new['target']+F(18,10**8)
assert new['residual_upper']<F(1,10**18)
assert new['energy'][1]<new['target']+F(181,10**9)
print('PASS: independent exact-rational twelve-fibre checks.')
print('Profile is exactly normalized, even, and pointwise > 7/10.')
print('Coercivity is > 1/10 on every zero-mean real L2 perturbation.')
print('At q=337017/500000: positive profile disproves extending this periodic test unchanged; margin > 5/10^6.')
print('At q=16851/25000: all-profile variational minimum exceeds 2-2q by > 18/10^8.')
print('The minimum is less than 2-2q+181/10^9; residual norm squared < 1/10^18.')
print('This certifies a 67.404% joint-affine METHOD CEILING, not a zero record.')
