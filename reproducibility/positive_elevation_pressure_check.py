"""Exact, continuum-valid signed-Taylor pressure refinement; no numeric search."""
from fractions import Fraction as F
import json


def require(condition,label):
    if not condition:
        raise ArithmeticError(label)


def check():
    coeff=[F(v,10**6) for v in (122290,-45012,16241,1887,2711,3731,-5307,2982,-5153,5358)]
    def derivative(k):
        # parity works exactly for negative integer n as well.
        def sd(n):
            return F(1 if n%2==0 else -1,n) if n else F(0)
        return sd(k)+sum((v*(sd(k-j)+sd(k+j))/2
                          for j,v in enumerate(coeff,1)),F(0))
    p,c,C=F(17,4000),F(1891,200000),F(19,4)
    require(sum(map(abs,coeff))<1,'Positive profile, integral one')
    require((1+sum(map(abs,coeff)))*F(22,7)**2/3<4,'Global second derivative bound')
    require(derivative(1)==F(-260321359,264000000),'Signed derivative at one')
    require(derivative(2)==F(2577793683,6160000000),'Signed derivative at two')
    require(F(5,2)*p>c,'All nearest integers at least three paid')
    proofs=[]
    for k,r in ((1,F(1,20)),(2,F(1,25))):
        signed_value=coeff[k-1]/2
        A,B=abs(signed_value),abs(derivative(k))
        require(signed_value*derivative(k)<0,'Kernel magnitude decreases to the right')
        require(A-B*r-2*r*r>0,'Right Taylor lower bound positive on inner interval')
        q0=2*A*A+p*k-c
        q1=-4*A*B+p
        q2=2*B*B-8*A+C
        margin=4*q0*q2-q1*q1
        require(q0>0 and q2>0 and margin>0,'Right quadratic is positive on real line')
        left_slope=4*A*(B-2*r)-p
        require(B-2*r>0 and left_slope>0,'Left Taylor and pressure lower bound increases')
        require(2*C*r-p>0,'Outer pressure derivative positive')
        outer=p*(k-r)+C*r*r-c
        require(outer>0,'Outer interval starts positive')
        proofs.append(dict(integer=k,cutoff=str(r),A=str(A),B=str(B),
                           right_discriminant_margin=str(margin),left_slope=str(left_slope),
                           outer_margin=str(outer)))
    a,b=F(8722,10**6),F(4389,10**6)
    M=F(350)
    Q=1+M*C
    na,nb=a+(c-a)/Q,b+(p-b)/Q
    beta=F(13707941,2**28)
    require(0<na<beta and 0<nb<b,'Pair and span guards')
    B0=7*a+F(3882,10**6)
    Bstar=(M*C*B0+c+6*C)/Q
    smooth=(M*C*270+6)/Q
    print(json.dumps(dict(status='PASS: exact signed-Taylor pressure refinement',
          p=str(p),c=str(c),C=str(C),proofs=proofs,
          rigidity=str(M),bootstrap_denominator=str(Q),a=str(na),b=str(nb),
          affine_boundary=str(Bstar),smoothing_coefficient=str(smooth),
          acceptance='Every acceptance comparison used exact Fraction arithmetic'),indent=2))


if __name__=='__main__':
    check()
