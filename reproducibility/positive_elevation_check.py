"""Exact phase, tail and aggregate frame-payment certificate; standard library.

All inequalities here use Fraction. The nonnegative local slacks and their
span lower bound are inherited premises, not replayed by this checker.
The companion pressure and independent checkers verify the final count.
"""
from fractions import Fraction as F
from math import comb
import json


def require(condition, label):
    if not condition:
        raise ArithmeticError(label)


def split(row):
    left, right = [row[0]], [row[-1]]
    while len(row) > 1:
        row = [(a+b)/2 for a, b in zip(row, row[1:])]
        left.append(row[0])
        right.append(row[-1])
    return left, right[::-1]


def check():
    alpha, d0, theta = F(627,10**6), F(10016,10**6), F(17,20)
    h0, pi_upper = F(104,125), F(22,7)
    lags = list(range(12,50))
    scale = F(196)
    q = {l:(l+6)//7 for l in lags}
    cuts = {l:F(201*q[l]+17*l,20) for l in lags}
    weights = {l:scale/cuts[l]**2 for l in lags}
    for l in lags:
        offsets = list(range(0,7*(q[l]-1),7))+[l-7]
        require(len(offsets)==len(set(offsets))==q[l], 'Distinct cover starts')
        coverage = [sum(j<=k<j+7 for j in offsets) for k in range(l)]
        require(min(coverage)>=1, 'Every lag gap is covered')
        require(sum(x-1 for x in coverage)==7*q[l]-l, 'Exact overlap count')
        natural_cut = q[l]*d0/alpha-(7*q[l]-l)*theta
        require(cuts[l]>=natural_cut>0, 'Separation credited span cutoff')
        # Multiplicity proof: each relative start gives one possible anchor.
        for m in (l+1, l+7, 2*l+1):
            incidence = [0]*(m-7)
            for i in range(m-l):
                for start in offsets:
                    incidence[i+start]+=1
            require(max(incidence)<=q[l], 'Finite cover incidence bound')

    coeff = [F(v,10**6) for v in (122290,-45012,16241,1887,2711,3731,-5307,2982,-5153,5358)]
    z0 = F(51,5)**2
    tail_lower = 1+sum((min((-1)**j*v,(-1)**j*v*z0/(z0-j*j))
                        for j,v in enumerate(coeff,1)),F(0))
    require(tail_lower==F(3123731672835001,3750959784000000)>h0,
            'Uniform tail factor includes positive limiting terms')

    cheb = [[F(1)],[F(0),F(1)]]
    for l in range(2,max(lags)+1):
        row = [F(0)]*(l+1)
        for k,v in enumerate(cheb[-1]):
            row[k+1]+=2*v
        for k,v in enumerate(cheb[-2]):
            row[k]-=v
        cheb.append(row)
    numerator=[F(0)]*(max(lags)+1)
    numerator[0], numerator[1] = F(-1), F(1)
    for l,w in weights.items():
        numerator[0]+=w
        for k,v in enumerate(cheb[l]):
            numerator[k]-=w*v
    require(sum(numerator)==0, 'Exact phase quotient root')
    polynomial=[]
    running=F(0)
    for v in numerator[:-1]:
        running+=v
        polynomial.append(running)
    require(numerator[-1]==-polynomial[-1], 'Exact quotient identity')
    require(sum(polynomial)==sum(weights[l]*l*l for l in lags)-1,
            'Quotient endpoint cross-check')
    n=len(polynomial)-1
    power=[2**k*sum((polynomial[j]*comb(j,k)*(-1)**(j-k)
                    for j in range(k,n+1)),F(0)) for k in range(n+1)]
    bern=[sum((power[k]*F(comb(j,k),comb(n,k)) for k in range(j+1)),F(0))
          for j in range(n+1)]
    stack=[(bern,F(-1),F(1),0)]
    leaves=[]
    while stack:
        row,lo,hi,depth=stack.pop()
        if min(row)>0:
            lower = F((min(row)*10**9).numerator//(min(row)*10**9).denominator,10**9)
            leaves.append({'interval':[str(lo),str(hi)],'lower':str(lower),'depth':depth})
            continue
        require(depth<16 and len(stack)+len(leaves)<4096,'Bounded Bernstein subdivision')
        a,b=split(row)
        require(a[-1]>0,'Positive exact subdivision midpoint')
        mid=(lo+hi)/2
        stack.extend([(b,mid,hi,depth+1),(a,lo,mid,depth+1)])
    require(min(F(row['lower']) for row in leaves)>F(1,10000),'Positive uniform phase margin')

    # Smooth residual: 1-v^(-2) <= 2(v-1), v>=1, because
    # 2v^3-3v^2+1 = (v-1)^2(2v+1) >= 0.
    # c<=4 and c<=4pi^2*t^2*K^2/h0^2 imply
    # c<=4pi^2*R^2*K^2/h0^2 + 8(t-R)_+/R.
    # L>=alpha(t-R) and L>=0 give (t-R)_+<=L/alpha.
    good = max(weights[l]*pi_upper**2*cuts[l]**2/(8*h0**2) for l in lags)
    bad = sum((weights[l]*q[l]/(2*alpha*cuts[l]) for l in lags),F(0))
    boundary = sum((weights[l]*l/8 for l in lags),F(0))
    require(boundary<F(6),'Rounded finite boundary')
    rigidity=F(350)
    require(good<rigidity and bad<rigidity, 'Exact energy and aggregate slack prices')
    require(bad<F(112),'Displayed aggregate slack bound')
    result=dict(status='PASS: exact smooth aggregate payment certificate',
                lags=[min(lags),max(lags)],scale=str(scale),
                cut='(201*ceil(l/7)+17*l)/20',
                weight='78400/(201*ceil(l/7)+17*l)^2',
                tail_lower=str(tail_lower),h0=str(h0),
                polynomial_degree=n,bernstein_pieces=leaves,
                good_price=str(good),bad_price=str(bad),
                good_display=float(good),bad_display=float(bad),
                rigidity=str(rigidity),boundary=str(boundary),boundary_display=float(boundary),
                rounded_distance_boundary='6',
                acceptance='All rational; floating point values are display only.')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    check()
