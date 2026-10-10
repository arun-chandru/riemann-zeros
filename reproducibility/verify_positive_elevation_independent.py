"""Independent exact audit of the M350 phase and affine pressure certificate.

No import of either first implementation. Chebyshev coefficients use a
closed factorial formula; division is from highest degree; each Bernstein
interval is transformed directly. Standard library, no float acceptance.
"""
from fractions import Fraction as Q
from math import factorial, comb


def need(statement, label):
    if not statement:
        raise ArithmeticError(label)


def chebyshev_closed(n):
    if n == 0:
        return [Q(1)]
    out = [Q(0)] * (n + 1)
    for j in range(n // 2 + 1):
        degree = n - 2*j
        out[degree] = (Q(n, 2) * (-1)**j * factorial(n-j-1)
                       * 2**degree / (factorial(j)*factorial(degree)))
    return out


LAGS = tuple(range(12, 50))
alpha, d0, theta = Q(627, 10**6), Q(10016, 10**6), Q(17, 20)
h0 = Q(104, 125)
q = {l: (l+6)//7 for l in LAGS}
radius = {l: 16*q[l] - (7*q[l]-l)*theta for l in LAGS}
weight = {l: Q(196)/radius[l]**2 for l in LAGS}

# The start list itself proves the multiplicity bound for every length:
# a fixed complete frame start j has at most one anchor j-offset for each
# distinct offset. Each offset is valid inside its anchored lag edge.
for l in LAGS:
    starts = tuple(7*k for k in range(q[l]-1)) + (l-7,)
    need(len(set(starts)) == len(starts) == q[l], "distinct cover starts")
    need(min(starts) == 0 and max(starts) == l-7, "cover endpoints")
    need(all(0 <= s <= l-7 for s in starts), "all complete frames exist")
    mult = [sum(s <= k < s+7 for s in starts) for k in range(l)]
    need(min(mult) >= 1 and sum(mult) == 7*q[l], "coverage and overlap")
    need(radius[l] == Q(201*q[l]+17*l,20), "displayed radius")
    need(radius[l] > q[l]*d0/alpha-(7*q[l]-l)*theta > 0, "span excess payment")

# Closed coefficients are independent of the original recurrence.
num = [Q(0)]*50
num[0], num[1] = -Q(1), Q(1)
for l in LAGS:
    num[0] += weight[l]
    for j, c in enumerate(chebyshev_closed(l)):
        num[j] -= weight[l]*c
need(sum(num) == 0, "exact root at one")
phase = [Q(0)]*49
phase[-1] = -num[-1]
for k in range(48, 0, -1):
    phase[k-1] = phase[k]-num[k]
need(phase[0] == num[0], "descending division by one minus t")
need(sum(phase) == sum(weight[l]*l*l for l in LAGS)-1, "phase endpoint identity")

pieces = (
    (Q(-1), Q(-3,4), Q(164977,50000000)),
    (Q(-3,4), Q(-1,2), Q(255867387,10**9)),
    (Q(-1,2), Q(0), Q(370473121,10**9)),
    (Q(0), Q(1,2), Q(207778923,200000000)),
    (Q(1,2), Q(3,4), Q(743826529,250000000)),
    (Q(3,4), Q(1), Q(1129786507,125000000)),
)
for lo, hi, displayed_lower in pieces:
    # Direct substitution t=lo+(hi-lo)v, with no de Casteljau subdivision.
    mono = [sum((phase[j]*comb(j,k)*lo**(j-k)*(hi-lo)**k
                 for j in range(k,49)), Q(0)) for k in range(49)]
    bern = [sum((mono[k]*Q(comb(j,k),comb(48,k))
                 for k in range(j+1)),Q(0)) for j in range(49)]
    need(min(bern) > displayed_lower > Q(1,10000), "all direct Bernstein coefficients")
print("PASS: independent closed-form/direct-Bernstein phase comparison on all six pieces.")

# Strengthened tail bound. Signs decide which endpoint is used.
coeff = [Q(v,10**6) for v in
         (122290,-45012,16241,1887,2711,3731,-5307,2982,-5153,5358)]
zmin = Q(2601,25)
tail = Q(1)
for j,a in enumerate(coeff,1):
    c = (-1)**j*a
    tail += c*zmin/(zmin-j*j) if c < 0 else c
need(tail == Q(3123731672835001,3750959784000000) > h0, "positive-tail contribution")
need(12*theta == Q(51,5), "every chosen edge lies in tail")

ME = max(weight[l]*Q(22,7)**2*radius[l]**2/(8*h0*h0) for l in LAGS)
MH = sum((weight[l]*q[l]/(2*alpha*radius[l]) for l in LAGS),Q(0))
boundary = sum((weight[l]*l/8 for l in LAGS),Q(0))
need(ME == Q(1890625,5408) < 350, "ordered energy price")
need(MH < 112 < 350 and boundary < 6, "slack price and boundary")
print("PASS: aggregate cover payment constants M=350 and B=6.")

# Signed continuum pressure proof, independently reconstructing derivatives.
def dsinc(k):
    return Q(0) if not k else Q(1 if k % 2 == 0 else -1,k)

def deriv(k):
    return dsinc(k)+sum((a*(dsinc(k-j)+dsinc(k+j))/2
                        for j,a in enumerate(coeff,1)),Q(0))

p, c, C = Q(17,4000), Q(1891,200000), Q(19,4)
need(sum(map(abs,coeff)) < 1, "profile positive")
need(Q(22,7)**2*(1+sum(map(abs,coeff)))/3 < 4, "global second derivative bound")
need(Q(5,2)*p > c, "large nearest integers")
for k,r in ((1,Q(1,20)),(2,Q(1,25))):
    value, slope = coeff[k-1]/2, deriv(k)
    A, B = abs(value), abs(slope)
    need(value*slope < 0, "signed direction")
    need(A-B*r-2*r*r > 0, "right Taylor positivity")
    a0 = 2*A*A+p*k-c
    a1 = p-4*A*B
    a2 = 2*B*B-8*A+C
    need(a0 > 0 and a2 > 0 and 4*a0*a2-a1*a1 > 0, "right inner quadratic")
    displayed_discriminants = {
        1: Q(53039685134812231,193600000000000000000),
        2: Q(23907572162979771241,948640000000000000000000),
    }
    need(4*a0*a2-a1*a1 == displayed_discriminants[k], "published discriminant margin")
    need(B-2*r > 0 and 4*A*(B-2*r) > p, "left inner bound")
    displayed_slopes = {
        1: Q(2804524299211,13200000000000),
        2: Q(1835448537709,70000000000000),
    }
    need(4*A*(B-2*r)-p == displayed_slopes[k], "published left slope")
    need(2*C*r > p and p*(k-r)+C*r*r > c, "outer monotone pressure")
    need(p*(k-r)+C*r*r-c == (Q(2583,400000) if k==1 else Q(259,40000)),
         "published outer pressure margin")
need(deriv(1) == -Q(260321359,264000000), "first derivative display")
need(deriv(2) == Q(2577793683,6160000000), "second derivative display")
print("PASS: independent signed-Taylor proof covers every real gap g>=17/20.")

# The final frozen scalars; these replace all earlier positive increments.
a, b = Q(8722,10**6), Q(4389,10**6)
M = Q(350)
mix = 1+M*C
new_a = (M*C*a+c)/mix
new_b = (M*C*b+p)/mix
B0 = 7*a+Q(3882,10**6)
new_boundary = (M*C*B0+c+6*C)/mix
smoothing = (M*C*270+6)/mix
need(mix == Q(3327,2), "bootstrap denominator")
need(new_a == Q(725489,83175000), "final a")
need(new_b == Q(584077,133080000), "final b")
need(new_boundary == Q(27293111,332700000), "affine endpoint loss")
need(smoothing == Q(299254,1109), "45-band smoothing loss")
need(0 < new_a < Q(13707941,2**28) and 0 < new_b < b, "close-pair and span guards")
print("PASS: exact affine floor a=725489/83175000, b=584077/133080000 and transfer guards.")
print("All acceptance computations used exact rational arithmetic.")
