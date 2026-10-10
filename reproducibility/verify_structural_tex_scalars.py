"""Exact checks of the scalar expressions printed in the structural TeX.

Standard library only. The two-lattice determinant, elementary endpoint
positivity, 15/16 union certificate and 9/10 root signs are reconstructed
without importing the original checkers.
"""
from fractions import Fraction as F
from math import factorial


def require(condition, text):
    if not condition:
        raise ArithmeticError(text)


lo, hi = F(157, 50), F(22, 7)
aL = F(3, 16) + 1/(8*hi) - 2/lo**2
bL = F(1, 8) - 2/(9*lo**2)
cU2 = (F(1, 4)-1/hi**2)/64
det = (aL-F(1,1000))*(bL-F(1,1000))-cU2
require(aL > F(1,1000) and bL > F(1,1000), "positive shifted diagonals")
require(det == F(17310471944257,330823607944500000) > 0,
        "displayed shifted determinant")
require(F(1,8)-1/(2*lo**2) > F(1,1000), "even coercivity")
endpoint = F(389,606)-F(2646,4379)
require(endpoint == F(99955,2653674) > 0, "elementary positive endpoint")
require(F(384,359)+F(50,231) < F(4,3), "root uniqueness sup norm bound")


def cos_sum(z, n):
    return sum(((-z)**k / factorial(2*k) for k in range(n+1)),F(0))


def sinc_sum(z, n):
    return sum(((-z)**k / factorial(2*k+1) for k in range(n+1)),F(0))


def minimum_bounds(u, lower_index):
    z = 1/(16*u)
    require(0 < 2*z <= 1, "alternating-series argument")
    n = lower_index
    require(n % 2 == 1, "odd lower index")
    cL, cU = cos_sum(z,n), cos_sum(z,n+1)
    sL, sU = sinc_sum(z,n), sinc_sum(z,n+1)
    c2L, c2U = cos_sum(2*z,n), cos_sum(2*z,n+1)
    s2L, s2U = sinc_sum(2*z,n), sinc_sum(2*z,n+1)
    require(min(cL,sL,c2L,s2L) > 0, "positive Taylor denominator bounds")
    lower = F(1,2)+u*(c2L/s2U-2*z*sU/cL)
    upper = F(1,2)+u*(c2U/s2L-2*z*sL/cU)
    require(lower < upper, "ordered minimum interval")
    return lower, upper


uU, eU = F(13711914,10**8),F(5251083389778,10**14)
require(F(1,8) < uU < F(1,3) and 0 < 2*eU < uU, "union count guards")
mL, _ = minimum_bounds(uU,15)
require(mL-2*uU*uU/(uU-eU) > F(1,10**14), "printed strict union margin")
require(1-2*eU == F(22374458305111,25000000000000), "union ceiling")
left, right = F(32349386085,10**11), F(32349386086,10**11)
require(minimum_bounds(left,9)[0] > 2*left, "root left endpoint sign")
require(minimum_bounds(right,9)[1] < 2*right, "root right endpoint sign")

qold,qnew = F(337017,500000),F(16851,25000)
require(F(77,250)-2/lo**2 > F(1,10), "all-mode coercivity")
require(F(651926252580,10**12) < 2-2*qold-F(1,200000),
        "coarse old test failure")
require(F(651920180779,10**12)-F(160,10**24) > 2-2*qnew+F(18,10**8),
        "coarse residual lower bound")
eps=F(1,10**11)
require(4*(eps+10*eps+10*eps**2) < F(1,10**9), "minimizer uniform error")

print("PASS: structural TeX determinant, endpoint positivity, union margin > 10^-14,")
print("9/10 root signs, and coarse twelve-fibre residual arithmetic.")
print("Every acceptance comparison used exact rational arithmetic.")
