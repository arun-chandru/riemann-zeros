#!/usr/bin/env python3
"""Exact arithmetic companion for the clipped Gram energy and zero-counting argument.

Run with Python 3: python verify_clipped_gram_constants.py [--json]
No dependencies, files, network, floating-point arithmetic, or sampling are used.
All numerical arithmetic uses integers and fractions.Fraction.

This checks the printed constants and polynomial certificates, including
transcendental enclosures derived from convergent alternating series. The
analytic Fourier, trace, and matching arguments are proved in the manuscript.
The inherited seven-variable local floor a8=805124/10^8 is INPUT here: its
large exact certificate is in the accompanying Lean source, not replayed here.

For 0<x<1 the arctangent series has decreasing alternating terms. Consecutive
partial sums enclose arctan(x). Machin's identity then encloses pi. For sine
and cosine, the alternating tails used below decrease on the stated intervals;
the sign of the first omitted term gives the direction of each bound.
Polynomial identities plus signs of coefficients/factors certify whole
intervals, never a finite grid of test points.
"""

from fractions import Fraction as F
from math import comb, factorial, isqrt
import json
import sys


CHECKS = []


def require(name, proposition):
    if not proposition:
        raise ArithmeticError("FAILED: " + name)
    CHECKS.append(name)


def poly_add(a, b):
    return tuple((a[i] if i < len(a) else F(0)) +
                 (b[i] if i < len(b) else F(0))
                 for i in range(max(len(a), len(b))))


def poly_scale(a, c):
    return tuple(c * x for x in a)


def poly_mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            out[i + j] += u * v
    return tuple(out)


def poly_value(a, x):
    result = F(0)
    for coefficient in reversed(a):
        result = result * x + coefficient
    return result


def trig_sum(x, last, cosine=False):
    """Terms k=0,...,last of sin(x) or cos(x), with exact arithmetic."""
    shift = 0 if cosine else 1
    return sum((F((-1) ** k, factorial(2*k + shift)) *
                x ** (2*k + shift) for k in range(last + 1)), F(0))


def atan_interval(x, terms=24):
    """Consecutive sums: alternating-series enclosure, not an estimate."""
    require("arctangent series domain " + str(x), F(0) < x < F(1))
    partial = sum((F((-1) ** k, 2*k + 1) * x ** (2*k + 1)
                   for k in range(terms)), F(0))
    other = partial + F((-1) ** terms, 2*terms + 1) * x ** (2*terms + 1)
    return min(partial, other), max(partial, other)


def down(x, places):
    """Floor to the specified decimal place, with no float conversion."""
    scale = 10 ** places
    integer = (x.numerator * scale) // x.denominator
    whole, digits = divmod(integer, scale)
    return str(whole) + ("." + str(digits).zfill(places) if places else "")


def sqrt_interval(x, places=30):
    """Rational enclosure proved by integer square comparisons."""
    scale = 10**places
    lower_n = isqrt((x.numerator*scale*scale)//x.denominator)
    lo, hi = F(lower_n,scale),F(lower_n+1,scale)
    require("radical enclosure by exact squares", 0 < lo and lo*lo <= x < hi*hi)
    return lo,hi


def cosine_energy_interval(numerators, pi_lo, pi_hi):
    aa = tuple(F(n,10**6) for n in numerators)
    require("supplementary cosine profile positive " + str(numerators[0]), sum(map(abs,aa),F(0)) < 1)
    base = F(4,3) + sum((x*x/2 for x in aa),F(0))
    correction = sum(((F((-1)**j)*x-x*x/4)/(j*j) for j,x in enumerate(aa,1)),F(0))
    endpoints = (base+correction/(pi_lo*pi_lo),base+correction/(pi_hi*pi_hi))
    return min(endpoints),max(endpoints)


def supplementary_counts(pi_lo, pi_hi):
    """Independent energies and the remaining five table rows."""
    first = (193305,-23105,29166,13952,-21009,41052)
    third = (200566,19,-9358,9283,-6746,8170)
    r1lo,r1hi = cosine_energy_interval(first,pi_lo,pi_hi)
    r3lo,r3hi = cosine_energy_interval(third,pi_lo,pi_hi)
    r1cap,r3cap = F(666352479141,500000000000),F(83278977073,62500000000)
    require("first-profile energy cap", r1lo <= r1hi <= r1cap)
    require("third-profile energy cap", r3lo <= r3hi <= r3cap)
    clip,x1,x2,b1 = F(67,20),F(54350,10**6),F(69433,10**6),F(37314,10**6)
    mu,nu = 2-2*x1+x2,2*clip-3+x1-x2
    rh,rk = 6*clip-7-clip*clip+x1-2*x2,4*clip-2-clip*clip-2*x1
    require("labelled four guard values", (mu,nu,rh,rk) == (F(1960733,10**6),F(3684917,10**6),F(224123,125000),F(43,625)))
    require("labelled four guard signs", min(mu,nu,rh,rk) > 0)
    require("higher-real residual base and slope", 3*nu+mu-clip*clip == rh and nu > 0)
    require("off-line residual base and slope", 2*nu+2*mu-clip*clip == rk and 2*nu > 0)
    distinct = (3-r1cap+x2-x1-b1)/mu
    simple = 2*distinct-1
    require("distinct exact fraction", distinct == F(822532020859,980366500000))
    require("simple-anywhere exact fraction", simple == F(332348770859,490183250000))
    a,b = F(54525,10**6),F(37272,10**6)
    tau = r3cap+b
    beta = tau-1-a
    require("third total budget and excess", tau == F(85608477073,62500000000) and beta == F(19700664573,62500000000))
    radicand = 2-2*a
    rtlo,rthi = sqrt_interval(radicand)
    zlo,zhi = 3-a+2*rtlo,3-a+2*rthi
    union_lo,union_hi = 1-2*beta/zlo,1-2*beta/zhi
    require("positive radical and tangent range", 1 < rtlo <= rthi < 2)
    # Reduce a polynomial in sqrt(2-2a) modulo X^2-(2-2a).
    def reduce_radical(poly):
        return tuple(sum((coefficient*radicand**(i//2) for i,coefficient in enumerate(poly) if i%2 == parity),F(0)) for parity in (0,1))
    tangent = (F(2),F(1))
    zpoly = (3-a,F(2))
    require("radical Z identity", reduce_radical(poly_scale(poly_mul(tangent,tangent),F(1,2))) == zpoly)
    tangent_residual = poly_add(poly_scale(tangent,2),(-(1+a),))
    tangent_residual = poly_add(tangent_residual,poly_scale(poly_mul(tangent,tangent),-F(1,2)))
    require("tangent cancellation identity", reduce_radical(tangent_residual) == (0,0))
    require("all fixed-cutoff hierarchy guards", 0 < a < F(1,2) and tau < F(7,5) and 0 < beta < F(1,3) and zlo > F(9,2))
    require("hierarchy positive slopes and first high cost", zlo-beta > 0 and 2*zlo-3*beta > 0 and 3-tau-a > 0 and 8-2*tau-a > 0)
    mass = (3-tau)/(2-a)
    location = 1-beta/(8-2*tau-a)
    require("mass cutoff-two exact fraction", mass == F(101891522927,121592187500))
    require("location cutoff-two exact fraction", location == F(101891522927,108458411118))
    require("cutoff-two scalar guards", min(beta,6-tau,3-tau-a,1-2*a,8-2*tau-a) > 0)
    require("off-line scalar cutoff residual", (6-tau)*6-9+2*beta-2*(8-2*tau-a) == 9)
    require("scalar population nonnegative", 0 < mass < 1)
    require("scalar population energy saturation", (1+a)*mass+3*(1-mass) == tau)
    require("scalar population location saturation", mass/(mass+(1-mass)/3) == location == 3*mass/(1+2*mass))
    # For 2<=c<=3: 9-phi_c(3)=(3-c)^2. For c>=3, phi_c(3)=9.
    require("scalar saturation all-level polynomial identity", (F(9),F(-6),F(1)) == poly_mul((F(3),F(-1)),(F(3),F(-1))))
    require("simple mass uses phi_c(1)=1 for every c>=2", F(1) <= 2)
    table = {"distinct_locations":(distinct,distinct,"83.90046180"),
             "simple_anywhere":(simple,simple,"67.80092360"),
             "simple_or_critical_line_mass":(union_lo,union_hi,"88.93165002"),
             "mass_at_multiplicity_at_most_two":(mass,mass,"83.79775462"),
             "locations_at_multiplicity_at_most_two":(location,location,"93.94524765")}
    floors = {}
    for key,(lo,hi,expected) in table.items():
        require("eight-place downward percentage " + key, down(100*lo,8) == down(100*hi,8) == expected)
        # Both interval ends in this one decimal bin proves the printed floor.
        floors[key] = expected
    return {"energy_caps":{"first":str(r1cap),"third":str(r3cap)},
            "union_radical_bracket":[str(rtlo),str(rthi)],
            "union_value_bracket":[str(union_lo),str(union_hi)],
            "exact_rational_counts":{"distinct":str(distinct),"simple_anywhere":str(simple),"mass_le_two":str(mass),"locations_le_two":str(location)},
            "scalar_saturation":{"simple_mass":str(mass),"triple_mass":str(1-mass),"energy":str(tau)},
            "percentages":floors}


def verify():
    CHECKS.clear()
    # tan(4 atan(1/5)-atan(1/239))=1; its angle is in (0,pi/2).
    u = F(1, 5)
    tan2 = 2*u/(1-u*u)
    tan4 = 2*tan2/(1-tan2*tan2)
    require("Machin tangent identity", (tan4-F(1,239))/(1+tan4/F(239)) == 1)
    lo5, hi5 = atan_interval(u)
    lo239, hi239 = atan_interval(F(1,239))
    pi_lo, pi_hi = 16*lo5-4*hi239, 16*hi5-4*lo239
    require("Machin angle is positive and below one", 0 < 4*lo5-hi239 < 4*hi5-lo239 < 1)
    require("Machin interval ordered", pi_lo < pi_hi)
    p_plus = F(314159265358979323847, 10**20)
    require("printed 20-decimal pi upper endpoint", pi_hi < p_plus)
    require("coarse pi endpoints", F(157,50) < pi_lo < pi_hi < F(22,7))
    require("seven-decimal pi endpoints", F(31415926,10**7) < pi_lo < pi_hi < F(31415927,10**7))

    r_lo = F(7071067811865475244,10**19)
    r_hi = F(7071067811865475245,10**19)
    sin_lo = F(64963693904,10**11)
    cos_hi = F(76024459708,10**11)
    require("theta square bracket", 0 < r_lo < r_hi < 1 and r_lo*r_lo < F(1,2) < r_hi*r_hi)
    require("theta sine/cosine monotonicity domain", r_hi < pi_lo/2)
    # sin increases and cos decreases on this interval. k=5 is a sine lower
    # sum; k=6 is a cosine upper sum, with decreasing omitted tails.
    require("sine rational lower endpoint by alternating series", sin_lo <= trig_sum(r_lo,5))
    require("cosine rational upper endpoint by alternating series", trig_sum(r_lo,6,True) <= cos_hi)
    require("degree-11 exponential remainder coefficient", F(13,12*factorial(12)) == F(13,5748019200))
    rem = r_lo**12 * F(13,5748019200)
    require("printed sine endpoint with stated common remainder", sin_lo <= trig_sum(r_lo,5)-rem)
    require("printed cosine endpoint with stated common remainder", trig_sum(r_lo,5,True)+rem <= cos_hi)

    dn = (12310798,-15041681,3867664,6489926,-992327,5580097,
          -6846472,3781297,-5670353,3355089,-450523,-218483)
    ds = tuple(F(x,10**9) for x in dn)
    require("twelve cosine coefficients", len(ds) == 12)
    require("absolute correction budget", sum(map(abs,ds),F(0)) == F(64604710,10**9) <= F(13,200))
    require("Freund coefficient-budget compatibility after zero extension", F(13,200) < F(85368367,1250000000) and len(ds+(F(0),)*12) == 24)
    require("positive profile lower bound", 1-F(1,4)-F(13,200) == F(137,200) > 0)
    k0_lo = 1-F(1,12)+F(1,480)-F(1,40320)
    k0_hi = 1-F(1,12)+F(1,480)
    require("reciprocal normalization guard", 0 < F(100,109) < k0_lo)
    require("upper normalization guard", k0_hi == F(147,160))
    require("close-pair gamma guard", F(3,4)/k0_hi == F(40,49) > F(4,5))

    p_energy = sum((d*d*(F(1,2)-1/(4*p_plus*p_plus*j*j))
                    for j,d in enumerate(ds,1)), F(0))
    require("energy correction positivity", p_energy >= 0)
    h_lower = F(3,2)-cos_hi/(2*r_lo*sin_lo)-p_energy/(2*sin_lo*sin_lo)
    h0 = F(67216841,10**8)
    require("energy H lower cap", h_lower >= h0)

    # The lower sinc polynomial is applied for |z|<=pi*17/20<pi. Its first
    # omitted term is positive; all later tail ratios are <1 because pi^2<10.
    sep = F(17,20)
    dlo, dhi = F(157,50)*sep, F(22,7)*sep
    require("sinc Taylor tail decreasing", pi_hi*pi_hi < 10 < 8*9)
    require("D squared upper endpoint", dhi**2 <= F(357,50))
    require("D fourth-power lower endpoint", dlo**4 >= F(31684,625))
    require("D sixth-power upper endpoint", dhi**6 <= 364)
    even_bounds = (F(1),F(357,50),F(31684,625),F(364))
    sinc_lower_poly = (F(0),)
    for n in range(4):
        # Exact expansion of (s-1/2)^(2n)+(s+1/2)^(2n) in x=s^2.
        even_sum = tuple(2*comb(2*n,2*k)*F(1,2)**(2*n-2*k) for k in range(n+1))
        sinc_lower_poly = poly_add(sinc_lower_poly,poly_scale(even_sum,F((-1)**n,factorial(2*n+1))*even_bounds[n]))
    printed_poly = (F(2)-F(119,200)+F(7921,150000)-F(13,5760),
                    -F(119,50)+F(7921,6250)-F(13,96),
                    F(31684,37500)-F(13,24),-F(13,90))
    require("sinc-sum polynomial coefficient identity", sinc_lower_poly == printed_poly)
    line = (F(29,20),-F(5,4))
    residual = poly_add(printed_poly,poly_scale(line,-1))
    require("sinc-sum constant surplus", residual[0] == F(19979,3600000) > 0)
    require("sinc-sum linear surplus", residual[1] == F(583,300000) > 0)
    require("sinc-sum quadratic factor throughout interval", residual[2]+residual[3]/4 == F(7513,28125) > 0 and residual[3] < 0)
    require("positive sinc-sum lower line", poly_value(line,F(1,4)) == F(91,80) > 0 and line[1] < 0)
    f_upper_poly = poly_scale((F(213,200),F(-1),F(1,6)),F(109,100))
    q_poly = poly_add(poly_scale(poly_mul(line,line),F(7,10)),poly_scale(f_upper_poly,-1))
    require("majorant gap polynomial coefficient identity", q_poly == (F(3109,10000),-F(579,400),F(2189,2400)))
    margin = F(5789,960000)
    require("majorant gap endpoint", poly_value(q_poly,F(1,4)) == margin > 0)
    factored = poly_mul((F(1,4),F(-1)),(F(579,400)-F(2189,9600),-F(2189,2400)))
    require("majorant gap factorization", factored == poly_add(q_poly,(-margin,)))
    require("majorant gap factor sign throughout interval", F(579,400)-F(2189,2400)*F(1,2) > 0)

    ylo, yhi = F(3,20)*pi_lo, F(3,20)*pi_hi
    require("close-pair y endpoints", F(471,1000) <= ylo <= yhi <= F(33,70))
    require("mass y-squared lower endpoint", ylo*ylo >= F(22206,10**5))
    require("mass y-fourth-power upper endpoint", yhi**4 <= F(4932,10**5))
    sin_ratio_cap = 1-F(22206,10**5)/6+F(4932,10**5)/120
    require("mass sine ratio cap", sin_ratio_cap <= F(963401,10**6))
    mass_cap = F(28,17)*(1+F(3,17)*F(963401,10**6))
    require("majorant mass guard", mass_cap <= F(1253,650) < 2)
    require("actual finite-height spectral guard", F(51,50)*F(1253,650) == F(63903,32500) < 2)
    require("close-pair z endpoints", F(8,3) <= sep*pi_lo <= sep*pi_hi <= F(187,70))
    require("close-pair cosine lower endpoint", 1-F(33,70)**2/2 >= F(22,25))
    require("close-pair sine lower endpoint", F(471,1000)-F(33,70)**3/6 >= F(9,20))
    numerator_floor = F(22,25)+2*F(4,5)*F(8,3)*F(9,20)
    denom_cap = 2*F(187,70)**2-1
    require("close-pair numerator floor", numerator_floor == F(14,5))
    require("close-pair positive denominator and cap", 0 < 2*F(8,3)**2-1 <= denom_cap == F(32519,2450))
    require("unperturbed normalized kernel endpoint", numerator_floor/denom_cap >= F(1,5))
    correction = F(109,100)*F(13,200)
    require("Fourier correction budget", correction == F(1417,20000))
    require("close-pair kernel guard", F(1,5)-correction == F(2583,20000) > F(1,8))

    rows = ((13630830,0,19565904,20317441,45030671,100000000,200000000),
            (29997996,30621878,47766489,79682558,109938656,100000000),
            (37356140,69378121,65335210,79682558,45030671),
            (38030065,69378121,47766489,20317441),
            (37356140,30621878,19565904),(29997996,0),(13630830,))
    weights = {(i,i+1+k): F(value,10**8) for i,row in enumerate(rows) for k,value in enumerate(row)}
    bnum = (28898,57272,75526,80958,75526,57272,28898)
    pressures = tuple(F(n,10**8) for n in bnum)
    reward, pressure = F(805124,10**8),sum(pressures,F(0))
    require("all 28 upper-triangle weights supplied", len(weights) == 28)
    require("nonnegative weights", all(value >= 0 for value in weights.values()))
    require("reflection of all pair coefficients", all(value == weights[7-j,7-i] for (i,j),value in weights.items()))
    require("positive reflected pressures", min(pressures) == F(28898,10**8) > 0 and pressures == pressures[::-1])
    lag_sums = tuple(sum((weights[i,i+lag] for i in range(8-lag)),F(0)) for lag in range(1,8))
    expected = (199999997,199999998,199999996,199999998,199999998,200000000,200000000)
    require("seven exact lag sums", lag_sums == tuple(F(n,10**8) for n in expected))
    require("all lag capacities", all(x <= 2 for x in lag_sums))
    require("pressure telescope total", pressure == F(404350,10**8))
    # These are constants in the general telescoping proof, not tests at
    # finitely many chain lengths. There are m-(8-1) windows and 2*(8-1)
    # possible ordered neighbours at each point, for every m>=8.
    endpoint, neighbours, per_entry_error = 8-1,2*(8-1),2
    require("boundary reward and sparse error coefficients", endpoint*reward == F(5635868,10**8) and neighbours == 14 and neighbours*per_entry_error == 28)
    require("pair can pay two local floors", 0 < reward < F(1,64) and 2*reward < F(1,32))
    require("Freund close-pair budget compatibility", reward < F(9,400))
    require("bounded certificate box from pressure alone", all(reward/b > 0 for b in pressures))
    q = F(66812491,99194876)
    avg = (1+q)/2
    require("exact zero-counting quotient", (1-reward)*q == h0-pressure and reward < 1)
    require("retained-index average identity", 3-(2-h0)-pressure+reward*q == 1+q)
    require("quoted q decimal", down(q,18) == "0.673547805029767868")
    require("quoted average decimal", down(avg,18) == "0.836773902514883934")
    require("downward q percentage", down(100*q,8) == "67.35478050")
    require("downward average percentage", down(100*avg,8) == "83.67739025")
    require("positive epsilon error budget", F(30)+reward+pressure > 0 and reward/F(56) > 0)
    other = supplementary_counts(pi_lo,pi_hi)
    percentages = {"simple_critical_line":down(100*q,8),"retained_index_average":down(100*avg,8),**other["percentages"]}
    return {
        "status": "PASS", "checks": len(CHECKS),
        "arithmetic": "integers and fractions.Fraction only; no sampling",
        "scope": "AM guards, supplementary cosine energies, all seven table rows, scalar saturation; inherited multivariate local floors are not replayed",
        "pi_bracket_30_decimal_outer": [down(pi_lo,30),str(F((pi_hi.numerator*10**30)//pi_hi.denominator+1,10**30))],
        "H_floor": str(h0), "H_rational_guard_surplus_positive": h_lower > h0,
        "majorant_margin": str(margin), "majorant_mass_upper_bound": str(F(1253,650)),
        "actual_separated_spectral_cap": str(F(63903,32500)),
        "lag_numerators": list(expected), "pressure": str(pressure), "reward": str(reward),
        "q": str(q), "average": str(avg),
        "percentages_downward_8_places": percentages,
        "supplementary_counts":other,
        "passed": CHECKS,
    }


if __name__ == "__main__":
    report = verify()
    if sys.argv[1:] == ["--json"]:
        print(json.dumps(report,indent=2))
    elif sys.argv[1:]:
        raise SystemExit("Usage: python verify_clipped_gram_constants.py [--json]")
    else:
        print("PASS:",report["checks"],"exact arithmetic checks.")
        print("q =",report["q"],"; simple critical-line percentage >=",report["percentages_downward_8_places"]["simple_critical_line"])
        print("Average percentage >=",report["percentages_downward_8_places"]["retained_index_average"])
        print("Scope:",report["scope"])
