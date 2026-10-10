'use strict';

// Finite rational certificate for the new all-profile joint-affine union
// obstruction. No floating-point value participates in acceptance.
const abs = x => x < 0n ? -x : x;
function gcd(a,b) { a=abs(a); b=abs(b); while(b) { const c=a%b; a=b; b=c; } return a; }
function Q(n,d=1n) { if(!d) throw Error('Zero denominator'); if(d<0n) {n=-n;d=-d;} const g=gcd(n,d); return [n/g,d/g]; }
const add=(a,b)=>Q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);
const neg=a=>[-a[0],a[1]];
const sub=(a,b)=>add(a,neg(b));
const mul=(a,b)=>Q(a[0]*b[0],a[1]*b[1]);
const div=(a,b)=>Q(a[0]*b[1],a[1]*b[0]);
const sq=a=>mul(a,a);
const pos=a=>a[0]>0n;
const check=(p,s)=>{if(!p) throw Error(s);};
const ONE=Q(1n), TWO=Q(2n);

// C(z)=cos(sqrt(z)), S(z)=sin(sqrt(z))/sqrt(z).
// For 0<z<=1, the index-15 partial sum is lower, index-16 upper.
function series(z,offset,last) {
  let sum=ONE, term=ONE;
  for(let k=1;k<=last;k++) {
    term=div(mul(neg(term),z),Q(BigInt((2*k+offset-1)*(2*k+offset))));
    sum=add(sum,term);
  }
  return sum;
}
function bounds(u) {
  const z=div(ONE,mul(Q(16n),u)), z2=mul(TWO,z);
  check(pos(z)&&!pos(sub(z2,ONE)),'Taylor argument in (0,1]');
  const cL=series(z,0,15),cU=series(z,0,16),sL=series(z,1,15),sU=series(z,1,16);
  const c2L=series(z2,0,15),c2U=series(z2,0,16),s2L=series(z2,1,15),s2U=series(z2,1,16);
  for(const x of [cL,cU,sL,sU,c2L,c2U,s2L,s2U]) check(pos(x),'Positive Taylor factors');
  const acL=div(c2L,s2U), acU=div(c2U,s2L);
  const btL=mul(z,div(sL,cU)),btU=mul(z,div(sU,cL));
  return {
    endpointL:sub(acL,btU),
    mL:add(Q(1n,2n),mul(u,sub(acL,mul(TWO,btU)))),
    mU:add(Q(1n,2n),mul(u,sub(acU,mul(TWO,btL))))
  };
}

// The classical bounds 3.14 < pi < 22/7 suffice for the strengthened
// zero-integral coercivity F_{1/8}(h) >= (1/1000)||h||_2^2.
const piL=Q(157n,50n),piU=Q(22n,7n),base=Q(1n,8n),c=Q(1n,1000n);
const aL=sub(add(Q(3n,16n),div(base,piU)),div(TWO,sq(piL)));
const bL=sub(base,div(TWO,mul(Q(9n),sq(piL))));
const off2U=mul(sq(base),sub(Q(1n,4n),div(ONE,sq(piU))));
check(pos(sub(aL,c))&&pos(sub(bL,c)),'Positive shifted odd diagonals');
check(pos(sub(mul(sub(aL,c),sub(bL,c)),off2U)),'Positive shifted odd determinant');
check(pos(sub(sub(base,div(ONE,mul(TWO,sq(piL)))),c)),'Even coercivity');
check(pos(bounds(base).endpointL),'Positive minimizing profile at u=1/8');

// One common rational count point. Near-optimality is not an acceptance
// premise; the strict rational inequality alone certifies the ceiling.
const u=Q(13711914n,100000000n),e=Q(5251083389778n,100000000000000n);
check(pos(sub(u,base))&&pos(sub(Q(1n,3n),u)),'u in proven variational interval');
check(pos(e)&&pos(sub(u,mul(TWO,e))),'Combinatorial guard 0<2e<u');
const mb=bounds(u),need=div(mul(TWO,sq(u)),sub(u,e));
check(pos(sub(mb.mL,need)),'Universal row residual is strictly positive');
const unionCap=sub(ONE,mul(TWO,e));
check(pos(sub(Q(9n,10n),unionCap)),'Union ceiling below 90 percent');
check(sub(unionCap,Q(22374458305111n,25000000000000n))[0]===0n,'Exact displayed ceiling');

console.log('PASS: exact rational coercivity, profile positivity, and common-row residual.');
console.log('Common witness: u = 0.13711914, e = 0.05251083389778.');
console.log('Every covered union deduction is at most 0.89497833220444.');
console.log('This is a method ceiling, not an upper bound on the actual zero union mass.');
