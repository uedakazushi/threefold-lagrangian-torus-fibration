#!/usr/bin/env python3
"""Finite independent checks for the relative elliptic-to-pants construction.

The mathematical proofs, including all existence claims, are in the TeX files.
These standard-library calculations check constants, orientations and examples.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import cmath
import json
import math


def simpson(fn, a, b, n):
    assert n % 2 == 0
    step = (b-a)/n
    return step/3 * (fn(a)+fn(b)+sum(
        (4 if i % 2 else 2)*fn(a+i*step) for i in range(1,n)))


def adaptive(fn, a, b, tol=1e-10, depth=32):
    fa, fb, fm = fn(a), fn(b), fn((a+b)/2)
    whole = (b-a)*(fa+4*fm+fb)/6
    def sub(a,b,fa,fm,fb,whole,tol,depth):
        m=(a+b)/2
        fl,fr=fn((a+m)/2),fn((m+b)/2)
        left=(m-a)*(fa+4*fl+fm)/6
        right=(b-m)*(fm+4*fr+fb)/6
        error=left+right-whole
        if abs(error)<=15*tol or depth==0:
            return left+right+error/15
        return (sub(a,m,fa,fl,fm,left,tol/2,depth-1)
                +sub(m,b,fm,fr,fb,right,tol/2,depth-1))
    return sub(a,b,fa,fm,fb,whole,tol,depth)


def lobe_area(nu,t,nr=96,nt=256):
    radius=math.sqrt(t)
    def angular_average(r):
        return sum(r/abs(nu+radius*r*cmath.exp(2j*math.pi*k/nt))
                   for k in range(nt))/nt
    return t/2*simpson(angular_average,0,1,nr)


def outer_period(nu,t):
    radius=math.sqrt(t)
    difference=radius-nu
    fn=lambda u: 1/math.hypot(difference,2*math.sqrt(nu*radius)*math.sin(u/2))
    return adaptive(fn,0,math.pi,tol=1e-10)/(2*math.pi)


def cutoff(p,R=3):
    lo,hi=R/3,2*R/3
    if p<=lo: return 1.0,0.0
    if p>=hi: return 0.0,0.0
    t=(p-lo)/(hi-lo)
    logratio=-1/t+1/(1-t)
    if logratio>700: return 0.0,0.0
    if logratio< -700: return 1.0,0.0
    increase=1/(1+math.exp(-logratio))
    derivative=increase*(1-increase)*(1/t**2+1/(1-t)**2)/(hi-lo)
    return 1-increase,-derivative


def differential_checks():
    R,nu=3.0,0.02
    minimum=float('inf')
    derivative_error=0.0
    count=0
    for p in [1+i/80 for i in range(81)]:
        chi,dchi=cutoff(p,R)
        for theta in [i/32 for i in range(32)]:
            cosine=math.cos(4*math.pi*theta)
            exact=2*p-2*nu*(chi+p*dchi)*cosine
            def h(q): return q*q-2*nu*q*cutoff(q,R)[0]*cosine+nu*nu
            eps=1e-6
            fd=(h(p+eps)-h(p-eps))/(2*eps)
            derivative_error=max(derivative_error,abs(exact-fd))
            minimum=min(minimum,exact)
            count+=1
    assert minimum>0 and derivative_error<2e-8
    return {'sampled_derivatives':count,'minimum_radial_derivative':minimum,
            'maximum_finite_difference_error':derivative_error,
            'note':'The uniform derivative proof is analytic; this is one smooth-cutoff sample.'}


def lobe_checks():
    records=[]
    refinement=0.0
    for nu in [0.08,0.3,1.0]:
        for ratio in [0.0625,0.25,0.64]:
            t=ratio*nu*nu
            value=lobe_area(nu,t)
            refined=lobe_area(nu,t,192,512)
            refinement=max(refinement,abs(value-refined))
            lower=t/(4*(nu+math.sqrt(t)))
            upper=t/(4*(nu-math.sqrt(t)))
            assert lower<refined<upper
            records.append({'nu':nu,'t_over_nu_squared':ratio,'area':refined,
                            'lower_bound':lower,'upper_bound':upper})
    saddle_error=max(abs(simpson(lambda theta:2*nu*math.cos(4*math.pi*theta),
                                -0.125,0.125,1024)-nu/math.pi)
                     for nu in [0.08,0.3,1.0])
    assert refinement<1e-8 and saddle_error<1e-11
    periods=[]
    for nu in [0.1,0.4,1.0]:
        values=[outer_period(nu,nu*nu*(1+delta)) for delta in [1e-2,1e-3,1e-4,1e-5]]
        assert all(a<b for a,b in zip(values,values[1:]))
        slope=(values[-1]-values[-2])/math.log(10)
        expected=1/(2*math.pi*nu)
        assert abs(slope/expected-1)<0.002
        periods.append({'nu':nu,'delta':[1e-2,1e-3,1e-4,1e-5],
                        'periods':values,'log_slope':slope,'limiting_log_slope':expected})
    return {'area_examples':records,'max_quadrature_refinement_difference':refinement,
            'maximum_saddle_area_error':saddle_error,'outer_period_examples':periods}


def lattice_and_budget():
    lattice=0
    for p,ell,r in product(range(-4,5),repeat=3):
        coefficients=(p+ell,ell+r,-r)
        assert (coefficients[1]==0)==(ell+r==0)
        if ell+r==0:
            # (p+ell)e_x-r e_z = p gamma_U + ell(Gamma_H-Gamma_J).
            assert coefficients==(p+ell,0,ell)
        lattice+=1
    budgets=0
    for g in range(5):
        for b in range(1,6):
            chi=2-2*g-b
            areas=[F(i+1,50) for i in range(b)]
            for signs in product([-1,1],repeat=b):
                signed=sum(s*a for s,a in zip(signs,areas))
                chern=chi-sum(signs)
                area=signed-chern
                assert area+chi==sum(s*(a+1) for s,a in zip(signs,areas))
                if all(s==1 for s in signs):
                    assert area==sum(areas)+2*(g+b-1)
                budgets+=1
    thresholds=0
    for Ab in [F(1,5),F(1,2),F(4,5)]:
        A0=Ab*F(3,4)
        inner=A0*F(1,3)
        for eH,eJ in product(range(-5,6),repeat=2):
            bridge=inner+2+eH*(Ab+2)
            if eH<0: assert bridge<0
            else: assert bridge>=inner+2
            assert (inner+eH*Ab==bridge)==(eH==-1)
            thresholds+=1
    return {'extendable_lattice_checks':lattice,'signed_budget_checks':budgets,
            'natural_annulus_twist_checks':thresholds,
            'scope':'Integer and rational consistency; no assertion of a bridge embedding.'}


def membrane_checks():
    def vector(r,theta,x,y,e1=0.2):
        p=math.pi*(x*x+y*y)
        a=cmath.exp(2j*math.pi*theta)
        first=a*math.sqrt((r+e1+p)/math.pi)
        second=a.conjugate()*complex(x,y)
        return [first.real,first.imag,second.real,second.imag,0.0,0.0]
    error=0.0
    count=0
    for params in product([0.0,0.15],[0.0,0.13,0.39],[0.0,0.1],[0.0,0.08]):
        columns=[]
        for k in range(4):
            eps=1e-6
            lo=list(params); hi=list(params)
            lo[k]-=eps; hi[k]+=eps
            columns.append([(a-b)/(2*eps) for a,b in zip(vector(*hi),vector(*lo))])
        for i in range(4):
            for j in range(4):
                value=sum(columns[i][k]*columns[j][k+1]-columns[i][k+1]*columns[j][k]
                          for k in [0,2,4])
                expected=1 if (i,j) in [(0,1),(2,3)] else (-1 if (j,i) in [(0,1),(2,3)] else 0)
                error=max(error,abs(value-expected))
        count+=1
    assert error<2e-8
    return {'four_parameter_symplectic_samples':count,'maximum_pullback_error':error,
            'analytic_pullback':'d r wedge d theta + d Re(zeta) wedge d Im(zeta)',
            'scope':'A family of membranes in the exact chart; not a saturated corner fibration.'}


def main():
    report={'scope':'Finite consistency checks only; direct geometric proofs are in the manuscript.',
            'cutoff':differential_checks(),'lobe_actions_and_periods':lobe_checks(),
            'lattices_and_surface_budget':lattice_and_budget(),'membranes':membrane_checks(),
            'status':'passed'}
    path=Path(__file__).with_name('boundary_surgery_checks.json')
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print('Boundary surgery checks passed:',path.name)


if __name__=='__main__':
    main()
