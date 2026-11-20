#!/usr/bin/env python3
"""Independent numerical checks of the positive-winding lower-branch action.
These checks support explicit formulas, not the global fibration theorem.
"""
from fractions import Fraction
from math import pi, cos, sqrt
import json

def simpson(f, a, b, n=4096):
    step=(b-a)/n
    val=f(a)+f(b)
    val+=4*sum(f(a+(2*k-1)*step) for k in range(1,n//2+1))
    val+=2*sum(f(a+2*k*step) for k in range(1,n//2))
    return val*step/3

def lower_root(H,j,R):
    lo=max(j,0.0)
    if R<1e-16:
        return lo
    hi=H/2
    target=pi**3*R*R
    assert lo*(lo-j)*(H-lo)<=target
    assert hi*(hi-j)*(H-hi)>target
    for _ in range(90):
        u=(lo+hi)/2
        if u*(u-j)*(H-u)<target:
            lo=u
        else:
            hi=u
    return (lo+hi)/2

def action(H,c):
    return simpson(lambda ph:lower_root(H,0,2*c*max(0,cos(ph))),-pi/2,pi/2)/(2*pi)

positivity=0
for H in [Fraction(1),Fraction(2),Fraction(4)]:
    for u in [H/Fraction(16),H/Fraction(8),H/Fraction(4),3*H/Fraction(8)]:
        for j in [-u,0,u/2,u]:
            for s in [-1,Fraction(-1,2),0,Fraction(1,2),1]:
                D=(2*u-j)*(H+s*u)+s*u*(u-j)
                assert D>0
                positivity+=1

examples=[]
for H,c in [(pi,0.01),(2*pi,0.025),(3*pi,0.05),(5*pi,0.08)]:
    actual=action(H,c)
    leading=2*c*sqrt(pi/H)
    quadratic=pi**3*c*c/(2*H*H)
    cubic=10*pi**3.5*c**3/(3*H**3.5)
    half=action(H,c/2)
    residual=actual-leading-quadratic-cubic
    residual_half=half-leading/2-quadratic/4-cubic/8
    assert actual>leading
    assert 14<residual/residual_half<19
    examples.append({'H_over_pi':H/pi,'c':c,'actual_action':actual,'frozen_action':leading,'quadratic_correction':quadratic,'cubic_correction':cubic,'remainder':residual,'remainder_ratio_when_halving_c':residual/residual_half})

H=2*pi
c=0.03
j=0.2*c
q=0.4*c*c
radius=sqrt(q)
def integrand(angle,derivative=False,Hvalue=H):
    R=sqrt(c*c+q+2*c*radius*cos(angle))
    u=lower_root(Hvalue,j,R)
    arg_derivative=(q+c*radius*cos(angle))/(R*R)
    if derivative:
        D=(2*u-j)*(Hvalue-u)-u*(u-j)
        u=-u*(u-j)/D
    return u*arg_derivative
analytic=simpson(lambda ang:integrand(ang,True),0,2*pi)/(2*pi)
step=H*1e-5
plus=simpson(lambda ang:integrand(ang,Hvalue=H+step),0,2*pi)/(2*pi)
minus=simpson(lambda ang:integrand(ang,Hvalue=H-step),0,2*pi)/(2*pi)
finite=(plus-minus)/(2*step)
assert abs(analytic-finite)<1e-10
report={'scope':'Explicit local lower-branch formulas only; not a global fibration proof.','radial_derivative_exact_rational_cases':positivity,'nodal_action_examples':examples,'return_derivative':{'analytic':analytic,'centered_difference':finite,'absolute_error':abs(analytic-finite)}}
from pathlib import Path
Path(__file__).with_name("positive_end_checks.json").write_text(
    json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(json.dumps(report,indent=2))
