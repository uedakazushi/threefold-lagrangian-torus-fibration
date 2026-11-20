#!/usr/bin/env python3
"""Selected numerical checks for the exact Hamiltonian wall straightening.

The existence and relative claims are proved in the accompanying manuscript.
This checks the nontrivial action inversion and Hamiltonian sign numerically.
Standard library only.
"""
from functools import lru_cache
from pathlib import Path
import json
import math


TAU=2*math.pi


def radial(h,j,u):
    return math.sqrt(max(0.0,(h-u)*u*(u-j)))/math.pi**1.5


def maximum_point(h,j):
    return (h+j+math.sqrt(h*h-h*j+j*j))/3


def lower_root(h,j,r):
    left,right=max(0.0,j),maximum_point(h,j)
    target=math.pi**3*r*r
    for _ in range(54):
        u=(left+right)/2
        if (h-u)*u*(u-j)<target:
            left=u
        else:
            right=u
    return (left+right)/2


def graph(h,j,c,q,phi):
    co,si=math.cos(TAU*phi),math.sin(TAU*phi)
    disc=math.sqrt(q-c*c*si*si)
    r=c*co+disc
    u=lower_root(h,j,r)
    D=(2*u-j)*(h-u)-u*(u-j)
    uq=math.pi**3*r/(D*disc)
    rphi=-TAU*c*si-TAU*c*c*si*co/disc
    uphi=2*math.pi**3*r*rphi/D
    return u,uq,uphi


@lru_cache(maxsize=8192)
def action(h,j,c,q,n=128):
    values=[graph(h,j,c,q,k/n) for k in range(n)]
    return (sum(x[0] for x in values)/n,
            sum(x[1] for x in values)/n)


def simpson(fun,left,right,n=128):
    if left==right:
        return 0.0
    dx=(right-left)/n
    total=fun(left)+fun(right)
    total+=4*sum(fun(left+dx*k) for k in range(1,n,2))
    total+=2*sum(fun(left+dx*k) for k in range(2,n,2))
    return total*dx/3


def generator(h,j,c,q,phi):
    # Integral from 0 to phi of U-I; its period is zero.
    I=action(h,j,c,q)[0]
    return simpson(lambda t:graph(h,j,c,q,t)[0]-I,0.0,phi)


def q_from_u(h,j,c,t,u,phi,qa,qb):
    q=(qa+qb)/2
    left,right=qa,qb
    for _ in range(14):
        U,Uq,_=graph(h,j,c,q,phi)
        I,Iq=action(h,j,c,q)
        error=(1-t)*U+t*I-u
        if abs(error)<2e-14:
            return q
        if error>0:
            right=q
        else:
            left=q
        step=q-error/((1-t)*Uq+t*Iq)
        q=step if left<step<right else (left+right)/2
    assert abs((1-t)*graph(h,j,c,q,phi)[0]+t*action(h,j,c,q)[0]-u)<1e-10
    return q


def K(h,j,c,t,u,phi,qa,qb):
    q=q_from_u(h,j,c,t,u,phi,qa,qb)
    return generator(h,j,c,q,phi)


def wall_example(h,j,c,e1):
    b=(h-e1)/2
    F=radial(h,j,maximum_point(h,j))
    fb=radial(h,j,b)
    assert 2*c<min(fb,F-fb)
    qa=(c*c+(fb-c)**2)/2
    qb=((fb+c)**2+(F-c)**2)/2
    assert action(h,j,c,qa)[0]<b<action(h,j,c,qb)[0]
    left,right=qa,qb
    for _ in range(45):
        q=(left+right)/2
        if action(h,j,c,q)[0]<b:
            left=q
        else:
            right=q
    q0=(left+right)/2
    I,Iq=action(h,j,c,q0)
    periodic_error=abs(generator(h,j,c,q0,1.0))
    normal_errors=[]
    positive=[]
    for t in (0.0,0.4,1.0):
        for phi in (0.17,0.43,0.79):
            U,Uq,Uphi=graph(h,j,c,q0,phi)
            u=(1-t)*U+t*I
            eps=2e-5
            Ku=(K(h,j,c,t,u+eps,phi,qa,qb)-K(h,j,c,t,u-eps,phi,qa,qb))/(2*eps)
            Kphi=(K(h,j,c,t,u,phi+eps,qa,qb)-K(h,j,c,t,u,phi-eps,qa,qb))/(2*eps)
            normal=-Kphi-(1-t)*Uphi*Ku
            normal_errors.append(abs(normal-(I-U)))
            positive.append((1-t)*Uq+t*Iq)
    assert abs(I-b)<1e-12
    assert periodic_error<1e-11
    assert max(normal_errors)<2e-7
    assert min(positive)>0
    # The full canonical coordinates are unimodular and have lambda=H dα+J dβ+u dφ.
    return dict(h=h,j=j,c=c,e1=e1,wall_u=b,q0=q0,
                toric_Q_value_at_wall=fb*fb+c*c,
                original_wall_Q_range=[(fb-c)**2,(fb+c)**2],
                action_at_q0=I,action_derivative=Iq,
                periodic_generator_error=periodic_error,
                normal_velocity_max_error=max(normal_errors),
                minimum_interpolated_q_derivative=min(positive))


def main():
    examples=[wall_example(*p) for p in
              [(3.0,0.0,0.02,1.0),(4.0,0.2,0.02,1.0),
               (2.0,-0.15,0.008,0.6)]]
    # One action annulus contains both charged walls even when their positions coincide.
    h,e1,e2,c=5.0,1.0,0.4,0.01
    walls=[]
    for j in (-0.1,0.0,0.2,0.4,0.5):
        b1=(h-e1)/2
        b13=(h+j-e1-e2)/2
        L=max(0.0,j)
        um=maximum_point(h,j)
        F=radial(h,j,um)
        for b in [b1,b13]:
            assert L<b<um
            f=radial(h,j,b)
            assert 2*c<min(f,F-f)
        if j==e2:
            assert b1==b13
        walls.append(dict(j=j,b1=b1,b13=b13))
    result={'scope':'Numerical checks of action inversion, periodic Hamiltonian generator and normal velocity; not a global fibration certificate.',
            'exact_straightening_examples':examples,
            'normal_velocity_samples':27,
            'two_wall_examples':walls,
            'max_normal_velocity_error':max(x['normal_velocity_max_error'] for x in examples),
            'status':'passed'}
    path=Path(__file__).with_name('straightening_checks.json')
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
