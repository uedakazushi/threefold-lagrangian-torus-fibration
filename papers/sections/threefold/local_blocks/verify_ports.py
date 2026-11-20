#!/usr/bin/env python3
"""Local action checks for residual ports of W=bar(x)*y*z.

Independent comparisons:
* the raw upper action is a polygonal line integral of u_+ d(arg W),
  while the smooth cap action uses the stable identity p/R^2;
* the merged action is computed by its boundary integral and separately
  by integrating angular slices over the reduced sphere;
* the elliptic endpoint slope is derived from the local quadratic form.
This is numerical support for local formulas, not a global existence proof.
"""
from pathlib import Path
from math import pi, sqrt, cos, sin, acos, atan2
import json

N=4096

def simpson(fn,left,right,n=N):
    step=(right-left)/n
    total=fn(left)+fn(right)
    total+=4*sum(fn(left+(2*k-1)*step) for k in range(1,n//2+1))
    total+=2*sum(fn(left+2*k*step) for k in range(1,n//2))
    return total*step/3

def peak(h,j):
    u=(h+j+sqrt(h*h-h*j+j*j))/3
    F=sqrt((h-u)*u*(u-j))/pi**1.5
    f_second=-F*(1/(h-u)**2+1/u**2+1/(u-j)**2)/2
    return u,F,f_second

def roots(h,j,r):
    u0,F,_=peak(h,j)
    if r>=F*(1-2e-15):
        return u0,u0
    if r<1e-16:
        return max(0.0,j),h
    target=pi**3*r*r
    lo=max(0.0,j)
    hi=u0
    for _ in range(65):
        mid=(lo+hi)/2
        if (h-mid)*mid*(mid-j)<target:
            lo=mid
        else:
            hi=mid
    lower=(lo+hi)/2
    # Solve the small variable p=h-u_+ directly, avoiding cancellation.
    lo=0.0
    hi=h-u0
    for _ in range(65):
        mid=(lo+hi)/2
        if mid*(h-mid)*(h-j-mid)<target:
            lo=mid
        else:
            hi=mid
    p=(lo+hi)/2
    return lower,h-p

def upper_p(h,j,r):
    u0,F,_=peak(h,j)
    if r>=F*(1-2e-15):
        return h-u0
    target=pi**3*r*r
    lo=0.0
    hi=h-u0
    for _ in range(65):
        mid=(lo+hi)/2
        if mid*(h-mid)*(h-j-mid)<target:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2

def upper_action(h,j,c,q):
    radius=sqrt(q)
    def integrand(t):
        wx=c+radius*cos(t)
        wy=radius*sin(t)
        r=sqrt(wx*wx+wy*wy)
        p=upper_p(h,j,r)
        # p/r^2 = pi^3/[(h-p)(h-p-j)], also at r=0 by continuity.
        return pi**3*(q+c*radius*cos(t))/((h-p)*(h-j-p))
    return simpson(integrand,0,2*pi)/(2*pi)

def raw_upper_polygon(h,j,c,q,n=32768):
    radius=sqrt(q)
    total=0.0
    angle_total=0.0
    for k in range(n):
        t0=2*pi*k/n
        t1=2*pi*(k+1)/n
        tm=(t0+t1)/2
        w0=complex(c+radius*cos(t0),radius*sin(t0))
        w1=complex(c+radius*cos(t1),radius*sin(t1))
        product=w1*w0.conjugate()
        delta=atan2(product.imag,product.real)
        rm=sqrt((c+radius*cos(tm))**2+(radius*sin(tm))**2)
        p=upper_p(h,j,rm)
        total+=(h-p)*delta
        angle_total+=delta
    return total/(2*pi),angle_total/(2*pi)

def merged_boundary(h,j,c,q):
    _,F,_=peak(h,j)
    value=(F*F+c*c-q)/(2*c*F)
    phi0=acos(max(-1.0,min(1.0,value)))
    width=pi-phi0
    # phi=phi0+(pi-phi0)*sin(t)^2 regularizes the fold endpoint.
    def integrand(t):
        if t<=0 or t>=pi/2:
            return 0.0
        phi=phi0+width*sin(t)**2
        r=c*cos(phi)+sqrt(max(0.0,q-c*c*sin(phi)**2))
        lower,upper=roots(h,j,r)
        return (upper-lower)*2*width*sin(t)*cos(t)
    return simpson(integrand,0,pi/2)/pi

def merged_vertical(h,j,c,q):
    u0,F,_=peak(h,j)
    rmin=sqrt(q)-c
    ulo,uhi=roots(h,j,rmin)
    def angular_fraction(u):
        r=sqrt(max(0.0,(h-u)*u*(u-j)))/pi**1.5
        if r<=0:
            return 0.0
        value=(r*r+c*c-q)/(2*c*r)
        return 1-acos(max(-1.0,min(1.0,value)))/pi
    def piece(left,right):
        width=right-left
        def integrand(t):
            if t<=0 or t>=pi/2:
                return 0.0
            u=left+width*sin(t)**2
            return angular_fraction(u)*2*width*sin(t)*cos(t)
        return simpson(integrand,0,pi/2)
    return piece(ulo,u0)+piece(u0,uhi)

def lower_saddle_action(h,j,c):
    _,F,_=peak(h,j)
    q=(F-c)**2
    def integrand(phi):
        r=c*cos(phi)+sqrt(max(0.0,q-c*c*sin(phi)**2))
        return roots(h,j,r)[0]
    return simpson(integrand,0,2*pi)/(2*pi)

cases=[]
for h,j in [(pi,0),(2*pi,0),(pi,0.1*pi),(2*pi,-0.08*pi)]:
    _,F,f_second=peak(h,j)
    c=0.12*F
    qs=(F-c)**2
    qM=(F+c)**2
    qstar=(c*c+qs)/2
    qplus=(qs+qM)/2
    Bport=upper_action(h,j,c,qstar)
    Cboundary=merged_boundary(h,j,c,qplus)
    Cvertical=merged_vertical(h,j,c,qplus)
    assert abs(Cboundary-Cvertical)<1e-8
    Ls=lower_saddle_action(h,j,c)
    Bs=upper_action(h,j,c,qs)
    Cs=merged_vertical(h,j,c,qs)
    balance_error=abs(Ls+Bs+Cs-h)
    assert balance_error<2e-8
    if j==0:
        assert 0<Bport+Cboundary<h
    slope=1/(2*sqrt(c*F*(F+c)*(-f_second)))
    asymptotic=[]
    for ratio in [0.001,0.0005,0.00025]:
        epsilon=qM*ratio
        numerical=merged_boundary(h,j,c,qM-epsilon)/epsilon
        relative_error=abs(numerical-slope)/slope
        assert relative_error<0.003
        asymptotic.append({'epsilon_over_qM':ratio,'numerical_C_over_epsilon':numerical,'predicted_slope':slope,'relative_error':relative_error})
    assert asymptotic[-1]['relative_error']<asymptotic[0]['relative_error']/3
    cases.append({'h':h,'j':j,'c':c,'upper_port_action':Bport,'merged_port_action_boundary':Cboundary,'merged_port_action_vertical_area':Cvertical,'independent_merged_area_error':abs(Cboundary-Cvertical),'saddle_actions':{'lower':Ls,'upper_cap':Bs,'merged_cap':Cs,'sum_minus_h':Ls+Bs+Cs-h},'port_sum_over_h':(Bport+Cboundary)/h,'max_endpoint_asymptotics':asymptotic})

h=pi
j=0.07*pi
_,F,_=peak(h,j)
c=0.1*F
jumps=[]
for rel in [0.02,0.01,0.005]:
    qleft=c*c*(1-rel)
    qright=c*c*(1+rel)
    rawleft,windleft=raw_upper_polygon(h,j,c,qleft)
    rawright,windright=raw_upper_polygon(h,j,c,qright)
    Bleft=upper_action(h,j,c,qleft)
    Bright=upper_action(h,j,c,qright)
    errleft=abs(rawleft+Bleft)
    errright=abs(rawright+Bright-h)
    assert abs(windleft)<1e-10 and abs(windright-1)<1e-10
    assert max(errleft,errright)<3e-9
    assert abs((rawright-rawleft)-h)<0.001*h
    jumps.append({'relative_distance_from_c_squared':rel,'raw_I_left':rawleft,'raw_I_right':rawright,'raw_jump':rawright-rawleft,'expected_jump_limit_h':h,'upper_B_left':Bleft,'upper_B_right':Bright,'polygon_winding_left':windleft,'polygon_winding_right':windright,'transition_error_left':errleft,'transition_error_right':errright})
Bcritical=upper_action(h,j,c,c*c)

report={
 'scope':'Local residual-port actions, explicit cap geometry, and marked-period checks only; not a global boundary-free fibration proof.',
 'methods':{'upper_raw':'Midpoint density with exact increments of polygonal argument, 32768 edges.','upper_cap':'Simpson integration using p/R^2=pi^3/[(h-p)(h-p-j)], 4096 intervals.','merged_boundary':'Root-difference integral with a sin-squared endpoint substitution.','merged_area':'Independent angular-slice integral over u, split at the radial maximum.','max_endpoint':'Comparison with the reduced quadratic ellipse area.'},
 'upper_transition':{'h':h,'j':j,'c':c,'smooth_B_at_q_equals_c_squared':Bcritical,'examples':jumps},
 'parameter_examples':cases,
 'status':'passed'
}
output=Path(__file__).with_name('ports_action_checks.json')
output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':'passed','output':str(output),'parameter_examples':len(cases),'raw_jump_examples':len(jumps),'max_independent_area_error':max(x['independent_merged_area_error'] for x in cases),'max_saddle_balance_error':max(abs(x['saddle_actions']['sum_minus_h']) for x in cases),'max_upper_transition_error':max(max(x['transition_error_left'],x['transition_error_right']) for x in jumps)},indent=2))
