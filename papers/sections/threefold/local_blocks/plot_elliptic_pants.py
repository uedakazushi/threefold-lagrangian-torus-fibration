#!/usr/bin/env python3
"""Scientific figure for the unmodified inner Morse model and its Reeb graph."""
from pathlib import Path
import io
import os
import numpy as np
os.environ.setdefault('MPLCONFIGDIR',str(Path('/tmp')/'general_type_matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                     'mathtext.fontset':'dejavusans','pdf.fonttype':42,
                     'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(10.8,4.7),constrained_layout=False)
gs=fig.add_gridspec(1,2,width_ratios=[1.45,1],left=.07,right=.98,bottom=.27,top=.86,wspace=.34)
ax=fig.add_subplot(gs[0,0])
x=np.linspace(-1.88,1.88,800)
y=np.linspace(-1.22,1.22,600)
X,Y=np.meshgrid(x,y)
H=(X*X-Y*Y-1)**2+4*X*X*Y*Y
ax.contourf(X,Y,H,levels=[0,.25,4],colors=['#e1e4e8','#f2f7fa'])
ax.contour(X,Y,H,levels=[.25],colors=['#237b78'],linewidths=2.0)
ax.contour(X,Y,H,levels=[1],colors=['#ac3939'],linewidths=2.1)
ax.contour(X,Y,H,levels=[4],colors=['#245779'],linewidths=2.0)
ax.scatter([0],[0],marker='x',s=45,c='#ac3939',linewidths=2,zorder=5)
ax.text(-1,0,'removed',ha='center',va='center',fontsize=8.5,color='#45515c')
ax.text(1,0,'removed',ha='center',va='center',fontsize=8.5,color='#45515c')
ax.annotate('saddle',xy=(0,0),xytext=(0,.64),ha='center',color='#ac3939',
            arrowprops={'arrowstyle':'-','color':'#ac3939','lw':1})
ax.set(xlabel=r'$\operatorname{Re}w$',ylabel=r'$\operatorname{Im}w$',
       xlim=(-1.88,1.88),ylim=(-1.22,1.22))
ax.set_aspect('equal')
ax.set_title(r'Inner model: $h(w)=|w^2-1|^2$',pad=12,fontsize=12)
ax.grid(color='#d7dfe4',alpha=.5,lw=.5)

bx=fig.add_subplot(gs[0,1])
bx.plot([0,0],[1,2.4],c='#245779',lw=3)
bx.plot([-.78,0],[.15,1],c='#237b78',lw=3)
bx.plot([.78,0],[.15,1],c='#237b78',lw=3)
bx.scatter([0],[1],s=50,c='#ac3939',zorder=5)
bx.scatter([0,-.78,.78],[2.4,.15,.15],s=45,facecolors='white',
           edgecolors=['#245779','#237b78','#237b78'],linewidths=2,zorder=5)
bx.text(.12,2.4,r'$h=4$: one $T^3$',va='center',fontsize=10)
bx.text(-.12,1.15,r'$h=1$',ha='right',va='center',color='#ac3939',fontsize=10)
bx.text(.12,1.15,r'$T^2\times(S^1\vee S^1)$',va='center',fontsize=10)
bx.text(-.78,-.06,r'$\gamma_{-}$',ha='center',va='top',fontsize=12)
bx.text(.78,-.06,r'$\gamma_{+}$',ha='center',va='top',fontsize=12)
bx.text(0,-.43,r'$h=1/4$: two regular $T^3$ ports',ha='center',fontsize=10)
bx.text(0,2.9,'Connected-component quotient',ha='center',fontsize=12)
bx.set(xlim=(-1.3,1.65),ylim=(-.58,3.12))
bx.axis('off')

fig.legend(handles=[Line2D([0],[0],color='#237b78',lw=2,label='new regular ports'),
                    Line2D([0],[0],color='#ac3939',lw=2,label='hyperbolic singular level'),
                    Line2D([0],[0],color='#245779',lw=2,label='intermediate regular level')],
           loc='lower center',bbox_to_anchor=(.5,.065),ncol=3,frameon=False,fontsize=9.5)
fig.text(.5,.018,'The prescribed outer collar lies beyond this inner model and is unchanged.',
         ha='center',fontsize=9,color='#495862')
for suffix in ['pdf','png']:
    data=io.BytesIO()
    fig.savefig(data,format=suffix,dpi=180,facecolor='white')
    path=ROOT/('elliptic_pants_model.'+suffix)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_bytes(data.getvalue())
    os.replace(temporary,path)
plt.close(fig)
print('Inner-model figure generated.')
