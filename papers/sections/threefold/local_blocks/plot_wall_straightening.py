#!/usr/bin/env python3
"""Regenerate the two panels illustrating exact relative wall straightening.

Requires matplotlib.  The construction is proved in wall_extension.tex;
this plot is an illustration, not an existence certificate.
"""
from pathlib import Path
import math
import io
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'general-type-fibration-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from verify_straightening import action, graph, radial, maximum_point


def main():
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False,
                         'axes.spines.right':False, 'pdf.fonttype':42})
    h,j,c,e1,e2=3.0,0.0,0.02,1.0,0.4
    b=(h-e1)/2
    fb=radial(h,j,b)
    F=radial(h,j,maximum_point(h,j))
    lo=(c*c+(fb-c)**2)/2
    hi=((fb+c)**2+(F-c)**2)/2
    for _ in range(45):
        mid=(lo+hi)/2
        if action(h,j,c,mid)[0]<b:
            lo=mid
        else:
            hi=mid
    q=(lo+hi)/2
    phase=[k/400 for k in range(401)]
    u=[graph(h,j,c,q,p)[0] for p in phase]
    I=action(h,j,c,q)[0]
    fig,axs=plt.subplots(1,2,figsize=(10.4,3.8),layout='constrained')
    ax=axs[0]
    for t,col,style in [(0,'#b96b26','-'),(.5,'#748a93','--'),(1,'#1d6276','-')]:
        ax.plot(phase,[(1-t)*v+t*I for v in u],color=col,ls=style,
                lw=2 if t in (0,1) else 1.7,label=fr'$t={t:g}$')
    ax.axhline(b,ls=':',color='#17374b',lw=.8,zorder=0)
    ax.set(xlabel=r'period-one angle $\phi$',ylabel=r'$u=P_y$',
           title=r'(a) $U_t=(1-t)U+tI$: equal action')
    ax.set_xticks([0,.25,.5,.75,1])
    ax.legend(loc='upper right',frameon=False,ncol=3)
    ax.text(.04,.06,r'$I(q)=b_1=1$'+'\n'+fr'$q={q:.8f}$',
            transform=ax.transAxes,color='#17374b')
    ax.set_ylim(min(u)-.04,max(u)+.065)
    ax.grid(axis='y',alpha=.18)

    ax=axs[1]
    js=[-.1+.7*k/400 for k in range(401)]
    b1=[(h-e1)/2 for v in js]
    b13=[(h+v-e1-e2)/2 for v in js]
    delta=.035
    mids=[(a+b)/2 for a,b in zip(b1,b13)]
    half=[math.sqrt(((a-b)/2)**2+delta**2) for a,b in zip(b1,b13)]
    bottom=[a-b for a,b in zip(mids,half)]
    top=[a+b for a,b in zip(mids,half)]
    # The entire illustrated enclosing interval satisfies the theorem's gaps.
    for v,lower,upper in zip(js,bottom,top):
        peak=radial(h,v,maximum_point(h,v))
        assert max(0,v)<lower<upper<maximum_point(h,v)
        assert radial(h,v,lower)>2*c
        assert peak-radial(h,v,upper)>2*c
    ax.fill_between(js,bottom,top,color='#d6e7eb',label='one smooth action band')
    ax.plot(js,b1,color='#b96b26',lw=2,label=r'$b_1=(h-e_1)/2$')
    ax.plot(js,b13,color='#1d6276',lw=2,label=r'$b_{13}=(h+j-e_1-e_2)/2$')
    ax.axvline(e2,color='#657a84',lw=1,ls=':')
    ax.scatter([e2],[b],color='#17374b',s=22,zorder=5)
    ax.annotate(r'$j=e_2$',(e2,b),xytext=(.455,.86),
                arrowprops={'arrowstyle':'-','color':'#657a84'},color='#17374b')
    ax.set(xlabel=r'moment $j$',ylabel=r'wall position $u$',
           title='(b) Both charged walls in one band')
    ax.set_ylim(.67,1.27)
    ax.legend(loc='upper left',frameon=False,fontsize=8.8)
    ax.grid(axis='y',alpha=.18)
    path=Path(__file__).with_name('wall_straightening_geometry.pdf')
    buffer=io.BytesIO()
    fig.savefig(buffer,format='pdf',metadata={'Title':'Exact action straightening of regular walls',
                                            'Author':'Research note', 'CreationDate':None})
    data=buffer.getvalue()
    assert data.rstrip().endswith(b'%%EOF')
    stage=path.with_suffix('.pdf.tmp')
    stage.write_bytes(data)
    os.replace(stage,path)
    plt.close(fig)
    print(path)


if __name__=='__main__':
    main()
