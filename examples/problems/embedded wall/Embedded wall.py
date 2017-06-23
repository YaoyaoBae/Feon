# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from feon.sa import *
from feon.tools import pair_wise
import matplotlib.pyplot as plt
from feon.sa.draw2d import *
if __name__ == "__main__":
    E1 = 2.85e6
    E2 = 200e6
    k = 15000
    I = 0.0427
    A = 0.8
    A1 = 0.003
    ka = 0.6

    nds1 =[Node(0,-i) for i in xrange(10)]
    nds2 = [Node(0,-(i+20)*0.5) for i in xrange(11)]
    nds3 = [Node(-0.5,-(i+20)*0.5) for i in xrange(11)]
    nds4 = [Node(-1.5,-2),Node(-1.5,-6)]

    els=[]
    for nd in pair_wise(nds1+nds2):
        els.append(Beam2D11(nd,E1,A,I))

    for i in xrange(11):
        els.append(Spring2D11((nds2[i],nds3[i]),k))

    els.append(Link2D11((nds4[0],nds1[2]),E2,A1))
    els.append(Link2D11((nds4[1],nds1[6]),E2,A1))
        
    s = System()
    s.add_nodes(nds1,nds2,nds3,nds4)
    s.add_elements(els)

    nid1 = [nd.ID for nd in nds3]
    nid2 = [nd.ID for nd in nds4]
    s.add_fixed_sup(nid1,nid2)
    for i,el in enumerate(els[:10]):
        s.add_element_load(el.ID,"tri",-18*ka)
        s.add_element_load(el.ID,"q",-i*18*ka)
        
    for el in els[10:20]:
        s.add_element_load(el.ID,"q",-180*ka)

    for nd in nds1:
        nd.set_disp(Uy =0)

    for nd in nds2:
        nd.set_disp(Uy = 0)

    s.solve()


    disp = np.array([nd.disp["Ux"] for nd in nds1]+[nd.disp["Ux"] for nd in nds2])*1000
    Mz = [el.force["Mz"][0][0] for el in els[:20]]

    fig1,fig2,fig3 = plt.figure(),plt.figure(),plt.figure()
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)
    ax3 = fig3.add_subplot(111)
    
    Y1 = [-i for i in xrange(10)]+[-(i+20)*0.5 for i in xrange(11)]
    Y2 = [-i-0.5 for i in xrange(10)]+[-(i+20)*0.5-0.5 for i in xrange(10)]
    ax1.plot(disp,Y1,"r--")
    ax1.set_xlabel("$Ux/mm$")
    ax1.set_ylabel("$Height/m$")
    ax2.plot(Mz,Y2,"r-+")
    ax2.set_xlabel("$Mz/kN.m$")
    ax2.set_ylabel("$Height/m$")

    for el in els[:20]:
        draw_element(ax3,el,lw = 10,color = "g")
        
    for el in els[20:31]:
        draw_spring(ax3,el,color = "k")

    for el in els[31:]:
        draw_element(ax3,el,lw = 1.5,color = "k",marker = "s")

    for nd in nds3+nds4:
        draw_fixed_sup(ax3,nd,factor = (0.4,4),color ="k")
    
    ax3.set_xlim([-2,2])
    ax3.set_ylim([-16,1])
    plt.show()
        
