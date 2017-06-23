# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from feon.sa import *
from feon.tools import pair_wise
from feon.sa.draw2d import *
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
if __name__ == "__main__":

    #material property
    E = 210e6
    A1 = 31.2e-2
    A2 = 8.16e-2

    #create nodes and elements
    nds1 = []
    nds2 = []
    for i in xrange(13):
        nds1.append(Node(i,0))
    for i in xrange(11):
        nds2.append(Node(i+1,-1))
    els = []
    for e in pair_wise(nds1):
        els.append(Link2D11((e[0],e[1]),E,A1))
    for e in pair_wise(nds2):
        els.append(Link2D11((e[0],e[1]),E,A1))

    for i in xrange(6):
        els.append(Link2D11((nds1[i],nds2[i]),E,A2))
    for i in xrange(6):
        els.append(Link2D11((nds2[i+5],nds1[i+7]),E,A2))

    for i in xrange(11):
        els.append(Link2D11((nds1[i+1],nds2[i]),E,A2))

    #create system
    s = System()

    
    #add nodes and elements in to the system
    s.add_nodes(nds1,nds2)
    s.add_elements(els)

    #apply boundry condition
    s.add_node_force(nds1[0].ID,Fy = -1000)
    s.add_node_force(nds1[-1].ID,Fy = -1000)
    for i in xrange(1,12):
        s.add_node_force(nds1[i].ID,Fy = -1900)
    s.add_fixed_sup(nds1[0].ID)
    s.add_rolled_sup(nds1[-1].ID,"y")

    #solve 
    s.solve()

    #show results
    disp = [np.sqrt(nd.disp["Ux"]**2+nd.disp["Uy"]**2) for nd in s.get_nodes()]
    
    eforce = [el.force["N"][0][0] for el in s.get_elements()]
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.yaxis.get_major_formatter().set_powerlimits((0,1)) 
    ax2 = fig.add_subplot(212)
    ax2.yaxis.get_major_formatter().set_powerlimits((0,1)) 
    ax.set_xlabel(r"$Node ID$")
    ax.set_ylabel(r"$Disp/m$")
    ax.set_ylim([-0.05,0.05])
    ax.set_xlim([-1,27])
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.plot(range(len(disp)),disp,"r*-")
    ax2.set_xlabel(r"$Element ID$")
    ax2.set_xlim([-1,46])
    ax2.set_ylabel(r"$N/kN$")
    ax2.set_ylim(-40000,40000)
    ax2.xaxis.set_minor_locator(MultipleLocator(1))
    for i in xrange(len(eforce)):
        ax2.plot([i-0.5,i+0.5],[eforce[i],eforce[i]],"ks-",ms = 3)
    plt.show()
    draw_bar_info(els[5])
    
