## Feon
Feon is a python-based finite element analysis framework for education and research purpose by Dr. Pei Yaoyao at the sponsor of Hubei Univeristy of Technology!
## Version
Current:1.0.0
## Requirements
Numpy is a must. Matplotlib is needed for visualization. Mpmath is need for Matrix derivation.
## Installation
```
Using pip:
$pip install feon
```
Or
```
python install setup.py
```
## Packages
* sa---For structrual analysis
* ffa --- For fluid flow analysis
* derivation --- matrix derivation 

## Elements suported

* Spring1D11 
* Spring2D11
* Spring3D11

* Link1D11
* Link2D11
* Link3D11

* Beam1D11
* Beam2D11
* Beam3D11

* Tri2d11S---- Triange elements for plane stress problem
* Tri2D11 ---- Triange elements for plane strain problem
* Tetra3D11 
* Quad2D11S 
* Quad2D11
* Plate3D11 ---Midline plate
* Brick3D11

**We name the elements with "Name" + "dimensison" + 'order" + "type", type 1 means elastic.**

## Examples
**A Truss prombelem**
![image](https://github.com/YaoyaoBae/Feon/blob/master/examples/problems/truss/screenshot.png)
```python
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

```
**Embedded wall
![image](https://github.com/YaoyaoBae/Feon/blob/master/examples/problems/embedded%20wall/screenshot.png)
```python
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
        

```

