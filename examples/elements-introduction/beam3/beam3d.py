# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from feon.sa import *
from feon.tools import pair_wise
if __name__ == "__main__":

    E = 210e6
    G = 84e6
    A = 0.02
    I =[5e-5,10e-5,20e-5]

    n0 = Node(0,0,0)
    n1 = Node(0,4,0)
    n2 = Node(4,4,0)
    n3 = Node(0,4,0)

    n4 = Node(0,0,5)
    n5 = Node(0,4,5)
    n6 = Node(4,4,5)
    n7 = Node(0,4,5)

    n8 = Node(1,0,5)
    n9 = Node(3,0,5)

    nds1 = [n0,n3,n2,n1]
    nds2 = [n4,n7,n6,n5]
    nds3 = [n4,n8,n9,n7,n6,n5]
    els = []
    for nd in pair_wise(nds3,True):
       els.append(Beam3D11(nd,E,G,A,I))

    for i in xrange(4):
       els.append(Beam3D11((nds1[i],nds2[i]),E,G,A,I))


    s = System()
    s.add_nodes(nds1,nds3)
    s.add_elements(els)

    s.add_node_force(nds2[2].ID,Fx = -15)
    s.add_element_load(els[1].ID,"q",(0,-5))
    s.add_element_load(els[4].ID,"tri",(-10,0))

    s.add_fixed_sup([nd.ID for nd in nds1])
    s.solve()
   



    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.ticker import MultipleLocator
    noe = s.noe
    non = s.non

    N = np.array([el.force["N"] for el in els])
    N1 = N[:,0]
    N2 = N[:,1]
   
    Ty = np.array([el.force["Ty"] for el in els])
    Ty1 = Ty[:,0]
    Ty2 = Ty[:,1]

    Tz = np.array([el.force["Tz"] for el in els])
    Tz1 = Tz[:,0]
    Tz2 = Tz[:,1]

    Mx = np.array([el.force["Mx"] for el in els])
    Mx1 = Mx[:,0]
    Mx2 = Mx[:,1]

    My = np.array([el.force["My"] for el in els])
    My1 = My[:,0]
    My2 = My[:,1]

    Mz = np.array([el.force["Mz"] for el in els])
    Mz1 = Mz[:,0]
    Mz2 = Mz[:,1]


   
    fig1,fig2 = plt.figure(),plt.figure()
   
    ax1 = fig1.add_subplot(311)
    ax2 = fig1.add_subplot(312)
    ax3 = fig1.add_subplot(313)

    ax4 = fig2.add_subplot(311)
    ax5 = fig2.add_subplot(312)
    ax6 = fig2.add_subplot(313)
    
    ax3.set_xticks([-1,noe+1],1)
    ax3.set_xlabel(r"$Element ID$")
    ax1.set_ylabel(r"$N/kN$")
    ax2.set_ylabel(r"$Ty/kN$")
    ax3.set_ylabel(r"$Tz/kN$")
    ax3.xaxis.set_major_locator(MultipleLocator(1))
    ax1.xaxis.set_major_locator(MultipleLocator(1))
    ax2.xaxis.set_major_locator(MultipleLocator(1))
    ax2.set_ylim([-7,7])
   
    for i in xrange(noe):
       ax1.plot([i-0.5,i+0.5],[N1[i],N2[i]],"gs-")
       ax2.plot([i-0.5,i+0.5],[Ty1[i],Ty1[i]],"rs-")
       ax3.plot([i-0.5,i+0.5],[Tz1[i],Tz1[i]],"ks-")
       
    ax6.set_xticks([-1,noe+1],1)
    ax6.set_xlabel(r"$Element ID$")
    ax4.set_ylabel(r"$Mx/kNm$")
    ax5.set_ylabel(r"$My/kNm$")
    ax6.set_ylabel(r"$Mz/kNm$")
    ax6.xaxis.set_major_locator(MultipleLocator(1))
    ax4.xaxis.set_major_locator(MultipleLocator(1))
    ax5.xaxis.set_major_locator(MultipleLocator(1))
   
    for i in xrange(noe):
       ax4.plot([i-0.5,i+0.5],[Mx1[i],Mx2[i]],"gs-")
       ax5.plot([i-0.5,i+0.5],[My1[i],My2[i]],"rs-")
       ax6.plot([i-0.5,i+0.5],[Mz1[i],Mz2[i]],"ks-")

    plt.show()
   
                          
      

   
   
   

   
   
   
   
      

