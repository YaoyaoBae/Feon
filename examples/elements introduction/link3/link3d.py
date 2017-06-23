# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from feon.sa import *

if __name__ == "__main__":
    
    E = 200e6
    A1 = 0.001
    A2 = 0.002

    n0 = Node(0,0,0)
    n1 = Node(0,-4,-5)
    n2 = Node(-3,0,-5)
    n3 = Node(0,4,-5)

    e0 = Link3D11((n0,n1),E,A1)
    e1 = Link3D11((n0,n2),E,A2)
    e2 = Link3D11((n0,n3),E,A1)

    s = System()
    s.add_nodes(n0,n1,n2,n3)
    s.add_elements(e0,e1,e2)
    s.add_node_force(0,Fx = 12)
    s.add_fixed_sup(1,2,3)
    s.solve()


    from matplotlib.ticker import FuncFormatter
    import matplotlib.pyplot as plt
    import numpy as np
    def stresses(x,pos):
        return "$%1.1fMPa$"%(x*1e-3)
    x = np.arange(3)
    stress = [abs(el.sx[0][0]) for el in [e0,e1,e2]]


    formatter = FuncFormatter(stresses)
    fig,ax = plt.subplots()
    ax.yaxis.set_major_formatter(formatter)
   
    plt.bar(x,stress,0.2,color = ["r","b","g"])
    ax.set_xticks(x+0.1)
    ax.set_xticklabels(("$Bar 0$","$Bar 1$","$Bar 2$"))
    ax.set_ylabel("$N/kN$")
    ax.set_xlim([-0.5,3])
    plt.show()
    
   
   

