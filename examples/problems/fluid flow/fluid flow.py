# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from feon.ffa import *
from feon.tools import pair_wise
import numpy as np
if __name__ == "__main__":

    #permeability 
    Kxx = -2e-5

    #create nodes and elements
    A = np.pi*(np.linspace(0.06,0.15,7)[:-1]+0.0075)
    nds = [Node(-i*0.1,0) for i in xrange(7)]
    els = []
    for i in xrange(6):
        els.append(E1D((nds[i],nds[i+1]),Kxx,A[i]))

    #create FEA system

    s = System()
    s.add_nodes(nds)
    s.add_elements(els)
    
    s.add_node_head(0,0.2)
    s.add_node_head(6,0.1)
    s.solve()

    print [nd.head["H"] for nd in nds]
    print [el.velocity["Vx"] for el in els]
