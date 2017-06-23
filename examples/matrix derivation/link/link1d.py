# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------


from feon.sa import *
from feon.derivation import *
if __name__ == "__main__":

    #create element
    e = Line()

    #create nodes
    nds = [Node(0,0),Node(1,0)]

    #compute shape function
    e.calc_nbase(nds)

    #set L
    e.set_L(shape = (1,1),mapping = {(0,0):1})

    #set material
    E = 210e6
    A = 0.005

    #set constitutive matrix
    e.set_D(D=E*A)

    #compute element matrix
    e.calc_ke()


    #comparation
    e1 = Link1D11(nds,E,A)
    e1.calc_ke()
    
    print "stiffness matrix of e is:\n%r)"%(e.ke,)
    print "stiffness matrix of e1 is:\n%r)"%(e1.ke,)

    #quad link element
    nds1 = [Node(0,0),Node(1,0),Node(0.5,0)]
    e2 = Line()
    e2.calc_nbase(nds1)
    e2.set_D(D=E*A)
    e2.set_L(shape = (1,1),mapping = {(0,0):1})
    e2.calc_ke()
    print "stiffness matrix of e is:\n%r)"%(e2.ke,)
    
