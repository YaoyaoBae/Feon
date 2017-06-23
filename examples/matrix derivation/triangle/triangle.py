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
    e = Triangle()

    #create nodes
    nds = [Node(0,0),Node(1,0),Node(1,1)]

    #compute shape function
    e.calc_nbase(nds)

    #set L
    e.set_L(shape=(3,2),mapping={(0,0):(1,0),
                                  (1,1):(0,1),
                                  (2,0):(0,1),
                                  (2,1):(1,0)})


    #set material,for stress problem
    E = 210e6
    nu = 0.3
    t = 0.5
    a = t*E/(1-nu**2)
    D = a*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2.]])

    #set constitutive matrix
    e.set_D(D=D)

    #compute element matrix
    e.calc_ke()

    #comparation
    e1 = Tri2D11S(nds,E,nu,t)
    e1.calc_Ke()
    
    print "stiffness matrix of e is:\n%r)"%(e.ke,)
    print "stiffness matrix of e1 is:\n%r)"%(e1.Ke,)
    
