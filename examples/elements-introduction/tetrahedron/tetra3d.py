# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------
from feon.sa import *
if __name__ == "__main__":
    E = 210e6
    nu = 0.3

    n0 = Node(0,0,0)
    n1 = Node(0.025,0,0)
    n2 = Node(0,0.5,0)
    n3 = Node(0.025,0.5,0)
    n4 = Node(0,0,0.25)
    n5 = Node(0.025,0,0.25)
    n6 = Node(0,0.5,0.25)
    n7 = Node(0.025,0.5,0.25)

    e0 = Tetra3D11((n0,n1,n3,n5),E,nu)
    e1 = Tetra3D11((n0,n3,n2,n6),E,nu)
    e2 = Tetra3D11((n5,n4,n6,n0),E,nu)
    e3 = Tetra3D11((n5,n6,n7,n3),E,nu)
    e4 = Tetra3D11((n0,n5,n3,n6),E,nu)

    s = System()
    s.add_nodes(n0,n1,n2,n3,n4,n5,n6,n7)
    s.add_elements(e0,e1,e2,e3,e4)

    s.add_node_force(2,Fy = 3.125)
    s.add_node_force(7,Fy = 3.125)
    s.add_node_force(3,Fy = 6.25)
    s.add_node_force(6,Fy = 6.25)

    s.add_fixed_sup(0,1,4,5)
    
    s.solve()
