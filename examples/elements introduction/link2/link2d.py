# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------


from feon.sa import *
if __name__ == "__main__":
    E = 210e6
    A = 0.005
    K = 79e3
    n0 = Node(0,0)
    n1 = Node(5,0)
    n2 = Node(10,0)
    n3 = Node(15,0)
    n4 = Node(5,7)
    n5 = Node(10,7)
    n6 = Node(15,-1)

    e0 = Link2D11((n0,n1),E,A)
    e1 = Link2D11((n1,n2),E,A)
    e2 = Link2D11((n2,n3),E,A)
    e3 = Link2D11((n4,n0),E,A)
    e4 = Link2D11((n4,n1),E,A)
    e5 = Link2D11((n4,n2),E,A)
    e6 = Link2D11((n4,n5),E,A)
    e7 = Link2D11((n5,n2),E,A)
    e8 = Link2D11((n5,n3),E,A)
    e9 = Spring2D11((n3,n6),K)
    s = System()
    s.add_nodes(n0,n1,n2,n3,n4,n5,n6)
    s.add_elements(e0,e1,e2,e3,e4,e5,e6,e7,e8,e9)
    s.add_node_force(4,Fx = 30)
    s.add_fixed_sup(0,6)
    s.solve()

