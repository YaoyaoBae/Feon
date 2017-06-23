# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from feon.sa import *

if __name__ == "__main__":
    k = 200
    n0 = Node(0,0)
    n1 = Node(0,-3)
    n2 = Node(0,3)
    n3 = Node(3,0)
    n4 = Node(6,0)
    e0 = Spring2D11((n0,n3),k)
    e1 = Spring2D11((n1,n3),k)
    e2 = Spring2D11((n2,n3),k)
    e3 = Spring1D11((n3,n4),k)
    s = System()
    s.add_nodes(n0,n1,n2,n3,n4)
    s.add_elements(e0,e1,e2,e3)
    s.add_node_force(n4.ID,Fx = 5)
    s.add_fixed_sup(n0.ID,1,2)
    s.solve()


    from matplotlib.lines import Line2D
    import matplotlib.pyplot as plt
    from feon.sa.draw2d import *
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim([-3,9])
    ax.set_ylim([-4,4])
    for el in [e0,e1,e2,e3]:
        draw_element(ax,el,marker = "o",lw = 4,color = "g")
        draw_element_disp(ax,el,factor = 0.05,ms = 4)
        draw_element_ID(ax,el,dx = 0.2,dy = 0.2,color = "r")
    draw_fixed_sup(ax,n0,factor = (0.4,0.4))
    draw_fixed_sup(ax,n1,factor = (0.4,0.4))
    draw_fixed_sup(ax,n2,factor = (0.4,0.4))
    for nd in s.get_nodes():
        draw_node_ID(ax,nd,dx = 0.2,dy = 0.1,color = "b")
    plt.show()

    
    
    
