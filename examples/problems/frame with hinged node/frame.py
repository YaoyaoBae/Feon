# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from feon.sa import *
from feon.tools import pair_wise

#define beamlink element
class BeamLink2D11(StructElement):
    def __init__(self,nodes,E,A,I):
        StructElement.__init__(self,nodes)
        self.E = E
        self.A = A
        self.I = I

    #define node degree of freedom, left node has three dofs
    #while the right node has only two
    def init_unknowns(self):
        self.nodes[0].init_unknowns("Ux","Uy","Phz")
        self.nodes[1].init_unknowns("Ux","Uy")
        self._ndof = 3

    #transformative matrix
    def calc_T(self):
        TBase = _calc_Tbase_for_2d_beam(self.nodes)
        self._T = np.zeros((6,6))
        self._T[:3,:3] = self._T[3:,3:] = TBase

    #stiffness matrix
    def calc_ke(self):
        self._ke = _calc_ke_for_2d_beamlink(E = self.E,A = self.A,I = self.I,L = self.volume)

def _calc_ke_for_2d_beamlink(E = 1.0,A = 1.0,I = 1.0,L = 1.0):
    a00 = E*A/L
    a03 = -a00
    a11 = 3.*E*I/L**3
    a12 = 3.*E*I/L**2
    a14 = -a11
    a22 = 3.*E*I/L
    T = np.array([[a00,  0.,   0.,  a03,  0.,0.],
                  [ 0., a11,  a12,  0., a14, 0.],
                  [ 0., a12,  a22,  0.,-a12, 0.],
                  [a03,  0.,   0., a00,  0., 0.],
                  [ 0., a14, -a12,  0., a11, 0.],
                  [ 0.,  0.,    0.,  0., 0., 0.]])
    return T
    
def _calc_Tbase_for_2d_beam(nodes):
    
    x1,y1 = nodes[0].x,nodes[0].y
    x2,y2 = nodes[1].x,nodes[1].y
    le = np.sqrt((x2-x1)**2+(y2-y1)**2)

    lx = (x2-x1)/le
    mx = (y2-y1)/le
    T = np.array([[lx,mx,0.],
                  [-mx,lx,0.],
                  [0.,0.,1.]])
                  
    return T

if __name__ == "__main__":
    #materials
    E = 210e6
    A = 0.005
    I = 10e-5

    #nodes and elements
    n0 = Node(0,0)
    n1 = Node(0,3)
    n2 = Node(4,3)
    n3 = Node(4,0)
    n4 = Node(4,5)
    n5 = Node(8,5)
    n6 = Node(8,0)
    e0 = Beam2D11((n0,n1),E,A,I)
    e1 = BeamLink2D11((n1,n2),E,A,I)
    e2 = Beam2D11((n2,n3),E,A,I)
    e3 = Beam2D11((n2,n4),E,A,I)
    e4 = Beam2D11((n4,n5),E,A,I)
    e5 = Beam2D11((n5,n6),E,A,I)
    
    #system
    s = System()
    s.add_nodes([n0,n1,n2,n3,n4,n5,n6])
    s.add_elements([e0,e1,e2,e3,e4,e5])
    s.add_node_force(1,Fx = -10)
    s.add_node_force(5,Fx = -10)
    s.add_fixed_sup(0,3,6)
    s.solve()

    print n2.disp
    print e1.force
    
    
