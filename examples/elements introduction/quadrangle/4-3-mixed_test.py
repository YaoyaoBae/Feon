# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------


##from __future__ import division
from feon.sa import *
##import numpy as np
##from feon.tools import gl_quad2d

##class Quad2D11S(SoildElement):
##
##    def __init__(self,nodes,E,nu,t=1):
##        SoildElement.__init__(self,nodes)
##        self.E = E
##        self.nu = nu
##        self.t = t
##
##    def init_keys(self):
##        self.set_eIk(("sx","sy","sxy"))
##
##    def init_unknowns(self):
##        for nd in self.nodes:
##            nd.init_unknowns("Ux","Uy")
##
##        self._ndof = 2
##        
##    def calc_D(self):
##        self._D = _calc_D_for_quad2d11(self.E,self.nu)
##
##    def calc_B(self):
##        self._B,self.J = _calc_B_and_J_for_quad2d11(self.nodes,(0,0))
##
##    def calc_Ke(self):
##        self.calc_B()
##        self._Ke = gl_quad2d(self.func,2)
##        
##    def func(self,x):
##        self.calc_D()
##        B,J = _calc_B_and_J_for_quad2d11(self.nodes,x)
##        return self.t*np.dot(np.dot(B.T,self.D),B)*J
##
##def _calc_D_for_quad2d11(E = 1.,nu = 0.2):      
##    a = E/(1-nu**2)
##    D = a*np.array([[1.,nu,0.],
##                    [nu,1.,0.],
##                    [0.,0.,(1-nu)/2.]])
##    return D
##
##def _calc_B_and_J_for_quad2d11(nodes,x):
##    s = x[0]*1.0
##    t = x[1]*1.0
##    x1,y1 = nodes[0].x,nodes[0].y
##    x2,y2 = nodes[1].x,nodes[1].y
##    x3,y3 = nodes[2].x,nodes[2].y
##    x4,y4 = nodes[3].x,nodes[3].y
##
##    a = 1/4*(y1*(s-1)+y2*(-1-s)+y3*(1+s)+y4*(1-s))
##    b = 1/4*(y1*(t-1)+y2*(1-t)+y3*(1+t)+y4*(-1-t))
##    c = 1/4*(x1*(t-1)+x2*(1-t)+x3*(1+t)+x4*(-1-t))
##    d = 1/4*(x1*(s-1)+x2*(-1-s)+x3*(1+s)+x4*(1-s))
##    
##    B100 = -1/4*a*(1-t)+1/4*b*(1-s)
##    B111 = -1/4*c*(1-s)+1/4*d*(1-t)
##    B120 = B111
##    B121 = B100
##
##    B200 = 1/4*a*(1-t)+1/4*b*(1+s)
##    B211 = -1/4*c*(1+s)-1/4*d*(1-t)
##    B220 = B211
##    B221 = B200
##
##    B300 = 1/4*a*(1+t)-1/4*b*(1+s)
##    B311 = 1/4*c*(1+s)-1/4*d*(1+t)
##    B320 = B311
##    B321 = B300
##
##    B400 = -1/4*a*(1+t)-1/4*b*(1-s)
##    B411 = 1/4*c*(1-s)+1/4*d*(1+t)
##    B420 = B411
##    B421 = B400
##
##    B = np.array([[B100,   0,B200,   0,B300,   0,B400,  0],
##                  [0,   B111,0,   B211,0,   B311,0,  B411],
##                  [B120,B121,B220,B221,B320,B321,B420,B421]])
##
##    X = np.array([x1,x2,x3,x4])
##    Y = np.array([y1,y2,y3,y4]).reshape(4,1)
##    _J = np.array([[0,1-t,t-s,s-1],
##                  [t-1,0,s+1,-s-t],
##                  [s-t,-s-1,0,t+1],
##                  [1-s,s+t,-t-1,0]])
##    J = np.dot(np.dot(X,_J),Y)/8.
##
##    return B/J,J
    
    

if __name__ == "__main__":
    E = 210e6
    nu = 0.3
    t = 0.025
    
    n0 = Node(0,0)
    n1 = Node(0.25,0)
    n2 = Node(0.25,0.25)
    n3 = Node(0,0.25)
    n4 = Node(0.50,0)
    n5 = Node(0.50,0.25)
    e0 = Quad2D11S((n0,n1,n2,n3),E,nu,t)
    e1 = Quad2D11S((n1,n4,n5,n2),E,nu,t)

    s = System()
    s.add_nodes(n0,n1,n2,n3,n4,n5)
    s.add_elements(e0,e1)
    s.add_node_force(4,Fx = 9.75)
    s.add_node_force(5,Fx = 9.75)
    s.add_fixed_sup(0,3)
    s.solve()
    

