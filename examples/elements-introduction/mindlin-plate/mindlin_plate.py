# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from __future__ import division
from feon.sa import *
from feon.mesh import Mesh
from feon.tools import gl_quad2d
import numpy as np
class Plate11(StructElement):

    def __init__(self,nodes,E,G,nu,t):
        StructElement.__init__(self,nodes)
        self.E = E
        self.G = G
        self.nu = nu
        self.t = t
        

    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Uz","Phx","Phy")

        self._ndof = 6
        
    def calc_D(self):
        nu = self.nu
        a = self.E/(1-self.nu**2)
        self.Dc = np.array([[self.G,0],
                            [0,self.G]])
        self.Df = a*np.array([[1,nu,0],
                              [nu,1,0],
                              [0,0,(1-nu)/2.]])
    def calc_T(self):
        self._T = np.eye(24)

    def calc_B(self):
        self.calc_D()
        Bf,Bc,self._J= _calc_B_and_J_for_plate3d11(self.nodes,(0,0,0))
        DJ = np.linalg.det(self.J)
        B = self.t**3/12.*np.dot(np.dot(Bf.T,self.Df),Bf)*DJ + 5/6.*self.t*np.dot(np.dot(Bc.T,self.Dc),Bc)*DJ
        _B = np.zeros((24,24))
        _B[8:20,8:20] = B
        self._B = _B
    
    def calc_Ke(self):
        self.calc_T()
        self.calc_B()
        Ke = gl_quad2d(self.func,3)
        _Ke = np.zeros((24,24))
        _Ke[8:20,8:20] = Ke
        self._Ke = _Ke
        
        
    def func(self,x):
        self.calc_D()
        Bf,Bc,J = _calc_B_and_J_for_plate3d11(self.nodes,x)
        DJ = np.linalg.det(J)
        return self.t**3/12.*np.dot(np.dot(Bf.T,self.Df),Bf)*DJ + 5/6.*self.t*np.dot(np.dot(Bc.T,self.Dc),Bc)*DJ

    def func_jac(self,x):
        Bf,Bc,J = _calc_B_and_J_for_plate3d11(self.nodes,x)
        return np.linalg.det(J)

    @property
    def voa(self):
        return gl_quad2d(self.func_jac,2)
        

def _calc_B_and_J_for_plate3d11(nodes,x):
    s = x[0]
    t = x[1]

    x = [nd.x for nd in nodes]
    y = [nd.y for nd in nodes]

    N1,N2,N3,N4 = 1/4*(1-s)*(1-t),1/4*(1+s)*(1-t),1/4*(1+s)*(1+t),1/4*(1-s)*(1+t)

    N1s,N1t = 1/4*(t-1),1/4*(s-1)
    N2s,N2t = 1/4*(1-t),1/4*(-s-1)
    N3s,N3t = 1/4*(1+t),1/4*(1+s)
    N4s,N4t = 1/4*(-1-t),1/4*(1-s)

    Ns = [N1s,N2s,N3s,N4s]
    Nt = [N1t,N2t,N3t,N4t]
    
    xs = sum(Ns[i]*x[i] for i in xrange(4))
    xt = sum(Nt[i]*x[i] for i in xrange(4))
   
    ys = sum(Ns[i]*y[i] for i in xrange(4))
    yt = sum(Nt[i]*y[i] for i in xrange(4))
    
    J = np.array([[xs,ys],
                  [xt,yt]])
    
    J_v = np.linalg.inv(J)

    Nx = [J_v[0,0]*Ns[i] + J_v[0,1]*Nt[i] for i in xrange(4)]
    Ny = [J_v[1,0]*Ns[i] + J_v[1,1]*Nt[i] for i in xrange(4)]
    
    Bf = np.array([[0,Nx[0],0,0,Nx[1],0,0,Nx[2],0,0,Nx[3],0],
                   [0,0,Ny[0],0,0,Ny[1],0,0,Ny[2],0,0,Ny[3]],
                   [0,Ny[0],Nx[0],0,Ny[1],Nx[1],0,Ny[2],Nx[2],0,Ny[2],Nx[2]]])

    Bc = np.array([[Nx[0],N1,0,Nx[1],N2,0,Nx[2],N3,0,Nx[3],N4,0],
                   [Ny[0],0,N1,Ny[1],0,N2,Ny[2],0,N3,Ny[3],0,N4]])
    return Bf,Bc,J
if __name__ == "__main__":
    mesh = Mesh()
    mesh.build(mesh_type = "rect",x_lim = [0,1],y_lim = [0,1],size=[2,2])
    E = 210e6
    G = 84e6
    nu = 0.3
    t = 0.05
    nds = np.array([Node(p[0],p[1],0) for p in mesh.points])
    els = [Plate11(nds[c],E,G,nu,t) for c in mesh.elements]

    s= System()
    s.add_nodes(nds)
    s.add_elements(els)

    fixed_nds = [nd.ID for nd in nds if nd.x ==0 or nd.x ==1 or nd.y == 0 or nd.y ==1]
    s.add_fixed_sup(fixed_nds)
    s.add_node_force(4,Fz=-20)
    s.solve()
        

