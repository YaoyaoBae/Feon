# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from __future__ import division
from feon.sa import *
from feon.mesh import Mesh
from feon.tools import gl_quad3d
import numpy as np
class Brick3D11(SoildElement):

    def __init__(self,nodes,E,nu):
        SoildElement.__init__(self,nodes)
        self.E = E
        self.nu = nu
   
    def init_keys(self):
        self.set_eIk(("sx","sy","sz","sxy","syz","szx"))

    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Ux","Uy","Uz")

        self._ndof = 3
        
    def calc_D(self):
        self._D = _calc_D_for_brick3d11(self.E,self.nu)

    def calc_B(self):
        self._B,self._J= _calc_B_and_J_for_brick3d11(self.nodes,(0,0,0))

    def calc_Ke(self):
        self.calc_B()
        self._Ke = gl_quad3d(self.func,3)
        
    def func(self,x):
        self.calc_D()
        B,J = _calc_B_and_J_for_brick3d11(self.nodes,x)
        return np.dot(np.dot(B.T,self.D),B)*J


def _calc_D_for_brick3d11(E = 1.,nu = 0.3):      
    a = E/((1+nu)*(1-2*nu))
    D = a*np.array([[1-nu,nu,nu,0,0,0],
                    [nu,1-nu,nu,0,0,0],
                    [nu,nu,1-nu,0,0,0],
                    [0,0,0,(1-nu)/2.,0,0],
                    [0,0,0,0,(1-nu)/2.,0],
                    [0,0,0,0,0,(1-nu)/2.]])
    return D

def _calc_B_and_J_for_brick3d11(nodes,x):
    s = x[0]
    t = x[1]
    u = x[2]
    
    x = [nd.x for nd in nodes]
    y = [nd.y for nd in nodes]
    z = [nd.z for nd in nodes]

    N1s,N1t,N1u = (t-1)*(1+u)/8.,(s-1)*(1+u)/8.,(1-s)*(1-t)/8.
    N2s,N2t,N2u = (t-1)*(1-u)/8.,(s-1)*(1-u)/8.,(s-1)*(1-t)/8.
    N3s,N3t,N3u = (t+1)*(u-1)/8.,(1-s)*(1-u)/8.,(s-1)*(1+t)/8.
    N4s,N4t,N4u = (1+t)*(-u-1)/8.,(1-s)*(1+u)/8.,(1-s)*(1+t)/8.
    N5s,N5t,N5u = (1-t)*(1+u)/8.,(-s-1)*(1+u)/8.,(1+s)*(1-t)/8.
    N6s,N6t,N6u = (1-t)*(1-u)/8.,(-1-s)*(1-u)/8.,(1+s)*(t-1)/8.
    N7s,N7t,N7u = (1+t)*(1-u)/8.,(1+s)*(1-u)/8.,(-1-s)*(t+1)/8.
    N8s,N8t,N8u = (1+t)*(1+u)/8.,(1+s)*(1+u)/8.,(1+s)*(t+1)/8.

    
    Ns = [N1s,N2s,N3s,N4s,N5s,N6s,N7s,N8s]
    Nt = [N1t,N2t,N3t,N4t,N5t,N6t,N7t,N8t]
    Nu = [N1u,N2u,N3u,N4u,N5u,N6u,N7u,N8u]

    
    xs = sum(Ns[i]*x[i] for i in xrange(8))
    xt = sum(Nt[i]*x[i] for i in xrange(8))
    xu = sum(Nu[i]*x[i] for i in xrange(8))
    
    ys = sum(Ns[i]*y[i] for i in xrange(8))
    yt = sum(Nt[i]*y[i] for i in xrange(8))
    yu = sum(Nu[i]*y[i] for i in xrange(8))
    
    zs = sum(Ns[i]*z[i] for i in xrange(8))
    zt = sum(Nt[i]*z[i] for i in xrange(8))
    zu = sum(Nu[i]*z[i] for i in xrange(8))


    MJ = np.array([[xs,ys,zs],
                   [xt,yt,zt],
                   [xu,yu,zu]])
    
    J = np.linalg.det(MJ)
    J_v = np.linalg.inv(MJ)

    Nx = [J_v[0,0]*Ns[i] + J_v[0,1]*Nt[i] + J_v[0,2]*Nu[i] for i in xrange(8)]
    Ny = [J_v[1,0]*Ns[i] + J_v[1,1]*Nt[i] + J_v[1,2]*Nu[i] for i in xrange(8)]
    Nz = [J_v[2,0]*Ns[i] + J_v[2,1]*Nt[i] + J_v[2,2]*Nu[i] for i in xrange(8)]
    
    B =  np.array([[Nx[0],0,0,Nx[1],0,0,Nx[2],0,0,Nx[3],0,0,Nx[4],0,0,Nx[5],0,0,Nx[6],0,0,Nx[7],0,0],
                   [0,Ny[0],0,0,Ny[1],0,0,Ny[2],0,0,Ny[3],0,0,Ny[4],0,0,Ny[5],0,0,Ny[6],0,0,Ny[7],0],
                   [0,0,Nz[0],0,0,Nz[1],0,0,Nz[2],0,0,Nz[3],0,0,Nz[4],0,0,Nz[5],0,0,Nz[6],0,0,Nz[7]],
                   [Ny[0],Nx[0],0,Ny[1],Nx[1],0,Ny[2],Nx[2],0,Ny[3],Nx[3],0,Ny[4],Nx[4],0,Ny[5],Nx[5],0,Ny[6],Nx[6],0,Ny[7],Nx[7],0],
                   [0,Nz[0],Ny[0],0,Nz[1],Ny[1],0,Nz[2],Ny[2],0,Nz[3],Ny[3],0,Nz[4],Ny[4],0,Nz[5],Ny[5],0,Nz[6],Ny[6],0,Nz[7],Ny[7]],
                   [Nz[0],0,Nx[0],Nz[1],0,Nx[1],Nz[2],0,Nx[2],Nz[3],0,Nx[3],Nz[4],0,Nx[4],Nz[5],0,Nx[5],Nz[6],0,Nx[6],Nz[7],0,Nx[7]]])
    return B,J

if __name__ == "__main__":
    mesh = Mesh()
    mesh.build(mesh_type = "cube",x_lim = [0,0.5],y_lim = [0,0.025],z_lim=[0,0.25],size = [2,1,1])
    E = 210e6
    nu = 0.3
    nds = np.array([Node(p) for p in mesh.points])
    els = [Brick3D11(nds[c],E,nu) for c in mesh.elements]
    s = System()
    s.add_nodes(nds)
    s.add_elements(els)
    fixed_nds = [nd.ID for nd in nds if nd.x == 0.]
    add_force_nds = [nd.ID for nd in nds if nd.x == 0.5]
    s.add_fixed_sup(fixed_nds)
    for nid in add_force_nds:
       s.add_node_force(nid,Fx = 4.6875)
    s.solve()
    

