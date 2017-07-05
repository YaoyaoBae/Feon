# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from __future__ import division
import numpy as np
from ..base import ElementBase
from ..tools import gl_quad2d

#~ *********************************************************************
#~ ************************* CLASS STRUCTELEMENT ***********************
#~ *********************************************************************


class StructElement(ElementBase):

    def __init__(self,nodes):
        ElementBase.__init__(self,nodes)
        self.init_nodes(nodes)
        self.init_unknowns()
        self._force =  dict.fromkeys(self.eIk,0.)
        self.dens = 2.
        self.t = 1.
        
    def init_nodes(self,nodes):
        v = np.array(nodes[0].coord)-np.array(nodes[1].coord)
        le = np.linalg.norm(v)
        self._volume = le
        

    def init_keys(self):
        if self.dim == 2:
            self.set_eIk(("N","Ty","Mz"))
        if self.dim == 3:
            self.set_eIk(("N","Ty","Tz","Mx","My","Mz"))
        

    def init_unknowns(self):
        pass
    
    
    @property
    def T(self):
        return self._T

    @property
    def me(self):
        return self._me
    
    @property
    def ke(self):
        return self._ke

    
    @property
    def force(self):
        return self._force

    @property
    def sx(self):
        sx = self._force[self.eIk[0]]/self.A
        return sx
    
    def calc_ke(self):
        pass

    def calc_me(self):
        pass
        
    def calc_T(self):
        pass
        
    def calc_Ke(self):
        self.calc_T()
        self.calc_ke()
        self._Ke = self.t*np.dot(np.dot(self.T.T,self.ke),self.T)

    def calc_Me(self):
        self.calc_T()
        self.calc_me()
        self._Me = np.dot(np.dot(self.T.T,self.me),self.T)
        
    def evaluate(self):
        u = np.array([[nd.disp[key] for nd in self.nodes for key in nd.nAk[:self.ndof]]])
        self._undealed_force = np.dot(self.T,np.dot(self.Ke,u.T))
        self.distribute_force()

    def distribute_force(self):
        n = len(self.eIk)
        for i,val in enumerate(self.eIk):
            self._force[val] += self._undealed_force[i::n]
            
    def load_equivalent(self,ltype,val):
        raise NotImplementedError
    

class SoildElement(ElementBase):

    def __init__(self,nodes):
        ElementBase.__init__(self,nodes)
        self.init_unknowns()
        self._stress =  dict.fromkeys(self.eIk,0.)
        self.init_nodes(nodes)
        self.dens = 2
        self.t = 1
        
    def init_nodes(self,nodes):
        pass

    def init_unknowns(self):
        pass
    
    @property
    def B(self):
        return self._B

    @property
    def D(self):
        return self._D
    
    @property
    def ke(self):
        return self._ke
    
    @property
    def me(self):
        return self._me

    
    @property
    def stress(self):
        return self._stress

    def calc_B(self):
        pass
        
    def calc_Ke(self):
        self.calc_B()
        self.calc_D()
        self._Ke = self.t*self.volume*np.dot(np.dot(self.B.T,self.D),self.B)
        
    
    def evaluate(self):
        u = np.array([[nd.disp[key] for nd in self.nodes for key in nd.nAk[:self.ndof]]])
        self._undealed_stress = np.dot(np.dot(self.D,self.B),u.T)
        self.distribute_stress()

    def distribute_stress(self):
        n = len(self.eIk)
        for i,val in enumerate(self.eIk):
            self._stress[val] += self._undealed_stress[i::n]
            
    def load_equivalent(self,ltype,val):
        raise NotImplementedError
    

#~ *********************************************************************
#~ ************************* CLASS PLANEELEMENT ************************
#~ *********************************************************************
    
class Tri2D11S(SoildElement):
    
    def __init__(self,nodes,E,nu,t = 1,dens= 2):
        SoildElement.__init__(self,nodes)
        self.E = E
        self.nu = nu
        self.t = t
        self.dens = dens

    def init_nodes(self,nodes):
        v1 = np.array(nodes[0].coord[:2])-np.array(nodes[1].coord[:2])
        v2 = np.array(nodes[0].coord[:2])-np.array(nodes[2].coord[:2])
        area = np.cross(v1,v2)/2.
        self._volume = area
        
    def init_keys(self):
        self.set_eIk(("sx","sy","sxy"))

    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Ux","Uy")
        self._ndof = 2

    def calc_D(self):
        self._D = _calc_D_for_tri2d11(E = self.E,nu = self.nu)

    def calc_B(self):
        self._B = _calc_B_for_tri2d11(self.nodes,self.volume)
    
    
    
def _calc_D_for_tri2d11(E = 1.,nu = 0.2):      
    a = E/(1-nu**2)
    D = a*np.array([[1.,nu,0.],
                    [nu,1.,0.],
                    [0.,0.,(1-nu)/2.]])
    return D

def _calc_B_for_tri2d11(nodes,area):
    x1,y1 = nodes[0].x,nodes[0].y
    x2,y2 = nodes[1].x,nodes[1].y
    x3,y3 = nodes[2].x,nodes[2].y
    belta1 = y2 - y3
    belta2 = y3 - y1
    belta3 = y1 - y2
    gama1 = x3 - x2
    gama2 = x1 - x3
    gama3 = x2 - x1

    return 1./(2.*area)*np.array([[belta1,0,belta2,0,belta3,0],
                                 [0.,gama1,0,gama2,0,gama3],
                                 [gama1,belta1,gama2,belta2,gama3,belta3]])


class Quad2D11S(SoildElement):

    def __init__(self,nodes,E,nu,t=1,dens = 2):
        SoildElement.__init__(self,nodes)
        self.E = E
        self.nu = nu
        self.t = t

    def init_keys(self):
        self.set_eIk(("sx","sy","sxy"))

    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Ux","Uy")
        self._ndof = 2
        
    def calc_D(self):
        self._D = _calc_D_for_quad2d11(self.E,self.nu)

    def calc_B(self):
        self._B,self._J = _calc_B_and_J_for_quad2d11(self.nodes,(0,0))

    def calc_Ke(self):
        self.calc_B()
        self._Ke = gl_quad2d(self.func,3)
        
    def func(self,x):
        self.calc_D()
        B,J = _calc_B_and_J_for_quad2d11(self.nodes,x)
        return self.t*np.dot(np.dot(B.T,self.D),B)*J

def _calc_D_for_quad2d11(E = 1.,nu = 0.2):      
    a = E/(1-nu**2)
    D = a*np.array([[1.,nu,0.],
                    [nu,1.,0.],
                    [0.,0.,(1-nu)/2.]])
    return D

def _calc_B_and_J_for_quad2d11(nodes,x):
    s = x[0]*1.0
    t = x[1]*1.0
    x1,y1 = nodes[0].x,nodes[0].y
    x2,y2 = nodes[1].x,nodes[1].y
    x3,y3 = nodes[2].x,nodes[2].y
    x4,y4 = nodes[3].x,nodes[3].y

    a = 1/4*(y1*(s-1)+y2*(-1-s)+y3*(1+s)+y4*(1-s))
    b = 1/4*(y1*(t-1)+y2*(1-t)+y3*(1+t)+y4*(-1-t))
    c = 1/4*(x1*(t-1)+x2*(1-t)+x3*(1+t)+x4*(-1-t))
    d = 1/4*(x1*(s-1)+x2*(-1-s)+x3*(1+s)+x4*(1-s))
    
    B100 = -1/4*a*(1-t)+1/4*b*(1-s)
    B111 = -1/4*c*(1-s)+1/4*d*(1-t)
    B120 = B111
    B121 = B100

    B200 = 1/4*a*(1-t)+1/4*b*(1+s)
    B211 = -1/4*c*(1+s)-1/4*d*(1-t)
    B220 = B211
    B221 = B200

    B300 = 1/4*a*(1+t)-1/4*b*(1+s)
    B311 = 1/4*c*(1+s)-1/4*d*(1+t)
    B320 = B311
    B321 = B300

    B400 = -1/4*a*(1+t)-1/4*b*(1-s)
    B411 = 1/4*c*(1-s)+1/4*d*(1+t)
    B420 = B411
    B421 = B400

    B = np.array([[B100,   0,B200,   0,B300,   0,B400,  0],
                  [0,   B111,0,   B211,0,   B311,0,  B411],
                  [B120,B121,B220,B221,B320,B321,B420,B421]])

    X = np.array([x1,x2,x3,x4])
    Y = np.array([y1,y2,y3,y4]).reshape(4,1)
    _J = np.array([[0,1-t,t-s,s-1],
                  [t-1,0,s+1,-s-t],
                  [s-t,-s-1,0,t+1],
                  [1-s,s+t,-t-1,0]])
    J = np.dot(np.dot(X,_J),Y)/8.

    return B/J,J

#~ *********************************************************************
#~ ************************* CLASS TETRAELEMENT ***********************
#~ *********************************************************************
class Tetra3D11(SoildElement):

    def __init__(self,nodes,E,nu,dens = 2):
        SoildElement.__init__(self,nodes)
        self.E = E
        self.nu = nu
        self.dens = dens

    def init_nodes(self,nodes):
        V = np.ones((4,4))
        for i,nd in enumerate(nodes):
            V[i,1:] = nd.coord
        self._volume = abs(np.linalg.det(V)/6.)
        
    def init_keys(self):
        self.set_eIk(("sx","sy","sz","sxy","syz","szx"))

    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Ux","Uy","Uz")
        self._ndof = 3

    def calc_B(self):
        self._B = _calc_B_for_tetra3d11(self.nodes,self.volume)
        
    def calc_D(self):
        self._D = _calc_D_for_tetra3d11(self.E,self.nu)
    

def _calc_B_for_tetra3d11(nodes,volume):
    A = np.ones((4,4))
    belta = np.zeros(4)
    gama = np.zeros(4)
    delta = np.zeros(4)
    for i,nd in enumerate(nodes):
        A[i,1:] = nd.coord

    for i in range(4):
        belta[i] = (-1)**(i+1)*np.linalg.det(np.delete(np.delete(A,i,0),1,1))
        gama[i] = (-1)**(i+2)*np.linalg.det(np.delete(np.delete(A,i,0),2,1))
        delta[i] = (-1)**(i+1)*np.linalg.det(np.delete(np.delete(A,i,0),3,1))

    B =  1./(6.*volume)*np.array([[belta[0],0.,0.,belta[1],0.,0.,belta[2],0.,0.,belta[3],0.,0.],
                                  [0.,gama[0],0.,0.,gama[1],0.,0.,gama[2],0.,0.,gama[3],0.],
                                  [0.,0.,delta[0],0.,0.,delta[1],0.,0.,delta[2],0.,0.,delta[3]],
                                  [gama[0],belta[0],0.,gama[1],belta[1],0.,gama[2],belta[2],0,gama[3],belta[3],0.],
                                  [0.,delta[0],gama[0],0.,delta[1],gama[1],0.,delta[2],gama[2],0.,delta[3],gama[3]],
                                  [delta[0],0.,belta[0],delta[1],0.,belta[1],delta[2],0.,belta[2],delta[3],0,belta[3]]])
    return B
def _calc_D_for_tetra3d11(E = 1.,nu = 1.):
    a = E/(1+nu)/(1-2*nu)
    b = 1.-nu
    c = (1.-2*nu)/2.
    T = a*np.array([[b,nu,nu,0.,0.,0.],
                    [nu,b,nu,0.,0.,0.],
                    [nu,nu,b,0.,0.,0.],
                    [0.,0.,0.,c,0.,0.],
                    [0.,0.,0.,0.,c,0.],
                    [0.,0.,0.,0.,0.,c]])
    return T


#~ *********************************************************************
#~ *****************************CLASS BEAM *****************************
#~ *********************************************************************

class Beam1D11(StructElement):

    def __init__(self,nodes,E,A,I,dens = 2):
        StructElement.__init__(self,nodes)
        self.E = E
        self.A = A
        self.I = I
        self.dens = dens

    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Uy","Phz")

        self._ndof = 3
                     
    def calc_T(self):
        TBase = _calc_Tbase_for_2d_beam(self.nodes)
        self._T = np.zeros((6,6))
        self._T[:3,:3] = self._T[3:,3:] = TBase
        
    def calc_ke(self):
        self._ke = _calc_ke_for_2d_beam(E = self.E,A = self.A,I = self.I,L = self.volume)
        
    def load_equivalent(self,ltype,val):
        self.calc_T()
        A = _calc_element_load_for_2d_beam(self,ltype = ltype,val = val)
        n = len(self.eIk)
        for i,key in enumerate(self.eIk):
            self._force[key] += -A[i::n]   
        B = np.dot(self.T.T,A)
        count = 0
        for nd in self.nodes:
            for key in nd.nBk:
                nd._force[key] += B[count,0]
                count +=1
        return B



class Beam2D11(StructElement):

    def __init__(self,nodes,E,A,I,dens = 2):
        StructElement.__init__(self,nodes)
        self.E = E
        self.A = A
        self.I = I
        self.dens = dens

    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Ux","Uy","Phz")

        self._ndof = 3
                     
    def calc_T(self):
        TBase = _calc_Tbase_for_2d_beam(self.nodes)
        self._T = np.zeros((6,6))
        self._T[:3,:3] = self._T[3:,3:] = TBase
        
    def calc_ke(self):
        self._ke = _calc_ke_for_2d_beam(E = self.E,A = self.A,I = self.I,L = self.volume)
        
    def load_equivalent(self,ltype,val):
        self.calc_T()
        A = _calc_element_load_for_2d_beam(self,ltype = ltype,val = val)
        n = len(self.eIk)
        for i,key in enumerate(self.eIk):
            self._force[key] += -A[i::n]   
        B = np.dot(self.T.T,A)
        count = 0
        for nd in self.nodes:
            for key in nd.nBk:
                nd._force[key] += B[count,0]
                count +=1
        return B
    
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

def _calc_ke_for_2d_beam(E = 1.0,A = 1.0,I = 1.0,L = 1.0):
    a00 =  E*A/L
    a03 = -a00
    a11 = 12*E*I/L**3
    a12 = 6*E*I/L**2
    a14 = -a11
    a22 = 4*E*I/L
    a24 = -a12
    a25 = 2*E*I/L
    a45 = -a12
    
    T = np.array([[a00, 0.,  0.,  a03,  0., 0.],
                  [ 0., a11, a12,  0., a14,a12],
                  [ 0., a12, a22,  0., a24,a25],
                  [a03,  0.,  0., a00,  0., 0.],
                  [ 0., a14, a24,  0.,a11, a45],
                  [ 0., a12, a25,  0.,a45, a22]])
    return T

def _calc_element_load_for_2d_beam(el,ltype = "q",val = 1.):
    le = el.volume
    A = np.zeros((6,1))
    if ltype in ["uniform","Uniform","UNIFORM","q","Q"]:
        A[1][0] = 1/2.*val*le
        A[4][0] = 1/2.*val*le
        A[2][0] = 1/12.*val*le**2
        A[5][0] = -1/12.*val*le**2
                    
    elif ltype in ["triangle","Triangle","TRIAGNLE","Tri","tri"]:
        A[1][0] = 3/20.*val*le
        A[4][0] = 7/20.*val*le
        A[2][0] = 1/30.*val*le**2
        A[5][0] = -1/20.*val*le**2
    else:
        raise AttributeError,"Unkown load type(%r)"%(ltype,)

    return A



##
class Beam3D11(StructElement):

    def __init__(self,nodes,E,G,A,I,dens = 2):
        StructElement.__init__(self,nodes)
        self.E = E
        self.G = G
        self.A = A
        self.Ix = I[0]
        self.Iy = I[1]
        self.Iz = I[2]
        self.dens = dens

    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Ux","Uy","Uz","Phx","Phy","Phz")

        self._ndof = 6
        
    def calc_T(self):
        TBase = _calc_Tbase_for_3d_beam(self.nodes)
        self._T = np.zeros((12,12))
        n = 3
        m = 4
        for i in range(m):
            self._T[n*i:n*(i+1),n*i:n*(i+1)] = TBase
            
    def calc_ke(self):
        self._ke = _calc_ke_for_3d_beam(E = self.E,
                                        G = self.G,
                                        A = self.A,
                                        I = [self.Ix,self.Iy,self.Iz],
                                        L = self.volume)
    def calc_me(self):
        self._me = _calc_me_for_3d_beam(E=self.E,
                                        G = self.G,
                                        A = self.A,
                                        I = self.Ix,
                                        L = self.volume,
                                        dens = self.dens)
        

    def load_equivalent(self,ltype,val):
        self.calc_T()
        A = _calc_element_load_for_3d_beam(self,ltype = ltype,val = val)
        n = len(self.eIk)
        for i,key in enumerate(self.eIk):
            self._force[key] += -A[i::n]
        B = np.dot(self.T.T,A)
        count = 0
        for nd in self.nodes:
            for key in nd.nBk:
                nd._force[key] += B[count,0]
                count +=1
        return B
    

##
def _calc_Tbase_for_3d_beam(nodes):
    x1,y1,z1 = nodes[0].x,nodes[0].y,nodes[0].z
    x2,y2,z2 = nodes[1].x,nodes[1].y,nodes[1].z
    if x1 == x2 and y1 == y2:
        if z2 > z1:
            return np.array([[0.,0.,1.],
                             [0.,1.,0.],
                             [-1.,0.,0.]])
        else:
            return np.array([[0.,0.,-1.],
                             [0.,1.,0.],
                             [1.,0.,0.]])
    else:
        le = np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        lx = (x2-x1)/le
        mx = (y2-y1)/le
        nx = (z2-z1)/le

        d = np.sqrt(lx**2+mx**2)

        ly = -mx/d
        my = lx/d
        ny = 0.

        lz = -lx*nx/d
        mz = -mx*nx/d
        nz = d
        return np.array([[lx,mx,nx],
                         [ly,my,ny],
                         [lz,mz,nz]])

def _calc_ke_for_3d_beam(E = 1.0,G = 1.0,A = 1.0,I = [1.,1.,1.],L = 1.):
    Ix = I[0]
    Iy = I[1]
    Iz = I[2]
    
    a00 = E*A/L
    a06 = -a00

    a11 = 12.*E*Iz/L**3
    a15 = 6.*E*Iz/L**2
    a17 = -a11

    a22 = 12.*E*Iy/L**3
    a24 = -6.*E*Iy/L**2
    a28 = -a22

    a33 = G*Ix/L
    a39 = -a33

    a44 = 4.*E*Iy/L
    a48 = 6.*E*Iy/L**2
    a410 = 2.*E*Iy/L

    a55 = 4.*E*Iz/L
    a57 = -a15
    a511 = 2.*E*Iz/L

    a711 = -a15
    a810 = a48

    
    K = np.array([[a00, 0,   0,  0,    0,   0, a06,    0,   0,   0,     0,     0],
                  [0, a11,   0,  0,    0, a15,   0,  a17,   0,   0,     0,   a15],
                  [0,   0, a22,  0,  a24,   0,   0,    0, a28,   0,   a24,     0],
                  [0,   0,   0,a33,   0,    0,   0,    0,   0, a39,     0,     0],
                  [0,   0, a24,  0,  a44,   0,   0,    0, a48,   0,  a410,     0],
                  [0, a15,   0,  0,    0, a55,   0,  a57,   0,   0,     0,  a511],
                  [a06, 0,   0,  0,    0,   0, a00,    0,   0,   0,     0,     0],
                  [0, a17,   0,  0,    0, a57,   0,  a11,   0,   0,     0,  a711],
                  [0,   0, a28,  0,  a48,   0,   0,    0, a22,   0,  a810,     0],
                  [0,   0,   0,a39,    0,   0,   0,    0,   0, a33,     0,     0],
                  [0,   0, a24,  0, a410,   0,   0,    0,a810,   0,   a44,     0],
                  [0, a15,   0,  0,    0,a511,   0, a711,   0,   0,     0,   a55]])
   
    return K

def _calc_me_for_3d_beam(E = 1.0,G = 1.0,A = 1.0,I = 1.,L = 1.,dens = 20.):
    a = dens*A*L/210.
    r = I/A
    

    
    M = a*np.array([[70.,     0.,     0.,    0.,        0.,        0.,    35.,     0.,      0.,      0.,         0.,         0.],
                    [ 0.,    78.,     0.,    0.,        0.,      11*L,     0.,    27.,      0.,      0.,         0.,     -6.5*L],
                    [ 0.,     0.,    78.,    0.,    -11.*L,        0.,     0.,     0.,     27.,      0.,      6.5*L,         0.],
                    [ 0.,     0.,     0., 70.*r,        0.,        0.,     0.,     0.,      0.,  -35.*r,         0.,         0.],
                    [ 0.,     0., -11.*L,    0.,    2*L**2,        0.,     0.,     0.,  -6.5*L,      0.,  -1.5*L**2,         0.],
                    [ 0.,   11*L,     0.,    0.,        0.,    2*L**2,     0.,  6.5*L,      0.,      0.,         0.,  -1.5*L**2],
                    [35.,     0.,     0.,    0.,        0.,        0.,    70.,     0.,      0.,      0.,         0.,         0.],
                    [ 0.,    27.,     0.,    0.,        0.,     6.5*L,     0.,    78.,      0.,      0.,         0.,     -11.*L],
                    [ 0.,     0.,    27.,    0.,    -6.5*L,        0.,     0.,     0.,     78.,      0.,       11*L,         0.],
                    [ 0.,     0.,     0.,-35.*r,        0.,        0.,     0.,     0.,      0.,   70.*r,         0.,         0.],
                    [ 0.,     0.,  6.5*L,    0., -1.5*L**2,        0.,     0.,     0.,   11.*L,      0.,    2.*L**2,         0.],
                    [ 0., -6.5*L,     0.,    0.,        0., -1.5*L**2,     0., -11.*L,      0.,      0.,         0.,    2.*L**2]])
   
    return M

def _calc_element_load_for_3d_beam(el,ltype = "q",val = (1,1)):
    
    assert type(val) is tuple or list,"3D frame element load must have two direction"
    if len(val) == 1:
        val = [val[0],0.]
    le = el.volume
    A = np.zeros((12,1))
    if ltype in ["uniform","Uniform","UNIFORM","q","Q"]:
        A[1][0] = 1/2.*val[0]*le
        A[7][0] = 1/2.*val[0]*le
        A[5][0] = 1/12.*val[0]*le**2
        A[11][0] = -1/12.*val[0]*le**2
        
        A[2][0] = 1/2.*val[1]*le
        A[8][0] = 1/2.*val[1]*le
        A[4][0] = 1/12.*val[1]*le**2
        A[10][0] = -1/12.*val[1]*le**2            
                    
    elif ltype in ["triangle","Triangle","TRIAGNLE1","Tri","tri"]:
        A[1][0] = 3/20.*val[0]*le
        A[7][0] = 7/20.*val[0]*le
        A[5][0] = 1/30.*val[0]*le**2
        A[11][0] = -1/20.*val[0]*le**2
        
        A[2][0] = 3/20.*val[1]*le
        A[8][0] = 7/20.*val[1]*le
        A[4][0] = 1/30.*val[1]*le**2
        A[10][0] = -1/20.*val[1]*le**2              


    else:
        raise AttributeError,"Unkown load type(%r)"%(ltype,)

    return A

#~ *********************************************************************
#~ *****************************  CLASS LINK ***************************
#~ *********************************************************************
class Link1D11(StructElement):
    def __init__(self,nodes,E,A,dens = 2):
        StructElement.__init__(self,nodes)
        self.E = E
        self.A = A
        self.dens = dens

        
    def init_keys(self):       
        self.set_eIk(["N"])

    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Ux")

        self._ndof = 1
        
    def calc_T(self):
        self._T = np.array([[1,0],
                            [0,1]])
        
    def calc_ke(self):
        self._ke = _calc_ke_for_link(E = self.E,A = self.A,L = self.volume)


        
class Link2D11(StructElement):

    def __init__(self,nodes,E,A,dens = 2):
        StructElement.__init__(self,nodes)
        self.E = E
        self.A = A
        self.dens = dens
        
    def init_keys(self):
        self.set_eIk(["N"])
        
    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Ux","Uy")

        self._ndof = 2

    def calc_T(self):
        self._T = _calc_Tbase_for_2d_link(self.nodes)
            
    def calc_ke(self):
        self._ke = _calc_ke_for_link(E = self.E,A = self.A,L = self.volume)


class Link3D11(StructElement):

    def __init__(self,nodes,E,A,dens = 2):
        StructElement.__init__(self,nodes)
        self.E = E
        self.A = A
        self.dens = dens
    
    def init_keys(self):
        self.set_eIk(["N"])

    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Ux","Uy","Uz")

        self._ndof = 3

    def calc_T(self):
        self._T = _calc_Tbase_for_3d_link(self.nodes)
            
    def calc_ke(self):
        self._ke = _calc_ke_for_link(E = self.E,
                                     A = self.A,
                                     L = self.volume)

def _calc_Tbase_for_2d_link(nodes):
    x1,y1 = nodes[0].x,nodes[0].y
    x2,y2 = nodes[1].x,nodes[1].y
    le = np.sqrt((x2-x1)**2+(y2-y1)**2)

    lx = (x2-x1)/le
    mx = (y2-y1)/le
    T = np.array([[lx,mx, 0.,0.],
                  [0., 0.,lx,mx]])
    return T

def _calc_Tbase_for_3d_link(nodes):
    x1,y1,z1 = nodes[0].x,nodes[0].y,nodes[0].z
    x2,y2,z2 = nodes[1].x,nodes[1].y,nodes[1].z

    le = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

    lx = (x2-x1)/le
    mx = (y2-y1)/le
    nx = (z2-z1)/le

    T = np.array([[lx,mx,nx,0.,0.,0.,],
                  [0.,0.,0.,lx,mx,nx]])

    return T

def _calc_ke_for_link(E = 1.,A = 1.,L = 1.):
    return np.array([[E*A/L,-E*A/L],
                     [-E*A/L,E*A/L]])

#~ *********************************************************************
#~ *****************************CLASS SPRING ***************************
#~ *********************************************************************
class Spring1D11(StructElement):
    def __init__(self,nodes,ke,dens = 2):
        StructElement.__init__(self,nodes)
        self.k = ke
        self.dens = dens

    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Ux")

        self._ndof = 1
            
    def init_keys(self):
        self.set_eIk(["N"])
        
    def calc_T(self):
        self._T = np.array([[1,0],
                            [0,1]])
    def sx(self):
        pass
            
    def calc_ke(self):
        self._ke = _calc_ke_for_spring(ke = self.k)

class Spring2D11(StructElement):

    def __init__(self,nodes,ke,dens = 2):
        StructElement.__init__(self,nodes)
        self.k = ke
        self.dens = dens

    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Ux","Uy")

        self._ndof = 2
            

    @property
    def sx(self):
        pass
        
    def init_keys(self):
        self.set_eIk(["N"])
    
    def calc_T(self):
        self._T = _calc_Tbase_for_2d_link(self.nodes)
            
    def calc_ke(self):
        self._ke = _calc_ke_for_spring(ke = self.k)
        
        
class Spring3D11(StructElement):

    def __init__(self,nodes,ke,dens = 20):
        StructElement.__init__(self,nodes)
        self.k = ke
        self.dens = dens

    def init_unknowns(self):
        for nd in self.nodes:
            nd.init_unknowns("Ux","Uy","Uz")

        self._ndof = 3
        
    def init_keys(self):
        self.set_eIk(["N"])

    @property
    def sx(self):
        pass
    
    def calc_T(self):
        self._T = _calc_Tbase_for_3d_link(self.nodes)
            
    def calc_ke(self):
        self._ke = _calc_ke_for_spring(ke = self.k)

def _calc_ke_for_spring(ke = 1.0):
    return np.array([[ ke,-ke],
                     [ -ke,ke]])

























#~ *********************************************************************
#~ *****************************  FUNCTIONS ****************************
#~ *********************************************************************



def _calcKeLocFor2DLinkBeam(E = 1.0,A = 1.0,I = 1.0,L = 1.0):
    a00 = E*A/L
    a03 = -a00
    a11 = 3*E*I/L**3
    a12 = 3*E*I/L**2
    a14 = -a11
    a22 = 3*E*I/L
    T = np.array([[a00,  0.,  0,  a03,  0.,  0.],
                  [ 0., a11,  0.,  0., a14, a12],
                  [ 0.,  0.,  0.,  0.,  0.,  0.],
                  [a03,  0.,  0., a00,  0.,  0.],
                  [ 0., a14,  0.,  0., a11, -a12],
                  [ 0., a12,  0.,  0., -a12, a22]])
    return T

def _calcKeLocFor2DBeamLink(E = 1.0,A = 1.0,I = 1.0,L = 1.0):
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



def _calcElementLoadFor2DLinkBeam(el,ltype = "q",val = 1.):
    le = el.volume
    A = np.array([[0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.]])
    if ltype in ["uniform","Uniform","UNIFORM","q","Q"]:
        A[1][0] = 1/2.*val*le - 3/4.*val
        A[4][0] = 1/2.*val*le + 3/4.*val
        A[5][0] = -1/24.*val*le**2
                    
    elif ltype in ["triangle","Triangle","TRIAGNLE","Tri","tri"]:
        A[1][0] = 3/20.*val*le - 9/40.*val
        A[4][0] = 7/20.*val*le + 21/40.*val
        A[5][0] = -1/40.*val*le**2
            

    else:
        raise AttributeError,"Unkown load type(%r)"%(ltype,)

    return A


def _calcElementLoadFor2DBeamLink(el,ltype = "q",val = 1.):
    le = el.volume
    A = np.array([[0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.]])
    if ltype in ["uniform","Uniform","UNIFORM","q","Q"]:
        A[1][0] = 1/2.*val*le - 3/4.*val
        A[4][0] = 1/2.*val*le + 3/4.*val
        A[2][0] = 1/24.*val*le**2
                    
    elif ltype in ["triangle","Triangle","TRIAGNLE","Tri","tri"]:
        A[1][0] = 3/20.*val*le - 9/40.*val
        A[4][0] = 7/20.*val*le + 21/40.*val
        A[2][0] = 1/40.*val*le**2
            

    else:
        raise AttributeError,"Unkown load type(%r)"%(ltype,)

    return A




def _calcMeLocFor2DBeam(dens = 20.,A = 1.0,I = 1.0,L = 1.0):
    a = dens*L/210.
    
    T = a*np.array([[70.,       0.,          0.,  35.,       0.,       0.],
                    [ 0.,      78.,       11.*L,   0.,      27.,   -6.5*L],
                    [ 0.,    11.*L,      2*L**2,   0.,    6.5*L,-1.5*L**2],
                    [35.,       0.,          0.,  70.,       0.,       0.],
                    [ 0.,      27.,       6.5*L,   0.,      78.,   -11.*L],
                    [ 0.,   -6.5*L,   -1.5*L**2,   0.,   -11.*L,   2*L**2]])
    return T




def _calcKeLocFor3DLinkBeam(E = 1.0,G = 1.0,A = 1.0,I = [1.,1.,1.],L = 1.):
    Ix = I[0]
    Iy = I[1]
    Iz = I[2]
    
    a00 = E*A/L
    a06 = -a00

    a11 = 3.*E*Iz/L**3
    a15 = 3.*E*Iz/L**2
    a17 = -a11

    a22 = 3.*E*Iy/L**3
    a24 = -3.*E*Iy/L**2
    a28 = -a22

    a33 = G*Ix/L

    a44 = 3.*E*Iy/L

    a55 = 3.*E*Iz/L

    a711 = -a15

    a810 = -a24

    
    T = np.array([[a00, 0,   0,  0,    0,   0, a06,    0,   0,   0,     0,     0],
                  [0, a11,   0,  0,    0,   0,   0,  a17,   0,   0,     0,   a15],
                  [0,   0, a22,  0,    0,   0,   0,    0, a28,   0,   a24,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [a06, 0,   0,  0,    0,   0, a00,    0,   0,   0,     0,     0],
                  [0, a17,   0,  0,    0,   0,   0,  a11,   0,   0,     0,  a711],
                  [0,   0, a28,  0,    0,   0,   0,    0, a22,   0,  a810,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0, a33,     0,     0],
                  [0,   0, a24,  0,    0,   0,   0,    0,a810,   0,   a44,     0],
                  [0, a15,   0,  0,    0,   0,   0, a711,   0,   0,     0,   a55]])
   
    return T

def _calcKeLocFor3DBeamLink(E = 1.0,G = 1.0,A = 1.0,I = [1.,1.,1.],L = 1.):
    Ix = I[0]
    Iy = I[1]
    Iz = I[2]
    
    a00 = E*A/L
    a06 = -a00

    a11 = 3.*E*Iz/L**3
    a15 = 3.*E*Iz/L**2
    a17 = -a11

    a22 = 3.*E*Iy/L**3
    a24 = 3.*E*Iy/L**2
    a28 = -a22

    a33 = G*Ix/L

    a44 = 3.*E*Iy/L
    a48 = 3.*E*Iy/L**2

    a55 = 3.*E*Iz/L
    a57 = -a15

    
    T = np.array([[a00, 0,   0,  0,    0,   0, a06,    0,   0,   0,     0,     0],
                  [0, a11,   0,  0,    0, a15,   0,  a17,   0,   0,     0,     0],
                  [0,   0, a22,  0,  a24,   0,   0,    0, a28,   0,     0,     0],
                  [0,   0,   0,a33,   0,    0,   0,    0,   0,   0,     0,     0],
                  [0,   0, a24,  0,  a44,   0,   0,    0, a48,   0,     0,     0],
                  [0, a15,   0,  0,    0, a55,   0,  a57,   0,   0,     0,     0],
                  [a06, 0,   0,  0,    0,   0, a00,    0,   0,   0,     0,     0],
                  [0, a17,   0,  0,    0, a57,   0,  a11,   0,   0,     0,     0],
                  [0,   0, a28,  0,  a48,   0,   0,    0, a22,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0]])
   
    return T

def _calcKeLocFor3DLink(E = 1.0,A = 1.0,L = 1.):
    a00 = E*A/L
    a06 = -a00



    
    T = np.array([[a00, 0,   0,  0,    0,   0, a06,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [a06, 0,   0,  0,    0,   0, a00,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0],
                  [0,   0,   0,  0,    0,   0,   0,    0,   0,   0,     0,     0]])
   
    return T





def _calcElementLoadFor3DLinkBeam(el,ltype = "q",val = (1.,1.)):
    assert type(val) is tuple or list,"3D frame element load must have two direction"
    if len(val) == 1:
        val = [val[0],0.]
    le = el.volume
    A = np.array([[0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.]])
    if ltype in ["uniform","Uniform","UNIFORM","q","Q"]:
        A[1][0] = 1/2.*val[0]*le - 3/4.*val[0]
        A[7][0] = 1/2.*val[0]*le + 3/4.*val[0]
        A[11][0] = -1/24.*val[0]*le**2
        
        A[2][0] = 1/2.*val[1]*le - 3/4.*val[1]
        A[8][0] = 1/2.*val[1]*le + 3/4.*val[1]
        A[10][0] = -1/24.*val[1]*le**2            
                    
    elif ltype in ["triangle","Triangle","TRIAGNLE","Tri","tri"]:
        A[1][0] = 3/20.*val[0]*le - 9/40.*val[0]
        A[7][0] = 7/20.*val[0]*le + 21/40.*val[0]
        A[11][0] = -1/40.*val[0]*le**2
        
        A[2][0] = 3/20.*val[1]*le - 9/40.*val[1]
        A[8][0] = 7/20.*val[1]*le + 21/40.*val[1]
        A[10][0] = -1/40.*val[1]*le**2              


    else:
        raise AttributeError,"Unkown load type(%r)"%(ltype,)

    return A

def _calcElementLoadFor3DBeamLink(el,ltype = "q",val = (1.,1.)):
    assert type(val) is tuple or list,"3D frame element load must have two direction"
    if len(val) == 1:
        val = [val[0],0.]
    le = el.volume
    A = np.array([[0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.],
                  [0.]])
    if ltype in ["uniform","Uniform","UNIFORM","q","Q"]:
        A[1][0] = 1/2.*val[0]*le - 3/4.*val[0]
        A[7][0] = 1/2.*val[0]*le + 3/4.*val[0]
        A[5][0] = 1/24.*val[0]*le**2
        
        A[2][0] = 1/2.*val[1]*le - 3/4.*val[1]
        A[8][0] = 1/2.*val[1]*le + 3/4.*val[1]
        A[4][0] = 1/24.*val[1]*le**2            
                    
    elif ltype in ["triangle","Triangle","TRIAGNLE","Tri","tri"]:
        A[1][0] = 3/20.*val[0]*le - 9/40.*val[0]
        A[7][0] = 7/20.*val[0]*le + 21/40.*val[0]
        A[5][0] = 1/60.*val[0]*le**2
        
        A[2][0] = 3/20.*val[1]*le - 9/40.*val[1]
        A[8][0] = 7/20.*val[1]*le + 21/40.*val[1]
        A[4][0] = 1/60.*val[1]*le**2              

    else:
        raise AttributeError,"Unkown load type(%r)"%(ltype,)

    return A



if __name__ == "__main__":
    E = 210e6
    nu = 0.3
    t = 0.025
    n1 = Node(0,0)
    n2 = Node(0.5,0)
    n3 = Node(0.5,0.25)
    n4 = Node(0,0.25)

    n5 = Node(0,0,0)
    n6 = Node(0.025,0,0)
    n7 = Node(0.025,0.5,0)
    n8 = Node(0.025,0,0.25)
    e1 = Tri2D11S((n1,n3,n4),E,nu,t)
    e2 = Tri2D11S((n1,n2,n3),E,nu,t)
    e3 = Beam3D11((n5,n6),1,1,1,[1,1,1])
    e4 = Tetra3D11((n5,n6,n7,n8),E,nu)
    
    
