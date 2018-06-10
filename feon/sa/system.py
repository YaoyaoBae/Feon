# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -----------------------------------


import numpy as np
from .solver import *
from ..base import SystemBase
from .post_process import PostProcess
##from draw2d import Drawer
#~ *********************************************************************
#~ ****************************  CLASS SYSTEM **************************
#~ *********************************************************************

class System(SystemBase):

    def __init__(self):
        SystemBase.__init__(self)
        self._Force = {}
        self._Disp = {}        
        self._is_inited = False
        self._is_force_added = False
        self._is_disp_added = False
        self._is_system_solved = False
    

    def __repr__(self):
        return "%dD System: \nNodes: %d\nElements: %d"\
               %(self.dim,self.non,self.noe,)
    
    @property
    def Force(self):
        return self._Force


    @property
    def Disp(self):
        return self._Disp

    @property
    def DispValue(self):
        return self._DispValue

    @property
    def ForceValue(self):
        return self._ForceValue
    
    @property
    def KG(self):
        return self._KG
    @property
    def MG(self):
        return self._MG

    @property
    def KG_keeped(self):
        return self._KG_keeped

    @property
    def MG_keeped(self):
        return self._MG_keeped

    @property
    def Force_keeped(self):
        return self._Force_keeped

    @property
    def Disp_keeped(self):
        return self._Disp_keeped

    @property
    def deleted(self):
        return self._deleted

    @property
    def keeped(self):
        return self._keeped

    @property
    def nonzeros(self):
        return self._nonzeros
    


    def init(self):
        self._mndof = max(el.ndof for el in self.get_elements())
        self._nAk = self.nodes[0].nAk[:self.mndof]
        self._nBk = self.nodes[0].nBk[:self.mndof]
        self._dim = self.nodes[0].dim

    def calc_KG(self):
        self.init()
        n = self.non
        m = self.mndof
        shape = n*m
        self._KG = np.zeros((shape,shape))
        for el in self.get_elements():
            ID = [nd.ID for nd in el.nodes]
            el.calc_Ke()
            M = el.ndof
            for N1,I in enumerate(ID):
                for N2,J in enumerate(ID):
                    self._KG[m*I:m*I+M,m*J:m*J+M] += el.Ke[M*N1:M*(N1+1),M*N2:M*(N2+1)]

        self._is_inited = True
        
    def calc_MG(self):
        self.init()
        n = self.non
        m = self.mndof
        shape = n*m
        self._MG = np.zeros((shape,shape))
        for el in self.get_elements():
            ID = [nd.ID for nd in el.nodes]
            el.calc_Me()
            M = el.ndof
            for N1,I in enumerate(ID):
                for N2,J in enumerate(ID):
                    self._MG[m*I:m*I+M,m*J:m*J+M] += el.Me[M*N1:M*(N1+1),M*N2:M*(N2+1)]
        
    
    def add_element_load(self,eid,ltype,val): #only for beam type
        if not self._is_inited:
            self.calc_KG()
        assert eid <= self.noe,"Element does not exist"

        B = self.elements[eid].load_equivalent(ltype = ltype,val = val)

        self._is_force_added = True
        
    def add_node_force(self,nid,**forces):
        if not self._is_inited:
            self.calc_KG()
            
        assert nid+1 <= self.non,"Element does not exist"
        for key in forces.keys():
            assert key in self.nBk,"Check if the node forces applied are correct"

        self.nodes[nid].set_force(**forces)   
        self._is_force_added = True
        

    def add_node_disp(self,nid,**disp):
        if not self._is_inited:
            self.calc_KG()
        assert nid+1 <= self.non,"Element does not exist"
        for key in disp.keys():
            assert key in self.nAk,"Check if the node disp applied are correct"  
        self.nodes[nid].set_disp(**disp)
        val = disp.values()
        if len(val):
            self._is_disp_added = True
        

            
    def add_fixed_sup(self,*nids):
        if not self._is_inited:
            self.calc_KG()
        for nid in nids:
            if isinstance(nid,list) or isinstance(nid,tuple) or isinstance(nid,np.ndarray):
                for n in nid:
                    for key in self.nAk:
                        self.nodes[n]._disp[key] = 0.
            else:
                for key in self.nAk:
                    self.nodes[nid]._disp[key] = 0.

    def add_hinged_sup(self,*nids):
        if not self._is_inited:
            self.calc_KG()
        for nid in nids:
            if isinstance(nid,list) or isinstance(nid,tuple) or isinstance(nid,np.ndarray):
                for n in nid:
                    for key in self.nAk[:-1]:
                        self.nodes[n]._disp[key] = 0.
            else:
                for key in self.nAk[:-1]:
                    self.nodes[nid]._disp[key] = 0.
                    
    def add_rolled_sup(self,nid,direction = "x"):
        if not self._is_inited:
            self.calc_KG()
        if self.dim == 2:
            assert direction in ["x","y"],"Support dirction is x,y"
            if direction is "x":
                self.nodes[nid].set_disp(Ux = 0.)
            
            elif direction is "y":
                self.nodes[nid].set_disp(Uy = 0.)
                
        elif self.dim == 3:
            assert direction in ["x","y","z"],"Support dirction is x,y,and z"
            if direction is "x":
                self.nodes[nid].set_disp(Ux = 0.)
            
            elif direction is "y":
                self.nodes[nid].set_disp(Uy = 0.)
                
            elif direction is "z":
                self.nodes[nid].set_disp(Uz = 0.)

            
    def calc_deleted_KG_matrix(self):
        self._Force = [nd.force for nd in self.get_nodes()]
        self._Disp = [nd.disp for nd in self.get_nodes()]
        self._ForceValue = [val[key] for val in self.Force for key in self.nBk]
        self._DispValue = [val[key] for val in self.Disp for key in self.nAk]
        self._deleted = [row for row,val in enumerate(self.DispValue) if val is not None]
        self._keeped = [row for row,val in enumerate(self.DispValue) if val is None]
        if self._is_disp_added:
            self.check_boundary_condition(self.KG)
            
        self._Force_keeped = np.delete(self._ForceValue,self._deleted,0)
        self._KG_keeped = np.delete(np.delete(self._KG,self._deleted,0),self._deleted,1)


    def calc_deleted_MG_matrix(self):
        self._MG_keeped = np.delete(np.delete(self._MG,self._deleted,0),self._deleted,1)

    def check_boundary_condition(self,KG):
        self._nonzeros = [(row,val) for row,val in enumerate(self.DispValue) if val]
        if len(self.nonzeros):
            for i,val in self.nonzeros:
                for j in self.keeped:
                    self._ForceValue[j] -= KG[i,j]*val
        
    def check_deleted_KG_matrix(self):
        count = 0
        shape = self.KG_keeped.shape
        for i in range(shape[0]):
            if np.all(self.KG_keeped[i,:] == 0.):
                count += 1
        assert count == 0,"Check your bound conditions or system make sure it can be solved"

    def check_deleted_MG_matrix(self):
        count = 0
        shape = self.MG_keeped.shape
        for i in range(shape[0]):
            if np.all(self.MG_keeped[i,:] == 0.):
                count += 1
        assert count == 0,"Check your bound conditions or system make sure it can be solved"
        
    def solve(self,model = "static_elastic"):
        eval("solve"+"_"+model)(self)

             
    def results(self):
        self.postp = PostProcess(self.get_elements(),self.get_nodes(),self.dim)
        self.postp.results()

     

if __name__ == "__main__":
    pass
