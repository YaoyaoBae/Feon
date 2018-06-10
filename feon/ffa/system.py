# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -----------------------------------


import numpy as np
from .solver import *
from ..base import SystemBase
##from draw2d import Drawer
#~ *********************************************************************
#~ ****************************  CLASS SYSTEM **************************
#~ *********************************************************************

class System(SystemBase):

    def __init__(self):
        SystemBase.__init__(self)
        self._FlowRate = {}
        self._Head = {}        
        self._is_inited = False
        self._is_flowrate_added = False
        self._is_head_added = False
        self._is_system_solved = False
    

    def __repr__(self):
        return "%dD System: \nNodes: %d\nElements: %d"\
               %(self.dim,self.non,self.noe,)
    
    @property
    def FlowRate(self):
        return self._FlowRate

    @property
    def Head(self):
        return self._Head

    @property
    def HeadValue(self):
        return self._HeadValue

    @property
    def FlowRateValue(self):
        return self._FlowRateValue
    
    @property
    def KG(self):
        return self._KG
    
    @property
    def KG_keeped(self):
        return self._KG_keeped


    @property
    def FlowRate_keeped(self):
        return self._FlowRate_keeped

    @property
    def Head_keeped(self):
        return self._Head_keeped

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
        self._mndof = 1
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
            M = 1
            for N1,I in enumerate(ID):
                for N2,J in enumerate(ID):
                    self._KG[m*I:m*I+M,m*J:m*J+M] += el.Ke[M*N1:M*(N1+1),M*N2:M*(N2+1)]

        self._is_inited = True
        
        
    def add_node_flowrate(self,nid,**flowrate):
        if not self._is_inited:
            self.calc_KG()
            
        assert nid+1 <= self.non,"Element does not exist"
        for key in flowrate.keys():
            assert key in self.nBk,"Check if the node flow rate applied are correct"

        self.nodes[nid].set_flowrate(**flowrate)   
        self._is_flowrate_added = True
        

    def add_node_head(self,nid,head):
        if not self._is_inited:
            self.calc_KG()
        assert nid+1 <= self.non,"Element does not exist"
        
        self.nodes[nid].set_head(head)
        if head:
            self._is_head_added = True

            
            
    def calc_deleted_KG_matrix(self):
        self._FlowRate = [nd.flowrate for nd in self.get_nodes()]
        self._Head = [nd.head for nd in self.get_nodes()]
        self._FlowRateValue = [val[key] for val in self.FlowRate for key in self.nBk]
        self._HeadValue = [val[key] for val in self.Head for key in self.nAk]
        self._deleted = [row for row,val in enumerate(self.HeadValue) if val is not None]
        self._keeped = [row for row,val in enumerate(self.HeadValue) if val is None]
        if self._is_head_added:
            self.check_boundary_condition(self.KG)
        self._FlowRate_keeped = np.delete(self._FlowRateValue,self._deleted,0)
        self._KG_keeped = np.delete(np.delete(self._KG,self._deleted,0),self._deleted,1)


    def check_boundary_condition(self,KG):
        self._nonzeros = [(row,val) for row,val in enumerate(self.HeadValue) if val]
        if len(self.nonzeros):
            for i,val in self.nonzeros:
                for j in self.keeped:
                    self._FlowRateValue[j] -= KG[i,j]*val
        
    def check_deleted_KG_matrix(self):
        count = 0
        shape = self.KG_keeped.shape
        for i in range(shape[0]):
            if np.all(self.KG_keeped[i,:] == 0.):
                count += 1
        assert count == 0,"Check your bound conditions or system make sure it can be solved"

        
    def solve(self,model = "simple"):
        eval("solve"+"_"+model)(self)


if __name__ == "__main__":
    pass
