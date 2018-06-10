# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

import numpy as np


#~ *********************************************************************
#~ **************************  CLASS NODE ******************************
#~ *********************************************************************
__all__ = ["NodeBase","ElementBase","SystemBase"]

class NodeBase(object):
    
    def __init__(self,*coord):
        self.init_coord(*coord)
        self._ID = None
        self._nAk = None
        self._nBk = None
        self.init_keys()
        
        
    def __repr__(self):
        return "Node:%r"%(self.coord,)
    
    def __getitem__(self,key):
        return self.coord[key]

    def __setitem__(self,key,val):
        l = self.coord
        l[key] = val
        self.coord = l

    def __eq__(self,other):
        assert issubclass(type(other),NodeBase),"Must be Node type"
        a = np.array(self.coord)
        b = np.array(other.coord)
        if np.isclose(a,b).all():
            return True
        else:
            return False
        
    def init_coord(self,*coord):
        self._x = 0.
        self._y = 0.
        self._z = 0.
        
        if len(coord) == 1:
            self.dim = len(coord[0])
            if self.dim == 2:
                self._x = coord[0][0]*1.
                self._y = coord[0][1]*1.
                self.coord = (self.x,self.y)
            elif self.dim == 3:
                self._x = coord[0][0]*1.
                self._y = coord[0][1]*1.
                self._z = coord[0][2]*1.
                self.coord = (self.x,self.y,self.z)
            else:
                raise AttributeError("Node dimension is 2 or 3")
            
        elif len(coord) == 2:
            self.coord = tuple(coord)
            self.dim = 2
            self._x = coord[0]*1.
            self._y = coord[1]*1.
            self.coord = (self.x,self.y)
        elif len(coord) == 3:
            self.coord = tuple(coord)
            self.dim = 3
            self._x = coord[0]*1.
            self._y = coord[1]*1.
            self._z = coord[2]*1.
            self.coord = (self.x,self.y,self.z)
        else:
            raise AttributeError("Node dimension is 2 or 3")


    @property
    def nBk(self):
        return self._nBk

    @property
    def nAk(self):
        return self._nAk

    def init_unknowns(self):
        pass

    def init_keys(self):
        pass
        
    def set_nAk(self,val):
        self._nAk = val

    def get_nAk(self):
        return self._nAk

    def set_nBk(self,val):
        self._nBk = val

    def get_nBk(self):
        return self._nBk
    
    @property
    def x(self):
        return self._x
    
    @x.setter
    def x(self,val):
        self._x = val

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self,val):
        self._y = val

    @property
    def z(self):
        return self._z

    @z.setter
    def z(self,val):
        self._z = val

    @property
    def ID(self):
        return self._ID

    @ID.setter
    def ID(self,val):
        self._ID = val
    
        


#~ *********************************************************************
#~ ***********************  CLASS ELEMENT ******************************
#~ *********************************************************************


class ElementBase(object):

    def __init__(self,nodes):
        for nd in nodes:
            assert issubclass(type(nd),NodeBase),"Must be Node type"
            
        self._nodes = nodes
        self._dim = 0
        self._etype = self.__class__.__name__
        self._ID = None
        self._eIk = None
        self._ndof = None #node dof
        self._Ke = None
        self._Me = None
        self._volume = None
        self.init_keys()
        
        
        
    def __repr__(self):
        return "%s Element: %r"%(self.etype,self.nodes,)
    
    def __getitem__(self,key):
        return self._nodes[key]

    @property
    def volume(self):
        return self._volume

    @property
    def B(self):
        return self._B

    @property
    def D(self):
        return self._D

    @property
    def J(self):
        return self._J
    
    @property
    def dim(self):
        return self.nodes[0].dim
    
    @property
    def eIk(self):
        return self._eIk
    
    @property
    def Ke(self):
        return self._Ke

    @property
    def Me(self):
        return self._Me
    
    @property
    def ndof(self):
        return self._ndof
    
        
    @property
    def nodes(self):
        return self._nodes
    
    @property
    def etype(self):
        return self._etype
        
    @property
    def ID(self):
        return self._ID
    

    @property
    def non(self):
        return len(self._nodes)

    def init_keys(self):
        pass

    def func(self,x):
        pass

    def func_jac(self,x):
        pass


    def calc_Ke(self):
        pass
    
    def set_eIk(self,val):
        self._eIk = val

    def get_eIk(self):
        return self._eIk

    def set_ndof(self,val):
        self._ndof = val

    def get_ndof(self):
        return self._ndof
    
    def get_element_type(self):
        return self._etype

    def get_nodes(self):
        return self._nodes





#~ *********************************************************************
#~ **************************  CLASS SYSTEM ****************************
#~ *********************************************************************
    
class SystemBase(object):

    def __init__(self):
        self.nodes = {}
        self.elements = {}
        self._mndof = None
        self._nAk = None
        self._nBk = None
        self._dim = 0
        
    @property
    def mndof(self): #max length of keys
        return self._mndof
    
    @property
    def dim(self):
        return self._dim
    
    @property
    def nAk(self):
        return self._nAk
    @property
    def nBk(self):
        return self._nBk
    @property
    def non(self): #number of nodes: non
        return len(self.nodes)

    @property
    def noe(self): #number of elements
        return len(self.elements)
    
    def add_node(self,node):
        assert issubclass(type(node),NodeBase),"Must be Node type"
        n = self.non
        if node.ID is None:
            node._ID = n
        self.nodes[node.ID] = node

    def add_element(self,element):
        assert issubclass(type(element),ElementBase),"Must be Element type"
        n = self.noe
        if element.ID is None:
            element._ID = n
        self.elements[element.ID] = element

    def add_nodes(self,*nodes):
        for nd in nodes:
            if isinstance(nd,list) or isinstance(nd,tuple) or isinstance(nd,np.ndarray):
                for n in nd:
                    self.add_node(n)
            else:
                self.add_node(nd)
                
    def add_elements(self,*els):
        for el in els:
            if isinstance(el,list) or isinstance(el,tuple) or isinstance(el,np.ndarray):
                for e in el:
                    self.add_element(e)
            else:
                self.add_element(el)

    def init(self):
        pass

    def calc_KG(self):
        pass

    def calc_MG(self):
        pass
                
    def get_nodes(self):
        return self.nodes.values()

    def get_elements(self):
        return self.elements.values()
        


    
if __name__ == "__main__":
    n = NodeBase(1,2,3)
    n2 = NodeBase(1,2,3)
    e = ElementBase((n,n2))
    
    
