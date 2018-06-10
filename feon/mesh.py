# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

import numpy as np
from .tools import pair_wise
class Mesh(object):

    def __init__(self):
        self.dim = 2
        self.mesh_type = ""
        self.x_lim = None
        self.y_lim = None
        self.z_lim = None

        self.points = None       
        self.p_num = 0
        
        self.elements = None   
        self.e_num = 0

    def build(self, **params):
        for key in params.keys():
            assert key in ['mesh_type',"x_lim","y_lim","z_lim","size"],"unknow value key/keys"
        self.mesh_type = params['mesh_type']
        self.x_lim = params['x_lim']
        self.y_lim = params['y_lim']
        if "z_lim" in params:
            self.dim = 3
            self.z_lim = params['z_lim']
            
        if self.mesh_type is 'rect':
            self.points, self.elements = create_rect(self.x_lim, self.y_lim, 
                                             params['size'])
        elif self.mesh_type is 'tri_from_rect':
            self.points, self.elements = create_tri_from_rect(self.x_lim, self.y_lim, 
                                                        params['size'])
        elif self.mesh_type is "cube":
            self.points, self.elements = create_cube(self.x_lim, self.y_lim,
                                                     self.z_lim,params['size'])
        else:
            raise AttributeError("unknown mesh type")
        self.p_num = self.points.shape[0]
        self.e_num = self.elements.shape[0]

    def write_vtk(self, filename):
        import pyvtk
        if self.mesh_type is "cube":
            vtkelements = pyvtk.VtkData(
                pyvtk.UnstructuredGrid(
                    self.points,
                    hexahedron=self.elements),
                "Mesh")
            vtkelements.tofile(filename)
        

    def __repr__(self):
        return show_info(self)

def create_rect(x_lim, y_lim, size):
    nx = int(size[0])
    ny = int(size[1])
    X = np.linspace(x_lim[0],x_lim[1],nx+1)
    Y = np.linspace(y_lim[0],y_lim[1],ny+1)
    p = np.array([(i,j) for i in X for j in Y])
    e_cell = np.array([((ny+1)*i[0]+j[0],
                        (ny+1)*i[1]+j[0],
                        (ny+1)*i[1]+j[1],
                        (ny+1)*i[0]+j[1])
                       for i in pair_wise(range(nx+1)) for j in pair_wise(range(ny+1))],dtype = int)
    return p, e_cell

def create_tri_from_rect(x_lim, y_lim, size):
    nx = int(size[0])
    ny = int(size[1])
    X = np.linspace(x_lim[0],x_lim[1],nx+1)
    Y = np.linspace(y_lim[0],y_lim[1],ny+1)
    p = np.array([(i,j) for i in X for j in Y])
    e_cell = np.array([((ny+1)*i[0]+j[0],
                        (ny+1)*i[1]+j[0],
                        (ny+1)*i[1]+j[1],
                        (ny+1)*i[0]+j[0],
                        (ny+1)*i[1]+j[1],
                        (ny+1)*i[0]+j[1])
                       for i in pair_wise(range(nx+1)) for j in pair_wise(range(ny+1))],dtype = int)
    return p, e_cell.reshape(nx*ny*2,3)

    

def create_cube(x_lim, y_lim,z_lim,size):
    nx = int(size[0])
    ny = int(size[1])
    nz = int(size[2])
    X = np.linspace(x_lim[0],x_lim[1],nx+1)
    Y = np.linspace(y_lim[0],y_lim[1],ny+1)
    Z = np.linspace(z_lim[0],z_lim[1],nz+1)
    p = np.array([(i,j,k) for i in X for j in Y for k in Z])
    e_cell = np.array([((nz+1)*(ny+1)*i[0]+(nz+1)*j[0]+k[0],
                        (nz+1)*(ny+1)*i[0]+(nz+1)*j[1]+k[0],
                        (nz+1)*(ny+1)*i[0]+(nz+1)*j[1]+k[1],
                        (nz+1)*(ny+1)*i[0]+(nz+1)*j[0]+k[1],
                        (nz+1)*(ny+1)*i[1]+(nz+1)*j[0]+k[0],
                        (nz+1)*(ny+1)*i[1]+(nz+1)*j[1]+k[0],
                        (nz+1)*(ny+1)*i[1]+(nz+1)*j[1]+k[1],
                        (nz+1)*(ny+1)*i[1]+(nz+1)*j[0]+k[1],)
                       for i in pair_wise(range(nx+1))
                       for j in pair_wise(range(ny+1))
                       for k in pair_wise(range(nz+1))],dtype = int)
    return p, e_cell


def show_info(Mesh):
    s = '------------------------------------'
    s+= '\n Mesh description  (type: ' + Mesh.mesh_type + ')'
    s+= '\n X / Y/ Z limits: '
    s+= Mesh.x_lim.__str__() +'/' + Mesh.y_lim.__str__()+'/'+Mesh.z_lim.__str__()
    s+= '\n Number of points  :  %d'
    s = s%Mesh.p_num
    s+= '\n Number of elements  :  %d'
    s = s%Mesh.e_num
    s+='\n------------------------------------'
    return s

if __name__ == "__main__":
    mesh = Mesh()
    mesh.build(mesh_type = "cube",x_lim = [0,10],y_lim = [0,2],z_lim = [0,5],size = [10,2,5])
