# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from feon.sa import *
import numpy as np
##from math import pi, cos, sin
from meshpy.tet import MeshInfo, build
from meshpy.geometry import \
            generate_surface_of_revolution, EXT_OPEN, \
            GeometryBuilder
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from matplotlib import ticker
##import pyvtk

def write_vtk(filename,points,cells):
    vtkelements = pyvtk.VtkData(
            pyvtk.UnstructuredGrid(
                points,
                tetra=cells),
            "Mesh")
    vtkelements.tofile(filename)


def area_of_tri(nodes):
    v1 = np.array(nodes[0].coord[:2])-np.array(nodes[1].coord[:2])
    v2 = np.array(nodes[0].coord[:2])-np.array(nodes[2].coord[:2])
    area = np.cross(v1,v2)/2.
    return area

if __name__ == "__main__":
    #mesh
    r = 1
    l = 2
    rz = [(0,0), (r,0), (r,l), (0,l)]
    geob = GeometryBuilder()
    geob.add_geometry(*generate_surface_of_revolution(rz,
            radial_subdiv=20, ring_markers=[1,2,3]))
    mesh_info = MeshInfo()
    geob.set(mesh_info)
    mesh = build(mesh_info, max_volume=0.01)

   
    E = 210e6
    nu = 0.3
    P = 2000
    points = np.array(mesh.points)
    cells = np.array(mesh.elements)
    face_cells = np.array(mesh.faces)
    z_face0_cells = [c for c in face_cells if np.isclose(points[c][:,2],0).all()]
    z_face2_cells = [c for c in face_cells if np.isclose(points[c][:,2],2).all()]

    nds = np.array([Node(pt) for pt in points])
    els = [Tetra3D11(nds[c],E,nu) for c in cells]

    
    s = System()
    s.add_nodes(nds)
    s.add_elements(els)

    s.add_fixed_sup(*z_face0_cells)
    
    area = [area_of_tri(nds[c]) for c in z_face2_cells]
    s_area = sum(area)
    
    for i,cs in enumerate(z_face2_cells):
        co = area[i]/s_area/3.
        for c in cs:
            s.add_node_force(c,Fz = -P*co)
    s.solve()


##    #post_process
##    mesh.write_vtk("model.vtk")
##    face_nodes = [nds[c] for c in mesh.faces]
##    coords = np.zeros((len(face_nodes),3))
##    for i,nd in enumerate(face_nodes):
##        coords[i] = (nd.x + nd.disp["Ux"]*8e4,
##                     nd.y + nd.disp["Uy"]*8e4,
##                     nd.z + nd.disp["Uz"]*2e4)
##    write_vtk("disp.vtk",coords,mesh.faces)

    fig = plt.figure()
    ax1 = fig.add_subplot(211,aspect = "equal")
    ax2 = fig.add_subplot(212,aspect = "equal")

    face_nds = [nd for nd in nds if nd.z == 2]
    x = [nd.x for nd in face_nds]
    y = [nd.y for nd in face_nds]
    Ux = np.array([nd.disp["Ux"] for nd in face_nds])
    Uy = np.array([nd.disp["Uy"] for nd in face_nds])
    Uz = np.array([nd.disp["Uz"] for nd in face_nds])
    U = np.array([np.sqrt(nd.disp["Uz"]**2+nd.disp["Ux"]**2+nd.disp["Uy"]**2) for nd in face_nds])
    tr = tri.Triangulation(x,y)

    sfmt=ticker.ScalarFormatter(useMathText=True) 
    tp1 = ax1.tricontourf(tr,Ux,color = "k", cmap="jet")
    plt.colorbar(tp1,ax = ax1,ticks = np.linspace(min(Ux)*1e6,max(Ux)*1e6,6)*1e-6)
    tp2 = ax2.tricontourf(tr,Uy,color = "k", cmap="jet")
    plt.colorbar(tp2,ax = ax2,ticks = np.linspace(min(Ux)*1e6,max(Ux)*1e6,6)*1e-6)

    plt.show()
    
    
    
    

    
