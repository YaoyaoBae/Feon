# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from feon.mesh import Mesh

if __name__ == "__main__":
    mesh = Mesh()
    mesh.build(mesh_type = "cube",x_lim = [0,10],y_lim = [0,5],z_lim = [0,5],size = [10,5,5])
    print mesh
    print mesh.points
    print mesh.elements
    
    
