# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from feon.mesh import Mesh

if __name__ == "__main__":
    mesh = Mesh()
    mesh.build(mesh_type = "rect",x_lim = [0,2],y_lim = [0,2],size = [2,2])
    print mesh
    print mesh.points
    print mesh.elements
    
    
