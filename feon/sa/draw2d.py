# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import numpy as np

def draw_right_arrow(ax,node,inc,**kwargs):
    return ax.arrow(node.x,node.y,inc,0,**kwargs)

def draw_left_arrow(ax,node,inc,**kwargs):
    return ax.arrow(node.x,node.y,-inc,0.,**kwargs)

def draw_up_arrow(ax,node,inc,**kwargs):
    return ax.arrow(node.x,node.y,0.,inc,**kwargs)

def draw_down_arrow(ax,node,inc,**kwargs):
    return ax.arrow(node.x,node.y,0.,-inc,**kwargs)

def draw_element(ax,el,**kwargs):
    coords = np.array([nd.coord for nd in el.nodes])
    line = Line2D(coords[:,0],coords[:,1],**kwargs)
    line2 = Line2D((coords[:,0][-1],coords[:,0][0]),(coords[:,1][-1],coords[:,1][0]),**kwargs)
                    
    ax.add_line(line)
    ax.add_line(line2)

def draw_node_ID(ax,node,dx = 0.,dy = 0.,**kwargs):
    x,y = node.x,node.y
    ID = node.ID 
    ax.text(x+dx,y+dy,"$%d$"%ID,**kwargs)

def draw_element_ID(ax,el,dx = 0.,dy = 0.,**kwargs):
    coords = np.array([nd.coord for nd in el.nodes])
    n = len(coords)
    x,y = sum(coords[:,0])/n,sum(coords[:,1])/n
    ID = el.ID 
    ax.text(x+dx,y+dy,"$%d$"%ID,**kwargs)

def draw_element_disp(ax,el,factor = 1.,ls = "dashed",**kwargs):
    coord = np.array([nd.coord for nd in el.nodes])
    n = len(coord)
    disp = np.array([nd.disp[key] for nd in el.nodes for key in ["Ux","Uy"]])
    disp = disp.reshape(n,2)
    new_coords = coord+disp*factor*1e3
    line = Line2D(new_coords[:,0],new_coords[:,1],linestyle = ls,**kwargs)
    line2 = Line2D((new_coords[:,0][-1],new_coords[:,0][0]),(new_coords[:,1][-1],new_coords[:,1][0]),ls = ls,**kwargs)
    ax.add_line(line)
    ax.add_line(line2)
    
def draw_fixed_sup(ax,node,factor = (1,1),**kwargs):
    width = 0.1*factor[0]
    height = 0.1*factor[1]
    patch = mpatches.Rectangle((node.x - width * 0.5, node.y - height * 0.5),
                                           width, height,**kwargs)
    ax.add_patch(patch)

def draw_rolled_sup(ax,node,factor = 1,**kwargs):
    r = 0.1*factor
    patch = mpatches.Circle((node.x, node.y),radius = r/2.,**kwargs)
    ax.add_patch(patch)
    
def draw_hinged_sup(ax,node,factor = 1,**kwargs):
    r = 0.1*factor
    patch = mpatches.RegularPolygon((node.x, node.y),numVertices = 3,radius = r,**kwargs)
    ax.add_patch(patch)

    
def draw_spring(ax,el,**kwargs):
    x1,y1 = el.nodes[0].x,el.nodes[0].y
    x2,y2 = el.nodes[1].x,el.nodes[1].y
    p2x,p2y = (x1+x2)/2,(y1+y2)/2
    le = np.sqrt((x1-x2)**2+(y1-y2)**2)
    lx = (y2-y1)/le
    mx = (x2-x1)/le
    d = 1/5.*le
    dxy = d*(lx+mx)/2.
    p1x,p1y = p2x-d*lx,p2y+d*mx
    p3x,p3y = p2x+d*lx,p2y-d*mx
    ax.plot([x1,p1x,p2x,p3x,x2],[y1,p1y,p2y,p3y,y2],**kwargs)
    

def draw_bar_info(el):
    b = 1.0
    p_coord = [1.0*b,3.2*b]
    wd = 3.*b
    ht = 0.2*b

    r_y = 3.3*b #arrow y coord
        
    hl = 0.2*b
    hw = 0.08*b

    N_lx,N_ly = 0.5*b,3.5*b
    N_rx,N_ry = 4.5*b,3.5*b

    Ty_lx,Ty_ly = 0.5*b,3.1*b
    Ty_rx,Ty_ry = 4.5*b,3.1*b
        
    Mz_lx,Mz_ly = 1.0*b,3.7*b
    Mz_rx,Mz_ry = 4.0*b,3.7*b

    E_rx,E_ry = 2.5*b,3.0*b
    L_rx,L_ry = 2.5*b,2.7*b

    e_rx,e_ry = 2.5*b,2.4*b
    A_rx,A_ry = 2.5*b,2.1*b
    I_rx,I_ry = 2.5*b,1.8*b
    n1_rx,n1_ry = 2.5*b,1.5*b
    n2_rx,n2_ry = 2.5*b,1.2*b

        
    co1_rx,co1_ry = 1.0*b,2.8*b
    co2_rx,co2_ry = 4.0*b,2.8*b
    fig = plt.figure()
    ax = fig.add_subplot(111,aspect = "equal")
    ax.set_xlim(0,5*b)
    ax.set_ylim(0,5*b)
    ax.set_xticks([])
    ax.set_yticks([])
        
    coord1 = el.nodes[0].coord
    coord2 = el.nodes[1].coord
    N1,N2 = el.force["N"][0][0],el.force["N"][-1][0]
    if "Ty" in el.eIk:
        Ty1,Ty2 = el.force["Ty"][0][0],el.force["Ty"][0]
    else:
        Ty1,Ty2 =0.,0.
    if "Mz" in el.eIk:
            Mz1,Mz2 = el.force["Mz"][0][0],el.force["Mz"][1][0]
    else:
        Mz1,Mz2 = 0.,0.
    E = el.E
    A = el.A
    if hasattr(el,"I"):
        I = el.I
    else:
        I = 0.
        
    n1_Ux = el.nodes[0].disp["Ux"]
    n1_Uy = el.nodes[0].disp["Uy"]
    n1_Phz = el.nodes[0].disp["Phz"]
    n2_Ux = el.nodes[1].disp["Ux"]
    n2_Uy = el.nodes[1].disp["Uy"]
    n2_Phz = el.nodes[1].disp["Phz"]
    ID = el.ID
    volume = el.volume
    patch = mpatches.Rectangle(p_coord,width = wd,height = ht,
                                          hatch = "x",color = "y")
    ax.add_patch(patch)
        
    if N1 > 0:
        ax.plot(0.5*b,3.3*b,marker = r"$\longrightarrow$",ms = 40,color = "r")
    if N1 < 0:
        ax.plot(0.5*b,3.3*b,marker = r"$\longleftarrow$",ms = 40,color = "r")

    if N2 > 0:
        ax.plot(4.5*b,3.3*b,marker = r"$\longrightarrow$",ms = 40,color = "r")
    if N2 < 0:
        ax.plot(4.5*b,3.3*b,marker = r"$\longleftarrow$",ms = 40,color = "r")

    if Ty1 > 0:
        ax.plot(0.5*b,3.3*b,marker = r"$\uparrow$",ms = 40,color = "r")
    if Ty1 < 0:
        ax.plot(0.5*b,3.3*b,marker = r"$\downarrow$",ms = 40,color = "r")

    if Ty2 > 0:
        ax.plot(4.5*b,3.3*b,marker = r"$\uparrow$",ms = 40,color = "r")
    if Ty2 < 0:
        ax.plot(4.5*b,3.3*b,marker = r"$\downarrow$",ms = 40,color = "r")

        
    if Mz1 > 0:
        ax.plot(1.0*b,3.5*b,marker = r"$\curvearrowleft$",ms = 30,color = "r")
    if Mz1 < 0:
        ax.plot(1.0*b,3.5*b,marker = r"$\curvearrowright$",ms = 30,color = "r")

    if Mz2 > 0:
        ax.plot(4.0*b,3.5*b,marker = r"$\curvearrowright$",ms = 30,color = "r")
    if Mz2 < 0:
        ax.plot(4.0*b,3.5*b,marker = r"$\curvearrowleft$",ms = 30,color = "r")
            
    ax.text(N_lx,N_ly,"$N=%.2e$" % abs(N1),fontsize=12, zorder=9,ha="center", va="center")
    ax.text(N_rx,N_ry,"$N=%0.2e$" % abs(N2),fontsize=12, zorder=9,ha="center", va="center")
    ax.text(Ty_lx,Ty_ly,"$Ty=%.2e$" % abs(Ty1),fontsize=12, zorder=9,ha="center", va="center")
    ax.text(Ty_rx,Ty_ry,"$Ty=%.2e$" % abs(Ty2),fontsize=12, zorder=9,ha="center", va="center")
    ax.text(Mz_lx,Mz_ly,"$Mz=%.2e$" % abs(Mz1),fontsize=12, zorder=9,ha="center", va="center")
    ax.text(Mz_rx,Mz_ry,"$Mz=%.2e$" % abs(Mz2),fontsize=12, zorder=9,ha="center", va="center")
    ax.text(E_rx,E_ry,"$Element(%d)$" %ID,color='k',fontsize=12, zorder=9,ha="center", va="center")
    ax.text(L_rx,L_ry,"$Length = %.2e$" %volume,color='k',fontsize=12, zorder=9,ha="center", va="center")
    ax.text(e_rx,e_ry,"$E= %.2e$" %E,color='k',fontsize=12, zorder=9,ha="center", va="center")
    ax.text(A_rx,A_ry,"$A = %.2e$" %A,color='k',fontsize=12, zorder=9,ha="center", va="center")
    ax.text(I_rx,I_ry,"$I = %.2e$" %I,color='k',fontsize=12, zorder=9,ha="center", va="center")
    ax.text(n1_rx,n1_ry,"$Node1: Ux= %.2e,Uy= %.2e,Phz = %.2e$" %(n1_Ux,n1_Uy,n1_Phz),
                color='k',fontsize=12, zorder=9,ha="center", va="center")
    ax.text(n2_rx,n2_ry,"$Node2: Ux= %.2e,Uy= %.2e,Phz = %.2e$" %(n2_Ux,n2_Uy,n2_Phz),
                color='k',fontsize=12, zorder=9,ha="center", va="center")

    ax.text(co1_rx,co1_ry,"$%r$" %(coord1,),color='k',fontsize=12, zorder=9,ha="center", va="center")
    ax.text(co2_rx,co2_ry,"$%r$" %(coord2,),color='k',fontsize=12, zorder=9,ha="center", va="center")
    plt.show()  

if __name__ == "__main__":
    pass
