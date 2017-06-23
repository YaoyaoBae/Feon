# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------
import numpy as np

class PostProcess(object):

    def __init__(self,els,nds,dim):
        self.dim = dim
        self.nodes = nds
        self.els = els

    def get_nodes_disp(self,key):
        return np.array([nd.disp[key] for nd in self.nodes if key in nd.nAk])

    def get_elements_force(self,key):
        info = []
        for el in self.els:
            if hasattr(el,"force"):
                if key in el.eIk:
                    info.append(el.force[key])
        return np.array(info)

    def get_elements_stress(self,key):
        info = []
        for el in self.els:
            if hasattr(el,"stress"):
                if key in el.eIk:
                    info.append(el.stress[key])
        return np.array(info)
        
    def get_max_sx(self):
        for el in self.els:
            if "sx" in el.eIk:
                sx = np.array([el.stress["sx"][0][0] for el in self.els])
                return sx[np.argmax(abs(sx))],np.argmax(abs(sx))
            else:
                return "nonexist","nonexist"

    def get_max_sy(self):
        for el in self.els:
            if "sy" in el.eIk:
                sy = np.array([el.stress["sy"][0][0] for el in self.els])
                return sy[np.argmax(abs(sy))],np.argmax(abs(sy))
            else:
                return "nonexist","nonexist"
            
    def get_max_sz(self):
        for el in self.els:
            if "sz" in el.eIk:
                sz = np.array([el.stress["sz"][0][0] for el in self.els])
                return sz[np.argmax(abs(sz))],np.argmax(abs(sz))
            else:
                return "nonexist","nonexist"
    def get_max_sxy(self):
        for el in self.els:
            if "sxy" in el.eIk:
                sxy = np.array([el.stress["sxy"][0][0] for el in self.els])
                return sxy[np.argmax(abs(sxy))],np.argmax(abs(sxy))
            else:
                return "nonexist","nonexist"

    def get_max_syz(self):
        for el in self.els:
            if "syz" in el.eIk:
                syz = np.array([el.stress["syz"][0][0] for el in self.els])
                return syz[np.argmax(abs(syz))],np.argmax(abs(syz))
            else:
                return "nonexist","nonexist"
            
    def get_max_szx(self):
        for el in self.els:
            if "szx" in el.eIk:
                szx = np.array([el.stress["szx"][0][0] for el in self.els])
                return szx[np.argmax(abs(szx))],np.argmax(abs(szx))
            else:
                return "nonexist","nonexist"
            
    def get_max_N(self):
        for el in self.els:
            if "N" in el.eIk:
                N = np.array([el.force["N"][0][0] for el in self.els])
                return N[np.argmax(abs(N))],np.argmax(abs(N))
            else:
                return "nonexist","nonexist"
            
                

    def get_max_Ty(self):
        for el in self.els:
            if "Ty" in el.eIk:
                Ty = np.array([el.force["Ty"][0][0] for el in self.els])
                return Ty[np.argmax(abs(Ty))],np.argmax(abs(Ty))
            else:
                return "nonexist","nonexist"
            
    def get_max_Tz(self):
        for el in self.els:
            if "Phx" in el.eIk:
                Tz = np.array([el.force["Tz"][0][0] for el in self.els])
                return Tz[np.argmax(abs(Tz))],np.argmax(abs(Tz))
            else:
                return "nonexist","nonexist"
            
    def get_max_Mz(self):
        for el in self.els:
            if "Mz" in el.eIk:
                Mz = np.array([el.force["Mz"][0][0] for el in self.els])
                return Mz[np.argmax(abs(Mz))],np.argmax(abs(Mz))
            else:
                return "nonexist","nonexist"

    def get_max_Mx(self):
        for el in self.els:
            if "Mx" in el.eIk:
                Mx = np.array([el.force["Mx"][0][0] for el in self.els])
                return Mx[np.argmax(abs(Mx))],np.argmax(abs(Mx))
            else:
                return "nonexist","nonexist"
            
    def get_max_My(self):
        for el in self.els:
            if "My" in el.eIk:
                My = np.array([el.force["My"][0][0] for el in self.els])
                return My[np.argmax(abs(My))],np.argmax(abs(My))
            else:
                return "nonexist","nonexist"
            
    def get_max_Ux(self):
        Ux = np.array([nd.disp["Ux"] for nd in self.nodes])
        return Ux[np.argmax(abs(Ux))],np.argmax(abs(Ux))

    def get_max_Uy(self):
        Uy = np.array([nd.disp["Uy"] for nd in self.nodes])
        return Uy[np.argmax(abs(Uy))],np.argmax(abs(Uy))
    
    def get_max_Uz(self):
        if self.dim == 3:
            Uz = np.array([nd.disp["Uz"] for nd in self.nodes])
            return Uz[np.argmax(abs(Uz))],np.argmax(abs(Uz))
        else:
            return "nonexist","nonexist"

    def get_max_Phx(self):
        if self.dim == 3:
            Phx = np.array([nd.disp["Phx"] for nd in self.nodes])
            return Phx[np.argmax(abs(Phx))],np.argmax(abs(Phx))
        else:
            return "nonexist","nonexist"
    
    def get_max_Phy(self):
        if self.dim == 3:
            Phy = np.array([nd.disp["Phy"] for nd in self.nodes])
            return Phy[np.argmax(abs(Phy))],np.argmax(abs(Phy))
        else:
            return "nonexist","nonexist"

    
    def get_max_Phz(self):
        Phz = np.array([nd.disp["Phz"] for nd in self.nodes])
        return Phz[np.argmax(abs(Phz))],np.argmax(abs(Phz))
            
    def get_max_disp(self):
        if self.dim == 2:
            m = [np.sqrt(nd.disp["Ux"]**2 + nd.disp["Uy"]**2)
                 for nd in self.nodes]
            max_v = max(m)
            ID = m.index(max_v)
            return max_v,ID
        if self.dim == 3:
            m = [np.sqrt(nd.disp["Ux"]**2 + nd.disp["Uy"]**2 + nd.disp["Uz"]**2)
                    for nd in self.nodes]
            max_v = max(m)
            ID = m.index(max_v)

        return max_v,ID

    def results(self):
        results = Results.format(
            st = str(self.dim) + "D System",
            n = len(self.nodes),
            
            m = len(self.els),
            max_sx_ID = self.get_max_sx()[1],
            max_sx = self.get_max_sx()[0],
            
            max_sy_ID = self.get_max_sy()[1],
            max_sy = self.get_max_sy()[0],
            
            max_sz_ID = self.get_max_sz()[1],
            max_sz = self.get_max_sz()[0],
            
            max_sxy_ID = self.get_max_sxy()[1],
            max_sxy = self.get_max_sxy()[0],
            
            max_syz_ID = self.get_max_syz()[1],
            max_syz = self.get_max_syz()[0],
            
            max_szx_ID = self.get_max_szx()[1],
            max_szx = self.get_max_szx()[0],
            
            max_N_ID = self.get_max_N()[1],
            max_N = self.get_max_N()[0],
            
            max_Ty_ID = self.get_max_Ty()[1],
            max_Ty = self.get_max_Ty()[0],
            
            max_Tz_ID = self.get_max_Tz()[1],
            max_Tz = self.get_max_Tz()[0],
            
            max_Mx_ID = self.get_max_Mx()[1],
            max_Mx = self.get_max_Mx()[0],
            
            max_My_ID = self.get_max_My()[1],
            max_My = self.get_max_My()[0],
            
            max_Mz_ID = self.get_max_Mz()[1],
            max_Mz = self.get_max_Mz()[0],
            
            max_Ux_ID = self.get_max_Ux()[1],
            max_Ux = self.get_max_Ux()[0],

            max_Uy_ID = self.get_max_Uy()[1],
            max_Uy = self.get_max_Uy()[0],

            max_Uz_ID = self.get_max_Uz()[1],
            max_Uz = self.get_max_Uz()[0],
            
            max_Phx_ID = self.get_max_Phx()[1],
            max_Phx = self.get_max_Phx()[0],

            max_Phy_ID = self.get_max_Phy()[1],
            max_Phy = self.get_max_Phy()[0],
            
            max_Phz_ID = self.get_max_Phz()[1],
            max_Phz = self.get_max_Phz()[0],
            
            max_disp_ID = self.get_max_disp()[1],
            max_disp = self.get_max_disp()[0])
        
        print (results)
        
            
Results = """
==========================
         Results
==========================
Type: {st}
Number of nodes: {n}
Number of elements: {m}

Max element sx ID: {max_sx_ID}
Max element sx: {max_sx}

Max element sy ID: {max_sy_ID}
Max element sy: {max_sy}

Max element sz ID: {max_sz_ID}
Max element sz: {max_sz}

Max element sxy ID: {max_sxy_ID}
Max element sxy: {max_sxy}

Max element syz ID: {max_syz_ID}
Max element syz: {max_syz}

Max element szx ID: {max_szx_ID}
Max element szx: {max_szx}

Max element N ID: {max_N_ID}
Max element N: {max_N}

Max element Ty ID:{max_Ty_ID}
Max element Ty:{max_Ty}

Max element Tz ID: {max_Tz_ID}
Max element Tz: {max_Tz}

Max element Mx ID: {max_Mx_ID}
Max element Mx: {max_Mx}

Max element My ID: {max_My_ID}
Max element My:{max_My}

Max element Mz ID: {max_Mz_ID}
Max element Mz: {max_Mz}

Max node Ux ID: {max_Ux_ID}
Max node Ux: {max_Ux}

Max node Uy ID: {max_Uy_ID}
Max node Uy: {max_Uy}

Max node Uz ID: {max_Uz_ID}
Max node Uz: {max_Uz}

Max node Phx ID: {max_Phx_ID}
Max node Phx: {max_Phx}

Max node Phy ID: {max_Phy_ID}
Max node Phy: {max_Phy}

Max node Phz ID: {max_Phz_ID}
Max node Phz: {max_Phz}

Max node disp ID: {max_disp_ID}
Max node disp: {max_disp}

"""
