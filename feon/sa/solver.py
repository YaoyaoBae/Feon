# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

import numpy as np

def solve_static_elastic(system):
    assert system._is_force_added is True or system._is_disp_added is True,"No forces or disp on the structure"
    system.calc_deleted_KG_matrix()
    system.check_deleted_KG_matrix()
    KG,Force = system.KG_keeped,system.Force_keeped
    system._Disp_keeped = np.linalg.solve(KG,Force)


    for i,val in enumerate(system.keeped):
        I = val%system.mndof
        J = int(val/system.mndof)
        system.nodes[J].disp[system.nAk[I]] = system.Disp_keeped[i]

                
    for el in system.get_elements():
        el.evaluate()

    system._is_system_solved = True

def solve_dynamic_eigen_model(system):
    from scipy import linalg as sl
    if not system._is_inited:
        system.calc_KG()
    system.calc_deleted_KG_matrix()
    system.calc_MG()
    system.calc_deleted_MG_matrix()
    
    w1,system.model = sl.eigh(system.KG_keeped,system.MG_keeped)
    system.w = np.sqrt(w1)
    T = 2*np.pi/system.w
    system.freq = 1/T
    
