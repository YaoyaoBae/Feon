# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

import numpy as np
def solve_simple(system):
    assert system._is_flowrate_added is True or system._is_head_added is True,"No flowrate or head added on the system"
    system.calc_deleted_KG_matrix()
    system.check_deleted_KG_matrix()
    KG,FlowRate = system.KG_keeped,system.FlowRate_keeped
    system._Head_keeped = np.linalg.solve(KG,FlowRate)


    for i,val in enumerate(system.keeped):
        I = val%system.mndof
        J = int(val/system.mndof)
        system.nodes[J].head[system.nAk[I]] = system.Head_keeped[i]

                
    for el in system.get_elements():
        el.evaluate()

    system._is_system_solved = True

    
