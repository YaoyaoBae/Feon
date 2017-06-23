# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from __future__ import division
import numpy as np

def tri_quad(func,vertex,method):
    xi = method.points.T
    x = np.multiply.outer(1.0 - xi[0] - xi[1], vertex[0]) + \
        np.multiply.outer(xi[0], vertex[1]) + \
        np.multiply.outer(xi[1], vertex[2])
    area = area_of_tri(vertex[0],vertex[1],vertex[2])
    w = method.weights
    return sum(func(v)*w[i]*area for i,v in enumerate(x))

def tetra_quad(func,vertex,method):
    xi = method.points.T
    x = np.multiply.outer(1.0 - xi[0] - xi[1] - xi[2], vertex[0]) + \
        np.multiply.outer(xi[0], vertex[1]) +\
        np.multiply.outer(xi[1], vertex[2]) +\
        np.multiply.outer(xi[2], vertex[3])
    x = x.T
    volume = volume_of_tetra(vertex)
    w = method.weights
    return sum(func(v)*w[i]*volume for i,v in enumerate(x))


def area_of_tri(p1,p2,p3):
    x,y = p1[0],p1[1]
    xj,yj = p2[0],p2[1]
    xm,ym = p3[0],p3[1]
    return 0.5*((xj*ym-yj*xm) +(yj-ym)*x +(xm-xj)*y)

def volume_of_tetra(vertex):
    V = np.ones((4,4))
    for i,vt in enumerate(vertex):
            V[i,1:] = vt
    return abs(np.linalg.det(V)/6.)



##----------------------------------------------------------------
#####################FOR TRI_QUAD#################################
##----------------------------------------------------------------
def _as3():
    return np.array([
        [1.0/3.0, 1.0/3.0, 1.0/3.0]
        ])


def _as21(a):
    b = 1.0 - 2*a
    return np.array([
        [a, a, b],
        [a, b, a],
        [b, a, a],
        ])


def _as111(a, b):
    c = 1.0 - a - b
    return np.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        [b, a, c],
        [c, b, a],
        [a, c, b],
        ])

class Hammer(object):
    def __init__(self,index = 5):
        self.name = 'HMS(%d)' % index
        
        if index == 1:
            self.weights = np.concatenate([
                1.0 * np.ones(1),
                ])
            bary = np.concatenate([
                _as3(),
                ])
            self.degree = 1
            
        elif index == 2:
            self.weights = np.concatenate([
                1.0/3.0 * np.ones(3),
                ])
            bary = np.concatenate([
                self._r(0.5),
                ])
            self.degree = 2
            
        elif index == 3:
            self.weights = np.concatenate([
                1.0/3.0 * np.ones(3),
                ])
            bary = np.concatenate([
                self._r(-0.5),
                ])
            self.degree = 2
            
        elif index == 4:
            self.weights = np.concatenate([
                -9.0/16.0 * np.ones(1),
                25.0/48.0 * np.ones(3),
                ])
            bary = np.concatenate([
                _as3(),
                self._r(0.4),
                ])
            self.degree = 3
            
        else:
            assert index == 5
            self.weights = np.concatenate([
                9.0/40.0 * np.ones(1),
                (155.0 - np.sqrt(15.0)) / 1200.0 * np.ones(3),
                (155.0 + np.sqrt(15.0)) / 1200.0 * np.ones(3),
                ])
            
            bary = np.concatenate([
                _as3(),
                self._r((1 + np.sqrt(15)) / 7.0),
                self._r((1 - np.sqrt(15)) / 7.0),
                ])
            self.degree = 5

        self.points = bary[:, 1:]

    def _r(self, r):
        a = r + (1.0-r) / 3.0
        b = 0.5 * (1.0 - a)
        return np.array([
            [a, b, b],
            [b, a, b],
            [b, b, a],
            ])

class Cubtri(object):
    def __init__(self):
        self.name = 'CUBTRI'
        self.weights = np.concatenate([
            0.0378610912003147 * np.ones(1),
            0.0376204254131829 * np.ones(3),
            0.0783573522441174 * np.ones(3),
            0.1162714796569659 * np.ones(3),
            0.0134442673751655 * np.ones(3),
            0.0375097224552317 * np.ones(6),
            ])

        bary = np.concatenate([
            _as3(),
            _as21(0.1012865073234563),
            _as21(0.4701420641051151),
            _as21(0.2321023267750504),
            _as21(0.0294808608844396),
            _as111(0.7384168123405100, 0.2321023267750504),
            ])
        self.points = bary[:, [1, 2]]
        self.degree = 8

        
class Triex(object):
    def __init__(self):
        self.name = 'TRIEX'
        self.weights = np.concatenate([
                0.08797730116222190 * np.ones(1),
                0.008744311553736190 * np.ones(3),
                0.03808157199393533 * np.ones(3),
                0.01885544805613125 * np.ones(3),
                0.07215969754474100 * np.ones(3),
                0.06932913870553720 * np.ones(3),
                0.04105631542928860 * np.ones(6),
                0.007362383783300573 * np.ones(6),
                ])
        bary = np.concatenate([
                _as3(),
                _as21(0.02598914092828833),
                _as21(0.09428750264792270),
                _as21(0.4946367750172147),
                _as21(0.2073433826145142),
                _as21(0.4389078057004907),
                _as111(0.6779376548825902, 0.04484167758913055),
                _as111(0.8588702812826364, 0.0),
                ])
        self.points = bary[:, [1, 2]]
        self.degree = 11
    

##----------------------------------------------------------------
#####################FOR TETRA_QUAD###############################
##----------------------------------------------------------------


def _bs4():
    return np.array([
        [0.25, 0.25, 0.25, 0.25]
        ])


def _bs31(a):
    b = 1.0 - 3*a
    return np.array([
        [a, a, a, b],
        [a, a, b, a],
        [a, b, a, a],
        [b, a, a, a],
        ])
def _bs22(a):
    b = 0.5 - a
    return np.array([
        [a, a, b, b],
        [a, b, a, b],
        [b, a, a, b],
        [a, b, b, a],
        [b, a, b, a],
        [b, b, a, a],
        ])


def _bs211(a, b):
    c = 1.0 - 2*a - b
    return np.array([
        [a, a, b, c],
        [a, b, a, c],
        [b, a, a, c],
        [a, b, c, a],
        [b, a, c, a],
        [b, c, a, a],
        [a, a, c, b],
        [a, c, a, b],
        [c, a, a, b],
        [a, c, b, a],
        [c, a, b, a],
        [c, b, a, a],
        ])


def _bs1111(a, b, c):
    d = 1.0 - a - b - c
    return np.array([
        [a, b, c, d],
        [a, b, d, c],
        [a, c, b, d],
        [a, c, d, b],
        [a, d, b, c],
        [a, d, c, b],
        [b, a, c, d],
        [b, a, d, c],
        [b, c, a, d],
        [b, c, d, a],
        [b, d, a, c],
        [b, d, c, a],
        [c, a, b, d],
        [c, a, d, b],
        [c, b, a, d],
        [c, b, d, a],
        [c, d, a, b],
        [c, d, b, a],
        [d, a, b, c],
        [d, a, c, b],
        [d, b, a, c],
        [d, b, c, a],
        [d, c, a, b],
        [d, c, b, a],
        ])


class Zctetra(object):
    def __init__(self):
        self.name = "Zienkiewicz"
        self.weights = np.concatenate([
                -0.8 * np.ones(1),
                0.45 * np.ones(4),
                ])
        bary = np.concatenate([
                _bs4(),
                _bs31(1.0/6.0),
                ])
        self.points = bary[:, [1, 2, 3]]
        self.degree = 3

        
if __name__ == "__main__":
    def f(x):
        return x[0]
    
    vertex = np.array([[0,0,0],[0,-1.,0.],[1.,-1.,0],[0,-1.,1.]])
    vt2 = np.array([[0,0],[1,0],[1,1]])
##    
##    import time
##    s = time.clock()
    a = Cubtri()
##    end = time.clock()
##    print end - s
##    s = time.clock()
    b = Triex()
##    d = Triex(19)
##    end = time.clock()
##    print end-s
##    s = time.clock()
    c = Hammer()
##    end = time.clock()
##    print end-s
    a = tetra_quad(f,vertex,Zctetra())
    b = tri_quad(f,vt2,Hammer(4))
    
