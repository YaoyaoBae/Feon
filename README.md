## Feon
Feon is a python-based finite element analysis framework for education and research purpose at the sponsor of Hubei Univeristy of Technology!
## Version
Current:1.0.0
## Requirements
Numpy is a must. Matplotlib is needed for visualization. Mpmath is need for Matrix derivation.
## Installation
```
Using pip:
$pip install feon
```
Or
```
python install setup.py
```
## Packages
* sa---For structrual analysis
* ffa --- For fluid flow analysis
* derivation --- matrix derivation 

## elements suported

* Spring1D11 
* Spring2D11
* Spring3D11

* Link1D11
* Link2D11
* Link3D11

* Beam1D11
* Beam2D11
* Beam3D11

* Tri2d11S---- Triange elements for plane stress problem
* Tri2D11 ---- Triange elements for plane strain problem
* Tetra3D11 
* Quad2D11S 
* Quad2D11
* Plate3D11 ---Midline plate
* Brick3D11

**We name the elements with "Name" + "dimensison" + 'order" + "type", type 1 means elastic.**

## Examples
![](examples/elements introduction/beam1/screenshot.PNG)


