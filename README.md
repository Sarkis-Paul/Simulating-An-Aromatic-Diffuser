# Simulating & Evaluating the Performance of an Aromatic Diffuser


## Intoroduction

The aim of this project is to develop a model for the diffusion of an aromatic diffuser placed at a position 
(x*,y*,z*)=(1.5, 2, 2) in a 5mx5mx5m room. The quantity of interest is φ, the density of the diffuser mixture (g/cm3) in the room space over time modelled i.e. φ = φ(x,y,z,t) using the diffusion equation based upon the conservation of mass. The diffuser is modelled as a 4cmx4cmx4cm container with 75g of diffuser mixture that vapourises and diffuses out of the bottle over time (this is not a source, as the mass reduces over time).

A 3-dimensional space model is used over time (instead of 2D) despite the greater computation as it is more representative of the diffusion for the particular required case study.


## Partial Differential Equation (PDE) Diffusion Model

  (∂^2 φ)/(∂x^2 )+(∂^2 φ)/(∂y^2 )+  (∂^2 φ)/(∂z^2 )  =1/D   ∂φ/∂t          This is a parabollic PDE

D is the diffusivity of the diffuser in air 106 m2 s1. This is an estimation – more information can be found in section G.



## Boundary Conditions and Initial Values
Boundary Conditions

├ ∂φ/∂x┤|_(x=0)=0           ├ ∂φ/∂y┤|_(y=0)=0           ├ ∂φ/∂z┤|_(z=0)=0         Perfume molecules cannot diffuse through the boundary wall (flux is 0)

├ ∂φ/∂x┤|_(x=L)=0           ├ ∂φ/∂y┤|_(y=H)=0           ├ ∂φ/∂z┤|_(z=Z)=0         L, H and Z are the dimensions of the room (same reason as above)

Initial Condition

φ(x,y,z,0)=0           φ(x^*,y^*,z^*,0) = density of diffuser = mass/volume = 75/(4*4*4) = 1.1719 g/cm3




