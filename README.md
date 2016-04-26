# Method

This algorithm implements the Fourier transform of Maxwell’s equations in tensor form, allowing for the forward calculation of induced field from a cylinder of known susceptibility in 3D. 

In MRI, the Fourier transform of Maxwell’s Equations (Marques and Bowtell, 2005 and Salomir et al., 2003) allows for the conversion between distributions of magnetic susceptibility and field maps. While most structures in human tissue exhibit isotropic magnetic susceptibility, some biomolecules have been suggested to have magnetic susceptibility anisotropy - for example, myelin in the brain. This anisotropy can be characterized by a rank-2 tensor. The Fourier transform of Maxwell’s equations accommodating this tensor formulation of magnetic susceptibility is derived in Equation 9 of Liu 2010, *Susceptibility Tensor Imaging*. This routine incorporates the tensor formulation of magnetic susceptibility in the forward calculations of field distortions from a solid cylinder in 3D. 

# Input & Output 
function [A]= sus_tensor(PHI,THETA,B0,sus, AA,BB,CC)

**Inputs** (in spherical coordinates where applicable):  
Phi, 0 --> pi/2, azimuthal angle of the applied field to the cylinder  
Theta, 0 --> pi/2, orientation or polar angle of the applied field to longitudinal axis of the cylinder  
B0 - applied field strength, Tesla  
Sus - magnetic susceptibility in ppb  
AA,BB,CC - diagonals entries to a 3x3 matrix, denoting the susceptibility tensor  
For example, [1,1,1] denotes isotropic susceptibility | [1,-½,-½ ] denotes radial anisotropy of myelin. 

**Outputs**  
Field perturbation from cylinder in 3D. Output image is an orthogonal cross section about longitudinal axis of the cylinder. 

# Results
In result.png  
A) cross section of cylinder  
B) sus_tensor(0,pi/2,7,5,1,-½,-½)  
C) sus_tensor(pi/4,pi/2,7,5,1,-½,-½)  
D) sus_tensor(pi/2,pi/2,7,5,1,-½,-½)
