# ReactionDiffusion
Python code to solve reaction-diffusion equation, with homogeneous and non-homogeneous diffusion coefficient, in evolving curves

The code allows solving reaction-diffusion equations in curves that grow and deform over time, with a non-homogeneous 
diffusion coefficient and Neumann (zero flow) boundary conditions. 

The code solves partial differential equations of the form

∂U/∂t = Dg(s, t) ∂/∂s[h(s, t)∂U/∂s]− a(s, t)U + γF(U)

With conditions
 
h(s, t)∂U/∂s=0, on the boundary

U(s,0)=Uo(s)

Where U=[u,v]^T, D is a 2x2 coefficient matrix, g(s,t), h(s,t) and a(s,t) are given functions, F(U) =[f1(u,v) f2(u,v)]^T are nonlinear functions corresponding 
to the chemical kinetics of the model and γ is a scale factor.

For the case of reaction-diffusion equations in deformable curves, we have following definitions:

g(s,t), h(s,t) y a(s,t)

g(s,t)=|| Xs(s,t)||^{-1}
h(s,t)= D(s) g(s,t)
a(s,t)= ∂/∂t [||X_s(s,t)||] g(s,t)

Where ||X_s(s,t)|| is the Euclidean norm of the tangent vectors Xs of the curve at position s and instant t, and D(s) is the non-homogeneous 
diffusion coefficient.

The system is solved using the spectral elements method for spatial discretization and the implicit-explicit Euler method for temporal discretization.


**Package contents**

* datos.py
* elementos.py
* lagrange.py
* legendre.py
* main2.py
* matg.py
* nodos_pesos.py

**Description**

* datos.py:       System parameters such as the definition of the boundary, the known functions, initial conditions, reaction kinetics, among others, 
                   are introduced.
* elementos.py:   Maps the different nodes (Legendre Gauss-Lobatto nodes) from a local domain (location in each element) to a global domain (location 
                   on the domain of the problem being worked on).
* lagrange.py:    Calculates the value of the Lagrange polynomial of degree N and its derivative at a point specified as input.
* legendre.py:    Calculates the value of the Legendre polynomials of degree N and its derivative at a point specified as input.
* matg.py:	  Calculates the matrix G at each element.
* nodos_pesos.py: Calculates the Legendre Gauss-Lobatto nodes and the Gauss Legendre-Lobatto quadrature weights.
* main2.py:       Main code where the previous codes are integrated. In addition, the spectral element method and the implicit-explicit Euler method 
                   are implemented for the solution of reaction-diffusion systems..

**Instructions**

* Open the file "datos.py".

* Go to the section "def datos_problema", where systems parameters are defined:
	- s_0: Left limit of s (left boundary).
	- s_1: Right limit of s (right boundary).
	-  N : Degree of the approximating polynomial (Lagrange). This value is unique for all elements.
	-  K : Number of domain divisions. This value together with s_0 and S_1 define the spatial step (d_s) or the width of each subdivision
                of the domain.
	- d_t: Initial time step. This value is an initial time step since internally the code modifies this value under a set error criterion.
	- t_f: Execution time. This value indicates the final execution time t of the computational simulation.
	- gam: Scale factor γ.
	-  d : Numerical value of the diffusion coefficient of the morphogen with concentration v. This value is the 22 entry of the diagonal matrix D.
	- d_1: Numerical value of the diffusion coefficient of the morphogen with concentration u. This value is the 11 entry of the diagonal matrix D.
	- capt_dat: Time values ​​where you want to capture the information corresponding to the numerical solution of the system at that time.

* Enter the reaction kinetics, these must be located in the "def func1(u,v)" section for morphogen with concentration u and in the "def func2(u,v)" 
   section for the morphogen with concentration v.

* Introduce the numerical values of the equilibrium point u_0 and v_0 in "def phi1_0" and "def phi2_0", respectively.
	-Although the amplitude of the disturbances around the equilibrium point is usually 10%, it is possible to modify this value in
           the variables "p" and "q". These disturbances represent the initial condition Uo(s) of the problem.

* In "def curve()", the information corresponding to the curve worked on in the problem is filled in, that is, the domain of interest of the problem.
      This information includes
	-  p : Domain growth function in the x coordinate.
	- p_1: Domain growth function in the y coordinate.
	- p_2: Domain growth function in the z coordinate.
	-  x : x value of the curve parametrization.
	-  y : y value of the curve parametrization.
	-  z : z value of the curve parametrization.
	- d_3: Non-homogeneous diffusion function D(s)for u and v morphogens.

* Open and run the file "main2.py".

NOTE: -For the parameterization coordinates we assume a curve of the form X=(x(s,t), y(s,t), z(s,t)). In the case of anisotropic growth
        we have that x(s,t)=p(t)*A(s), y(s,t)=p_1(t)*B(s), z(s,t)=p_2(t)*C(s). With A,B,C functions that define the parameterization of 
        the curve.
      -In the more general case X=(x(s,t), y(s,t), z(s,t)) it is not necessary to use the variables p, p_1 and p_2; it is enough to 
        define the functions x(s,t), (s,t) and z(s,t) in the variables x,y and z respectively.

