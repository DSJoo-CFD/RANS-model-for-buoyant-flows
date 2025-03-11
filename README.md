<header>

<!--
  <<< Author notes: Course header >>>
  Include a 1280×640 image, course title in sentence case, and a concise description in emphasis.
  In your repository settings: enable template repository, add your 1280×640 social image, auto delete head branches.
  Add your open source license, GitHub uses MIT license.
-->

The current code is built on OpenFOAM-v2406. 
Below is a brief guide on how to use it.

## 1. Install OpenFOAM 
Follow the guidelines provided at [https://www.openfoam.com/] to install OpenFOAM. 

Note: This code is not compatible with the OpenFOAM Foundation release (often labeled as v8, ..., v12) available at [openfoam.org]

## 2. Set Up Your User Directory
After installation, update your personal "~/.bashrc" file to include OpenFOAM environment variables refer to the installation guidelines. 

Then, create the necessary directories for running OpenFOAM by executing:

  > mkdir -p $WM_PROJECT_USER_DIR/{run,applications,src}

Next, overwrite the provided source codes by copying the provided files into the appropriate location.

## 3. Compile the Solver
In a Linux environment, compile the solver located in "applications/mybuoyantBoussinesqSimpleFoam" by running the following commands sequentially:

  > app

  > cd mybuoyantBoussinesqSimpleFoam
  
  > wclean
  
  > wmake

## 4. Run the Examples  
The "run" folder contains the problems validated in the paper. 

You can run these examples just like any conventional OpenFOAM tutorial. 

Simply navigate to the example directory and execute:
  
  > run
  
  > cd RayleighBenardConvection
  
  > sh ./Allclean
  
  > sh ./Allrun
 


# Model equations
The following turbulence models are implemented together:

1. Standard $k$ - $\varepsilon$ RANS model by Launder and Sharma (1974).
2. Algebraic turbulent heat flux model for buoyancy effects by Kenjereš and Hanjalić (1995, Int. J. Heat Fluid Flow, 16(5)).


All Reynolds-averaged equations for incompressible flow are given below.

$\dfrac{D U_i}{D t}	= - \dfrac{\partial P}{\partial x_i} 	- g_i \beta \left( T - T_0 \right)	+ \dfrac{\partial }{\partial x_j} 	\left( \nu \dfrac{\partial U_i}{\partial x_j} -\overline{u_i u_j} \right),$

$\dfrac{D T}{D t}
	=
	Q
	+ \dfrac{\partial }{\partial x_i} 
	\left( \alpha \dfrac{\partial T}{\partial x_i} - \overline{\theta u_i} \right).$

Here, $U_i$ and $u_i$ denote the mean and fluctuating velocities from the Reynolds decomposition, respectively, while $T$ and $\theta$ represent the mean and fluctuating temperatures. The overbar indicates Reynolds averaging, so $\overline{u_i u_j}$ and $\overline{\theta u_i}$ represent the Reynolds stress ($\mathrm{m^2/s^2}$) and the turbulent heat flux ($\mathrm{K \cdot m/s}$), respectively. Additional variables include the gravitational acceleration $g_i$ ($\mathrm{m/s^2}$) and the thermal expansion coefficient $\beta$ ($\mathrm{K^{-1}}$). The mean pressure per unit density is denoted by $P$. $T_0$ represents the reference temperature in the Boussinesq approximation. The rate of temperature increase due to a volumetric heat source is given by $Q$ ($\mathrm{K/s}$), and $\alpha$ denotes the thermal diffusivity.  

$\dfrac{D k}{D t} = -\overline{u_i u_j} \dfrac{\partial U_i}{\partial x_j} 
	- g_i \beta \overline{\theta u_i}
	- \left( \varepsilon + \varepsilon_0 \right) 
	+ \dfrac{\partial }{\partial x_i} 
	\left( \left(\nu + \dfrac{\nu_t}{\sigma_k}\right) \dfrac{\partial k}{\partial x_i}  \right),$

 $\dfrac{D \varepsilon}{D t} =
	C_{\varepsilon1}\dfrac{\varepsilon}{k} \left(-\overline{u_i u_j} \dfrac{\partial U_i}{\partial x_j} \right)
	+ C_{\varepsilon g}  \dfrac{\varepsilon}{k} \left(- g_i \beta \overline{\theta u_i}\right)
	- C_{\varepsilon2} f_\varepsilon \dfrac{\varepsilon^{2}}{k}
	+ \dfrac{\partial }{\partial x_i} 
	\left( \left(\nu + \dfrac{\nu_t}{\sigma_\varepsilon}\right) \dfrac{\partial \varepsilon}{\partial x_i}  \right),$
 
 $\varepsilon_0 = 2\nu \left( \dfrac{\partial k^{1/2}}{\partial x_i} \dfrac{\partial k^{1/2}}{\partial x_i} \right), \quad 
	\nu_t = C_{\mu} f_{\mu} \dfrac{k^{2}}{\varepsilon}, \quad
	\overline{u_{i} u_{j}}
	=
	\dfrac{2}{3}k \delta_{ij} - \nu_{t} \left( \dfrac{\partial U_i}{\partial x_j} + \dfrac{\partial U_j}{\partial x_i} \right).$
 
Here, $\nu_t$ is the turbulent viscosity, $k$ is the turbulent kinetic energy ($\mathrm{m^2/s^2}$), and $\varepsilon$ is the dissipation rate ($\mathrm{m^2/s^3}$)

The damping functions are 

$f_{\varepsilon} = 1.0 - 0.3 \exp(-\mathit{Re}_{t}^{2})$

$f_{\mu} = \exp \left( -3.4 / \left( 1+ \mathit{Re}_{t}/50 \right)^{2}\right)$

with 

$\mathit{Re}_{t} = {k^{2}}/({\nu \varepsilon})$.

The standard coefficients are $C_{\varepsilon1} = 1.44$, $C_{\varepsilon2} = 1.92$,	$C_{\mu} = 0.09$, $\sigma_k = 1.0$, $\sigma_\varepsilon = 1.3$.

For the buoyancy-related coefficient \( C_{\varepsilon g} \), there is no consensus on a standard value.  Here, \( C_{\varepsilon g} = 1.44 \) is set as the default.  

If the turbulent heat flux is modeled using the gradient diffusion hypothesis, which is the most commonly used approach, it is expressed as:

$\overline{\theta u_i} = - \dfrac{\nu_{t}}{\mathit{Pr}_{t}} \dfrac{\partial T}{\partial x_i},$

where $\mathit{Pr}_t \approx 1$ is often used for passive scalar problems by the Reynolds analogy.

Instead, this code employs an algebraic turbulent heat flux model, given by:

$\overline{\theta u_i} = -C_{\theta}\,\dfrac{k}{\varepsilon+\varepsilon_0}
	\left( 
	\overline{u_{i} u_{j}} \dfrac{\partial T}{\partial x_j}
	+ \xi \overline{\theta u_{j}} \dfrac{\partial U_i}{\partial x_j}
	+ \eta \beta g_i \overline{\theta^{2}}
	\right).$

Here, $\overline{\theta^2}$ is the temperature variance (i.e., the Reynolds average of the square of the temperature fluctuation), obtained from its own transport equation:

$\dfrac{D \overline{\theta^{2}}}{D t}
	= - 2 \overline{\theta u_{i}} \dfrac{\partial T}{\partial x_i}
	- \dfrac{1}{R} \dfrac{\overline{\theta^2}}{k}\left( \varepsilon + \varepsilon_0\right) 
	+ \dfrac{\partial }{\partial x_i} 
	\left( \left(\alpha + \dfrac{\nu_t}{\mathit{Pr}_t}\right) \dfrac{\partial \overline{\theta^{2}}}{\partial x_i}  \right).$

 The model constants are $C_{\theta} = 0.15$, $\xi=0.6$, $\eta=0.6$, $R=0.75$, $\mathit{Pr}_t=1.0$.

 
At no-slip walls, $k=\varepsilon=\overline{\theta^2}=0$ should be imposed as Dirichlet boundary conditions.


# Notes
This code is written based on OpenFOAM's solver, "buoyantBoussinesqSimpleFoam".

All model constants and material properties (except for gravity) are defined in  "constant/transportProperties".

For instance, if you want to change the model constant to $C_{\varepsilon g}=1$, add the following line to "constant/transportProperties":
> Cepsg 1;

The variables, material properties, and model constants are named in the code as follows:

$k$ (k);
$\varepsilon$ (epsilon);
$\overline{\theta^2}$ (Tvar);
$\overline{\theta u_i}$ (thf);
$\nu$ (nu); 
$\alpha$ (alpha);
$\beta$ (beta);
$Q$ (Q);
$C_\mu$ (Cmu);
$C_{\varepsilon 1}$ (Ceps1);
$C_{\varepsilon 2}$ (Ceps2);
$C_{\varepsilon g}$ (Cepsg);
$\sigma_k$ (sigmaK);
$\sigma_\varepsilon$ (sigmaEps);
$C_\theta$ (Ct);
$\xi$ (Cxi);
$\eta$ (Ceta);
$R$ (CR);
$\mathit{Pr}_t$ (Prt)


# Numerical methods for stability
Details can be found at [https://arxiv.org/abs/2501.12973]

<footer>

<!--
  <<< Author notes: Footer >>>
  Add a link to get support, GitHub status page, code of conduct, license link.
-->

---

Get help: [Post in our discussion board](https://github.com/orgs/skills/discussions/categories/introduction-to-github) &bull; [Review the GitHub status page](https://www.githubstatus.com/)

&copy; 2024 GitHub &bull; [Code of Conduct](https://www.contributor-covenant.org/version/2/1/code_of_conduct/code_of_conduct.md) &bull; [MIT License](https://gh.io/mit)

</footer>
