<header>

<!--
  <<< Author notes: Course header >>>
  Include a 1280×640 image, course title in sentence case, and a concise description in emphasis.
  In your repository settings: enable template repository, add your 1280×640 social image, auto delete head branches.
  Add your open source license, GitHub uses MIT license.
-->

The current code is built on OpenFOAM-v2406. 
Below is a brief guide on how to use it.

# 1. Install OpenFOAM 
Follow the guidelines provided at [https://www.openfoam.com/] to install OpenFOAM. 

Note: This code is not compatible with the OpenFOAM Foundation release (often labeled as v8, ..., v12) available at [openfoam.org]

# 2. Set Up Your User Directory
After installation, update your personal "~/.bashrc" file to include OpenFOAM environment variables refer to the installation guidelines. 

Then, create the necessary directories for running OpenFOAM by executing:

  > mkdir -p $WM_PROJECT_USER_DIR/{run,applications,src}

Next, overwrite the provided source codes by copying the provided files into the appropriate location.

# 3. Compile the Solver
In a Linux environment, compile the solver located in "applications/mybuoyantBoussinesqSimpleFoam" by running the following commands sequentially:

  > app

  > cd mybuoyantBoussinesqSimpleFoam
  
  > wclean
  
  > wmake

# 4. Run the Examples  
The "run" folder contains the problems validated in the paper. 

You can run these examples just like any conventional OpenFOAM tutorial. 

Simply navigate to the example directory and execute:
  
  > run
  
  > cd RayleighBenardConvection
  
  > sh ./Allclean
  
  > sh ./Allrun


<footer>

<!--
  <<< Author notes: Footer >>>
  Add a link to get support, GitHub status page, code of conduct, license link.
-->

---

Get help: [Post in our discussion board](https://github.com/orgs/skills/discussions/categories/introduction-to-github) &bull; [Review the GitHub status page](https://www.githubstatus.com/)

&copy; 2024 GitHub &bull; [Code of Conduct](https://www.contributor-covenant.org/version/2/1/code_of_conduct/code_of_conduct.md) &bull; [MIT License](https://gh.io/mit)

</footer>
