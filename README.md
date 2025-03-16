# LearnFEM

Learning the theory of the finite element method can be a daunting task for the practicing engineer or (under)graduate student.
The finite element method pulls together many different mathematical concepts and tools into a "super-algorithm" and in my experience it is common for students of the method to become lost in the details and losing sight of the bigger picture. 
Many of these fine details are rooted in an effort to demonstrate an efficient implementation of the finite element method and, not insignificantly, due to a historical necessity to work with floating-point arithmetic.

And so while there are numerous textbooks that dive immediately into the details of implementation, in this lecture set we will begin by neglecting the concerns of efficiency and pragmatism of floating-point arithmetic for as long as possible -- only turning towards them once we have demonstrated the "big picture" of the finite element method and then desire to improve computational performance.
Specifically, we will utilize symbolic calculations via Matlab's Symbolic Toolbox, which will allow us to initially avoid the topic (and complexity) of quadrature and instead use *symbolic* integration and effortlessly compute derivatives (e.g., of basis functions) using symbolic operations.
Instead of immediately introducing piecewise polynomials, we will begin with a single element and use degree refinement to increase accuracy (at the expense of computational cost).
And finally, rather than using mathematical manipulations common to the derivation of weak formulations used in the finite element method, we will instead present FEM from a change-of-basis point-of-view.
Only after demonstrating Galerkin's method on a variety of partial differential equations will we gradually apply these "best practices" into our method to improve computational efficiency.

## Lectures
  * Lecture 0: Mathematical Preliminaries
  * Lecture 1: Introducing Galerkin's Method via Scalar Projection
  * Lecture 2: Adding Boundary Conditions to Galerkin's Method
  * Lecture 3: Solving Partial Differential Equations with Galerkin's Method
  * Lecture 4: Piecewise Polynomials in Galerkin's Method
  * Lecture 5: Improving Performance through Improved Weak Formulations
  * Lecture 6: Improving Performance through Quadrature
  * Lecture 7: Improving Performance through Floating-Point Arithmetic

## Building the Documentation

### Install Quarto
[Quarto](https://quarto.org/docs/get-started/)

### Install a TeX distribution
* [TeXLive](https://www.tug.org/texlive/windows.html#install)
* or TinyTex: `quarto install tinytex`

### Python Instructions
#### Create a Virtual Environment
```
python3.12.exe -m venv env
python3.12.exe -m pip install -r mkernel_python_requirements.txt
```
##### Activate Python Environment (Powershell)
```
Set-ExecutionPolicy -ExecutionPolicy Unrestricted -Scope CurrentUser
env\Scripts\Activate.ps1
```

### MATLAB Instructions
#### Create a Virtual Environment
The Matlab Engine for version 24.1.2 requires 3.9<= ver <3.12.
```
python3.11.exe -m venv mkernel_python_env
python3.12.exe -m pip install -r requirements.txt
```

#### Activate Python Environment (Powershell)
```
Set-ExecutionPolicy -ExecutionPolicy Unrestricted -Scope CurrentUser
mkernel_python_env\Scripts\Activate.ps1
```

#### Install MKernel
See: [MKernel: A Jupyter Kernel for Matlab](https://github.com/allefeld/mkernel)
```
python.exe -m pip install jupyter
python.exe -m pip install git+https://github.com/allefeld/mkernel.git
```

#### Freeze Python Environment
```
python.exe -m pip freeze > mkernel_python_requirements.txt
```

#### Install Python Environment
```
python.exe -m pip install -r mkernel_python_requirements.txt
```