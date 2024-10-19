## Install Quarto
[Quarto](https://quarto.org/docs/get-started/)

## Install a TeX distribution
* [TeXLive](https://www.tug.org/texlive/windows.html#install)
* or TinyTex: `quarto install tinytex`

## Python Instructions
### Create a Virtual Environment
```
python3.12.exe -m venv env
python3.12.exe -m pip install -r mkernel_python_requirements.txt
```
#### Activate Python Environment (Powershell)
```
Set-ExecutionPolicy -ExecutionPolicy Unrestricted -Scope CurrentUser
env\Scripts\Activate.ps1
```

## MATLAB Instructions
### Create a Virtual Environment
The Matlab Engine for version 24.1.2 requires 3.9<= ver <3.12.
```
python3.11.exe -m venv mkernel_python_env
python3.12.exe -m pip install -r requirements.txt
```

### Activate Python Environment (Powershell)
```
Set-ExecutionPolicy -ExecutionPolicy Unrestricted -Scope CurrentUser
mkernel_python_env\Scripts\Activate.ps1
```

### Install MKernel
See: [MKernel: A Jupyter Kernel for Matlab](https://github.com/allefeld/mkernel)
```
python.exe -m pip install jupyter
python.exe -m pip install git+https://github.com/allefeld/mkernel.git
```

### Freeze Python Environment
```
python.exe -m pip freeze > mkernel_python_requirements.txt
```

### Install Python Environment
```
python.exe -m pip install -r mkernel_python_requirements.txt
```