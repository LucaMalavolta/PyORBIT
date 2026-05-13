# ✏️ Start writing
 
 If you have a Windows computer, I suggest using the WindowsSubsystem for Linux 
 
 In all cases, I suggest using the Python environment by conda

https://docs.conda.io/projects/conda/en/stable/index.html

Instead Google Colab

## Windows installation of conda

https://visualstudio.microsoft.com/it/visual-cpp-build-tools/

The program will ask you to install additional packages,you will need the C++ development tools for Windows (it shopuld be the first box in the upper left area). The download will require more than 6 GB

## Setting up an environment in `conda`

Before proceeding with the installation, I suggest creating an environment dedicated to `PyORBIT` using **Python 3.10**. This is the version I'm currently use to test the code.

With conda/anaconda:

```{code} bash
conda create --name lucchin_school python=3.10
```

If you are using a Mac with ARM architecture, the following command should force conda to install x86_64 versions of Python and all packages (many thanks to Jinglin Zhao for the tip):

```{code} bash
CONDA_SUBDIR=osx-64 conda create -n lucchin_school python=3.10
```

To list the available environments, do:

```{code} bash
conda env list
```

The active environment will be marked with a \*

To activate the `lucchin_school` environment:

```{code} bash
WINDOWS: activate lucchin_school
LINUX, macOS: conda activate lucchin_school
```


## Google Colab

You can use a free Google Colab account to run the code 

Inside a Google Colab cell, execute the following commands
```
!pip install pytransit
```


Write here. Publish when ready. Share the link.

Your draft stays private until you choose to publish.