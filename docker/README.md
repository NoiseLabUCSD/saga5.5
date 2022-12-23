# Docker container for SAGA geoacoustic inversion software

By [William Jenkins](https://github.com/NeptuneProjects)  
21 November 2022

SAGA and its required FORTRAN dependencies can be installed using Docker. This allows for portability and eliminates the need to install or alter libraries on your system directly.

## Dockerfile

The `Dockerfile` is used to build a container image. Its contents include:
- Line 1: Specify Redhat’s Linux distribution to match the remote server we normally work on.
- Line 3: Install the following packages: Fortran compiler `gcc-gfortran`, debugger `gdb`, `git` for version control, `make` to compile SAGA, and `vim` to edit files within the container environment.
- Line 4: Clone the SAGA repository into the container.
- Line 5: Specify the working directory for SAGA.
- Line 6: Specify the host architecture.
- Lines 7-13: Install SAGA using a setup script.

## setup.sh

The setup script is necessary to add the SAGA working and bin directories to the container’s system path.

## Building the Docker container
Here is a simple way to build the container image and run it:
```bash
docker build . -t IMAGENAME
docker run -v /HOSTDIR:/saga IMAGENAME
```
`IMAGENAME` is the name of the container image, and `HOSTDIR` is the [local] directory you want to mount to the Docker image as a volume.
