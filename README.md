# Code for RV mechanics in early-stage PAH: Right ventricular mechanics in early-stage pulmonary arterial hypertension: a computational study

This repository contains supporting code for th 

In this code we try in estimate passive and active properties for the RV assuming known passive and active properties for the LV.


## Install

### Using docker
You can run the code using the provided docker image. For this you need to have [Docker](https://docs.docker.com/get-docker/). Once docker is installed you can pull the image
```
docker pull ghcr.io/oscarodeigah/rv_pah_project:latest
```
and start a new container (sharing your current directory with the container)
```
docker run --name rv_pah_project -w /home/shared -v $PWD:/home/shared -it ghcr.io/oscarodeigah/rv_pah_project:latest
```
The scripts are located in `/app` inside the container

### Using pip
The code in this repository is pure python and can therefore be installed with `pip` using the command
```
PIP_NO_BINARY="h5py" python3 -m pip install git+https://github.com/oscarodeigah/rv_pah_project
```
However, to run the code you need to have FEniCS (see https://fenicsproject.org/download/archive/ for more info about how to install FEniCS).


## Getting started

There is an example dataset in the `data/sample_datafiles` folder which contains the following files
- `CNT.h5` - file containing the mesh, facet markers and fiber orientations
- `CNT.yml` - file containing pressure-volume data for the LV and RV
- `LV_active_data.json` - file containing information about the active properties of the LV throughout the cycle

There are two different scripts located in the [`scripts`](scripts) folder

To run a simple optimization to estimate the RV passive and active properties you can run the script [`run_optimization.py`](scripts/run_optimization.py).

For the sample datasets this takes about 2.5 hours on a regular laptop. 

Once this is done you can run the [run_postprocessing.py](scripts/run_postprocessing.py) which will extract different metrics from the simulation

## Reproducing figures in the paper

We have also collected all the results in the paper in the folder called ["postprocessed_results"](postprocessed_results). To recreate the figures, you first need to install some dependencies

```
python3 -m pip install matplotlib numpy PyYAM
```
and then you can run 
```
python3 make_paper_figures.py
```
inside the ["postprocessed_results"](postprocessed_results) folder.


## Issues
If you face any issues or have any question, feel free to open an issue in the [issue tracker](https://github.com/oscarodeigah/rv_pah_project/issues/new) or send an email to <oscaro@simula.no>