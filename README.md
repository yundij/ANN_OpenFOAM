
# ANN-Drag-Model

This project provides two examples on how to deploy ANN models trained using Keras in C++ based OpenFOAM and CFDEM 

The two examples represent two approaches to simulate gas-particle flows:
- Euler-Euler approach
- Euler-Languange approach

Euler-Euler approach treates two phases as interpenetrating continua, and we use `twoPhaseEulerFOAM` solver in `OpenFOAM-3.x` for this simulation. 
Euler-Language approach solves the continuous phase on eulerian grids and tracks each particle and solve their behavior using Newtons equations of motion, and we use `CFDEM` for this simulation. 

Large scale simulations requires coarse mesh grid, which needs sub-grid drag models. Detailed discussions can be found in https://www.sciencedirect.com/science/article/pii/S0032591018310192
We developed sub-grid drag models using Keras, and then deployed them into the simulation platforms. 
- WenYuDriftFlux is located in `/OpenFOAM-3.0.1/applications/solvers/multiphase/twoPhaseEulerFoamPU/interfacialModels/dragModels/`
- WenYuANNDrag is located in `CFDEM/CFDEMcoupling-PUBLIC-5.x/src/lagrangian/cfdemParticle/subModels/forceModel/`

The deployment is based on `keras2cpp` (`https://github.com/pplonski/keras2cpp`). 

