# Influenza-with-Vaccination-and-Diffusion PDE model
## Members:
- [Mariem Sherif](https://github.com/orgs/Influenza-PDE-solutions/people/mariamsherif04)
- [Menna Ashraf](https://github.com/orgs/Influenza-PDE-solutions/people/Menaashraf)
- [Nayera Sherif](https://github.com/Nayera5)
- [Samar Hatem](https://github.com/samar04052004)
- [Sarah Sameh](https://github.com/orgs/Influenza-PDE-solutions/people/sarah012-210)
- [Shahd Ayman](https://github.com/Shahd-Ayman5)
- [Nada Hassan](https://github.com/Nadahassan147)
- [Nada Hesham](https://github.com/Nada-Hesham249)

## Github Webpage [here](https://influenza-pde-solutions.github.io/Influenza-with-Vaccination-and-Diffusion/).

## Problem Statement & Solution Overview:
The original problem models the spread of influenza with vaccination and spatial diffusion, and was solved in the reference work using `the Method of Lines (MOL)`.
, as required to implement two numerical solution schemes in addition to the book,
### we implemented 4 numerical solutions:

- Explicit finite difference scheme

- Crank–Nicolson method

- Finite volume method

- Strang splitting
---
After benchmarking the above classical numerical approaches,
### we developed a `Physics-Informed Neural Network (PINN)` model,
to solve the same system in a deep learning framework, and compared its performance and accuracy against the traditional methods and the reference solution.

PINNs integrate the physical laws of the system—represented by differential equations—directly into the loss function of a neural network.

Inputs:
- The network takes spatial and temporal inputs (x, t).
- 6 hidden layers with each 256 nueron
- labled data from the reference solution (20 samples)
  
Outputs:
- The model outputs the predicted values for each of the five variables: S(x,t), V(x,t), E(x,t), I(x,t), R(x,t).

Loss Function:
Our total loss is a weighted sum of:
- PDE residuals (calculated via automatic differentiation)
- Initial condition (IC) errors
- Boundary condition (BC) errors.

Training Strategy:
-  used an adaptive residual-based refinement strategy to focus on high-error regions.
-  used mean squared error (MSE) for all components of the loss.
---
# Visualization of Results
To assess the model behavior and accuracy, we visualized the predicted dynamics of the compartments (S, V, E, I, R) in two different ways:

## 1. Spatiotemporal Profiles
For selected time points, we plotted the solution across the spatial domain 𝑥.
- This shows how each compartment evolves across space at different times.

### A) Numerical 
![S](https://github.com/user-attachments/assets/b5dc99fc-a3a4-49c0-8d86-5e126849adc3)

### B) Deep learning 
![Screenshot 2025-06-19 020751 - Copy (2)](https://github.com/user-attachments/assets/ce0869a4-5e8f-4165-908d-69a36d541514)


## 2. Temporal Evolution at a Fixed Location
We fixed 𝑥=0 and plotted the solution over time.
- This helps us observe how each variable changes over time at a specific spatial point.

### A) Numerical 
![S with t](https://github.com/user-attachments/assets/7e1966ff-0bda-4890-aea9-cd926a55b7b5)

### B) Deep learning
![Screenshot 2025-06-19 020851](https://github.com/user-attachments/assets/47a82b69-cf14-417d-bcb4-534997287cc0)




