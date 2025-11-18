# Reduced-Adjoint Variational Data Assimilation (RA-VDA)

This repository contains the MATLAB scripts used to implement the Reduced-Adjoint
Variational Data Assimilation (RA-VDA) framework for coupled estimation of
root-zone soil moisture and soil hydraulic parameters using HYDRUS-1D.

The workflows include:
- Generation of Proper Orthogonal Decomposition (POD) bases  
- Construction of the reduced-order forward model  
- Reduced adjoint computation  
- Optimization loops for coupled state–parameter estimation  
- Running HYDRUS-1D and reading model outputs  

---

## Repository Contents

### **1. POD Generation**
- **`P_M_computing_param1.m`**  
  Generates PODs and reduced-order model when **initial soil moisture (S0) is known**  
  and **parameters are the only unknowns**.

- **`P_M_computing_S0_param1.m`**  
  Generates PODs and reduced-order model when **S0 and parameters are jointly estimated**.

---

### **2. Reduced-Adjoint Optimization Scripts**
- **`Reduced_adjoint_script_S0_params.m`**  
  Main script for **coupled estimation** (S0 and parameters).  
  Implements the two-step alternating minimization loop:
  1. **`Reduced_Adjoint_fun_S0_params.m`** – computes J & ∇J for S0  
  2. **`Reduced_Adjoint_fun_S0_params2.m`** – computes J & ∇J for parameters  
  Repeats until convergence.

- **`Reduced_adjoint_script_param.m`**  
  Optimization when **S0 is known** and only parameters are estimated.  
  Uses **`Reduced_Adjoint_fun_param.m`** for cost function and gradient.

---

### **3. HYDRUS-1D Forward Model Interface**
- **`Hydrus_Run_param.m`**  
  Prepares input files and runs HYDRUS-1D.

- **`reading_flx_Hy.m`**  
  Reads fluxes and variables from HYDRUS-generated ASCII output files.

---

## Requirements

- MATLAB (R2020b or later recommended)  
- HYDRUS-1D (free software; see https://www.pc-progress.com/en/Default.aspx?hydrus-1d)

---

## How to Run

1. Configure HYDRUS-1D project files.
2. Generate POD bases using the appropriate `P_M_*` script.
3. Run either:
   - `Reduced_adjoint_script_param.m`, or  
   - `Reduced_adjoint_script_S0_params.m`
4. View results stored in the `/results/` folder (user-created).

---

## Data Availability

Synthetic datasets used in the associated manuscript are generated directly from
HYDRUS-1D. Example inputs and instructions are provided in this repository.  
Additional configuration files used in experiments can be shared upon request.

---

## License

This repository is released under the **MIT License** to support transparency
and reproducibility.  
See the `LICENSE` file for details.

---

## Citation

If you use this code, please cite:

Heidary, P., Farhadi, L., & Altaf, M.U. (2025).  
*Coupled Estimation of Root zone Soil Moisture and Soil Hydraulic Parameters with Reduced-Adjoint Variational Data Assimilation using Near-Surface Soil Moisture Observations*  
Environmental Modelling & Software.
