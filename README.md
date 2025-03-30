# Age of Information in Energy Harvesting Aided Massive Multiple Access Networks

This repository provides a MATLAB implementation of the methods described in the following paper:

> **Fang, Zhengru, et al.**  
> _Age of Information in Energy Harvesting Aided Massive Multiple Access Networks_.  
> *IEEE Journal on Selected Areas in Communications*, vol. 40, no. 5, pp. 1441–1456, 2022.

The code focuses on simulating and optimizing the Age of Information (AoI) in massive multiple access networks that incorporate energy harvesting. It implements optimization algorithms using the convex-concave procedure (CCP) combined with MATLAB’s `fmincon` solver, and considers different multiple access schemes including TDMA, FDMA, and NOMA.

---

## File Structure

The repository is organized as follows:

```
project_root/
├── main.m                              % Main entry script to run the simulation and plotting routines
├── CreateModel.m                       % Constructs the system model (generates model parameters and constants)
├── Channel_create.m                    % Channel generation function based on node distances, setting channel gains, etc.
├── objfun.m                            % Objective function used in the optimization (for CCP iterations)
├── Power_allocation.m                  % Function for power allocation in the NOMA scheme
├── confun.m                            % Nonlinear constraint function (for general optimization problems)
├── confun_TDMA_convex.m                % Convexified constraint function for the TDMA scheme
├── confun_NOMA.m                       % Nonlinear constraint function for the NOMA scheme
├── confun_NOMA_convex.m                % Convexified constraint function for the NOMA scheme
├── Main_TDMA_MV_Num.m                  % Simulation script for the TDMA multiple access scheme (simulations over different numbers of MTCDs)
├── Main_FDMA_MV_Num.m                  % Simulation script for the FDMA multiple access scheme
├── Main_NOMA_MV_Num.m                  % Simulation script for the NOMA multiple access scheme
├── Main_NOMA_MV.m                      % Alternate version of the NOMA simulation script
```

---

## Overview

The code implements simulation and optimization for energy harvesting aided massive multiple access networks. The key features include:

- **System Modeling:**  
  The `CreateModel.m` script generates the system model with simulation parameters, including node numbers, channel parameters, power constraints, and energy harvesting settings.

- **Channel Generation:**  
  The `Channel_create.m` function sets up the channel gains based on node distances using a standard path-loss model (inverse square law with scaling).

- **Optimization Algorithms:**  
  The optimization problem (minimizing AoI under various constraints) is solved using the convex-concave procedure (CCP). The objective function is defined in `objfun.m` and the problem is constrained by several nonlinear functions (`confun.m`, `confun_TDMA_convex.m`, `confun_NOMA.m`, and `confun_NOMA_convex.m`).

- **Multiple Access Schemes:**  
  Separate simulation scripts are provided for TDMA (`Main_TDMA_MV_Num.m`), FDMA (`Main_FDMA_MV_Num.m`), and NOMA (`Main_NOMA_MV_Num.m` and `Main_NOMA_MV.m`). For NOMA, power allocation is calculated using the `Power_allocation.m` function.

- **Data Extraction and Plotting:**  
  Helper functions (`Get_data_from_struct_CCP_nonVaca.m` and `Get_data_from_struct_CCP_nonVaca_NOMA.m`) are used to extract simulation results for plotting. The main script (`main.m`) orchestrates the simulation runs and calls plotting routines to visualize results such as AoI and energy consumption.

---

## How to Run

1. **Requirements:**
   - MATLAB (R2018b or later is recommended)
   - MATLAB Optimization Toolbox

2. **Setup:**
   - Clone or download the repository.
   - Add the project directory and its subdirectories (if any) to your MATLAB path.
   - Ensure that any required data folders (e.g., `Data`, `figures`) are created in the project directory.

3. **Execution:**
   - Run the `main.m` script from within MATLAB. This script will sequentially call the simulation routines (for example, `Main_TDMA_MV_Num.m`, `Main_FDMA_MV_Num.m`, and `Main_NOMA_MV_Num.m`) and subsequently generate plots from the simulation results.
   - The simulation results are automatically saved (typically in a `Data/` folder) and plotted for analysis.

---

## Citation

If you use this code in your research, please cite the following paper:

```
@article{fang2022age,
  title={Age of Information in Energy Harvesting Aided Massive Multiple Access Networks},
  author={Fang, Zhengru and Wang, Jingjing and Ren, Yong and Han, Zhu and Poor, H Vincent and Hanzo, Lajos},
  journal={IEEE Journal on Selected Areas in Communications},
  volume={40},
  number={5},
  pages={1441--1456},
  year={May 2022},
  publisher={IEEE}
}
```
