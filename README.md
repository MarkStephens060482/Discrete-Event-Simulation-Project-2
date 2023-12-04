# Factory X Production Line Simulation

## Overview

This repository contains simulation code and analysis conducted for optimizing the production line at Factory X, a manufacturing company specializing in lamp production. The simulation evaluates the impact of system modifications and an ongoing advertising campaign on Factory X's production schedules.

## Scenario

At Factory X, the production line consists of two critical stages, referred to as Machine 1 and Machine 2. Incoming orders for lamps queue up until Machine 1 initiates manufacturing. Following Machine 1 processing, lamps queue again, awaiting Machine 2. Notably, a limited queue space for lamps waiting at Machine 2 may impact overall production.

## Questions Addressed

The simulation aims to answer key operational questions:

1. **Sufficiency of Waiting Space:** Assessing the adequacy of the queue space for lamps awaiting Machine 2 processing.
2. **Optimization Potential:** Evaluating the impact of enhancing the speed of Machine 1 or Machine 2 within financial constraints.
3. **Impact of Advertising Campaign:** Assessing the system's preparedness for a 25% increase in order arrivals due to an ongoing advertising campaign.

## Files Included

- `factory_simulation_2_update.jl`: Julia code implementing the simulation model.
- `factory_simulation_2_update_simharness.jl`: Simulation Harness for simulation, Executable file
- `factory_simulation_2_update_simharness.ipynb`: Run the Jupyter notebook file
- `factory_simulation_2_update_analysis.ipynb`: Jupyter Notebook containing detailed analysis, visualizations, and insights derived from simulation results.
- `data/`: Directory containing relevant datasets used in the simulation and analysis.

## Usage

### Prerequisites

- Julia 1.9.3
- Jupyter Notebook
###  Instructions
1. Download the files into your working directory.
2. Run the simulation harness notebook 
3. Run the Analysis Notebook
