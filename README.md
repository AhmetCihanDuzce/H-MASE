# H-MASE: Hybrid Multi-Agent State Estimation for Stochastic Urban Networks

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This repository contains the official MATLAB simulation source code for the paper: **"Resilient Multi-Agent State Estimation for Smart City Traffic: A Systems Engineering Approach to Emission Mitigation"**.

## Overview

The **Hybrid Multi-Agent State Estimation (H-MASE)** protocol is a fully decentralized, fault-tolerant framework designed to monitor macroscopic traffic networks. Formulated as a distributed least-squares convex optimization problem, the protocol leverages a topology-induced cyber-physical graph comprising:
* **Physical Sensor Agents (PSAs):** Process local stochastic flow measurements.
* **Virtual Logic Agents (VLAs):** Act as fault-immune topological relays, preserving network-wide algebraic connectivity without contracting the state space.

This simulation validates the theoretical Input-to-State Stability (ISS) bounds, spatial variance reduction, and the autonomous fault-isolation mechanisms mathematically proven in the paper.

## Features of the Simulation

1.  **Synthetic Urban Topology Generation:** Automatically constructs a realistic, topologically heterogeneous ring-radial network (15 internal junctions, 10 boundary nodes, 55 directed edges).
2.  **Strict Operations Research Framework:** Simulates traffic evolution as bounded random walks under persistent stochastic measurement noise, without relying on exact state-transition matrices.
3.  **Distributed Projected Gradient Descent:** Implements the inner-loop consensus mechanism where agents reach global observability strictly via local hop-by-hop communications.
4.  **Autonomous Fault Isolation:** Injects a catastrophic structural sensor fault ($f = 80.0$) and demonstrates the real-time, autonomous constraint-dropping logic ($\sigma_i \to 0$) operating with a proven zero false-alarm rate.

## ⚙️ Requirements & Installation

* **Software:** MATLAB.
* **Toolboxes:** No specialized toolboxes are strictly required, though MATLAB's built-in Graph Theory functions (`digraph`, `shortestpath`, `plot`) are utilized for topology generation and visualization.

**To run the simulation:**
1. Clone this repository or download the `HMASE.m` file.
2. Open MATLAB and navigate to the folder containing the file.
3. Run the script from the command window:
   ```matlab
   run HMASE.m
