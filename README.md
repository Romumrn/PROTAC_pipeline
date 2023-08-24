# PROTAC_pipeline

Provide a nextflow pipeline for "Rational prediction of PROTAC compatible ligase-target interfaces" (https://github.com/GilbertoPPereira/PROTACability). 

## Prerequisites:

You need the following tools and scripts installed:
- Nextflow (https://www.nextflow.io/docs/latest/getstarted.html)
- Docker (https://docs.docker.com/engine/install/)

## Analyse
### Automatic Restraints (Step 1):
Use MDAnalysis to identify residues near the ligand.
Employ NACCESS to compute the change in surface area (dASA) for all residues.
Select residues with dASA > 25% and within 6Å of the ligand.

### Docking Simulation (Step 2):
Employ lightdock3_setup.py to prepare docking simulation.
Utilize ligand-free receptor and target proteins.
Specify parameters like the number of dockings and restraint file.
Perform the docking simulation.

### Conformation Generation (Step 3):
Utilize lgd_generate_conformations.py to generate potential solutions for each swarm.

### Clustering (Step 4):
Apply lgd_cluster_bsas.py to minimize redundancy in binding pose predictions within each swarm.

### Ranking (Step 5):
Use lgd_rank.py to rank binding poses based on LightDock scores.

### DockQ Evaluation (Step 6):
Utilize DockQ to assess predicted poses against reference structures.
Rank the poses using DockQ scores.

### Energy Rescoring (Step 7):
Use voronota-voromqa to perform energy-based rescoring of interface conformations.

### Filtering and Final Ranking (Step 8):
Transfer ligands from original structures to predicted poses.
Utilize jwalk to calculate Surface Accessible Surface Distance (SASD) for ranking.
Generate a comprehensive ranking and assess the accuracy of the predictions.

## Process Description: 
- Setup_calculation: Get the automatic restraints...,
- Get_conformation_and_cluster: Generate conformations and perform clustering...,
- Make_swarm_dir: Organize results into swarm directories...,
- Rank_LightDock: Rank poses using LightDock scores...,
- Extract_best_200_docking: Extract top 200 docking results...,
- Dock_Q_scores: Evaluate docking results using DockQ...,
- Voromqa_scores: Perform energy-rescoring using VoroMQA...,
- Align_Pymol: Align ligands and calculate SASD...,
- Run_Jwalk: Calculate SASD using Jwalk...

## Add input explanation 
