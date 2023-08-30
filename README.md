# PROTAC_pipeline

A Nextflow script for running the PROTAC pipeline from "Rational prediction of PROTAC compatible ligase-target interfaces" (https://github.com/GilbertoPPereira/PROTACability). This pipeline combines restraint-based LightDock,
energy-based rescoring, and minimal solvent-accessible surface distance filtering to produce PROTAC-compatible Protein-Protein Interactions (PPIs).

The Nextflow pipeline offers a seamless solution for executing complex scientific workflows. When combined with Singularity, the power of this pipeline is further amplified. Singularity's unique ability to directly pull Docker images from DockerHub ensures instant access to a diverse array of tools and resources. This amalgamation streamlines workflow management, enhances reproducibility, and fosters collaborative research by effortlessly bridging the gap between environment standardization and task orchestration. The Nextflow-Singularity integration facilitates efficient and secure execution, advancing scientific endeavors through simplified workflow deployment and robust containerization.
It's essential to have Singularity installed on your computer. If not, please follow this link for installation instructions: [Installing Singularity.]( https://docs.sylabs.io/guides/3.0/user-guide/installation.html) 

  Usage: 

```
  ./nextflow run pipeline.nf
```

Please modify the file `nextflow.config` to fit with your files, custom the number  and parametrize the number of CPUs, etc...

___                  

## Prerequisites:

You need the following tools and scripts installed:
- Nextflow (https://www.nextflow.io/docs/latest/getstarted.html)
- Singularity (https://docs.sylabs.io/guides/3.0/user-guide/installation.html)


## Analyze
### Automatic Restraints (Step 1): Get the automatic restraints
- Use `MDAnalysis` to identify residues near the ligand.
- Employ NACCESS to compute the change in surface area (dASA) for all residues.
- Select residues with dASA > 25% and within 6Å of the ligand.

### Docking Simulation (Step 2): Generate conformations and perform clustering
- Employ `lightdock3_setup.py` to prepare docking simulation.
- Utilize ligand-free receptor and target proteins.
- Specify parameters like the number of dockings and restraint file.
- Perform the docking simulation.

### Conformation Generation (Step 3):
- Utilize `lgd_generate_conformations.py` to generate potential solutions for each swarm.

### Clustering (Step 4):
- Apply `lgd_cluster_bsas.py` to minimize redundancy in binding pose predictions within each swarm.

### Ranking (Step 5):
- Use `lgd_rank.py` to rank binding poses based on LightDock scores.

### DockQ Evaluation (Step 6):
- Utilize `DockQ` to assess predicted poses against reference structures.
- Rank the poses using `DockQ` scores.

### Energy Rescoring (Step 7):
- Use `voronota-voromqa` to perform energy-based rescoring of interface conformations.

### Filtering and Final Ranking (Step 8):
- Transfer ligands from original structures to predicted poses.
- Utilize `jwalk` to calculate Surface Accessible Surface Distance (SASD) for ranking.
- Generate a comprehensive ranking and assess the accuracy of the predictions. 


## Output file :
The output is a csv file containing the final ranking of the ligands. The columns are as follows:
...

You will also have the best 200 pdb files.