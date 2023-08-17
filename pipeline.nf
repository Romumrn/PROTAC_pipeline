#!/usr/bin/env nextflow
 
params.ligase = "$baseDir/data/6hax/ligase.pdb"
params.target = "$baseDir/data/6hax/target.pdb"
params.refpdb = "$baseDir/data/6hax/ref.pdb"
params.restreint = "$baseDir/data/6hax/rest.dat"
params.outdir = "$baseDir/result"
params.glowworms = 100
params.simulationstep = 20
nextflow.enable.dsl=2


log.info """\
            Protac PIPELINE 

    https://github.com/Romumrn/PROTAC_pipeline
=================================================
ligase = $params.ligase
target = $params.target
reference = $params.refpdb
restreint = $params.restreint
outdir = $params.outdir

-------------------------------------------------
glowworms number = $params.glowworms
Simulation steps = $params.simulationstep
"""


// ./nextflow run pipeline.nf -with-docker pipeline_gilberto:v0

process Setup_calculation {
    input : 
        path ligase
        path target
        path restreint
        val glowworms
        val step

    output:
        path "lightdock_*" , emit: lightdock_files
        path "swarm*/*" , emit: swarm_directories
        
    """
    lightdock3_setup.py $ligase $target -g $glowworms --noxt --now --noh -rst $restreint --verbose_parser 
    
    lightdock3.py setup.json $step -s fastdfire -c 20
    """
}

process Get_conformation_and_cluster {
    //publishDir params.outdir
    tag "${meta}"

    input:
        path ligase
        path target
        path lighdock_file
        val step
        tuple val(meta), path(filepath)

    output:
        path "*" 
        
    """
    lgd_generate_conformations.py $ligase $target gso_${step}.out 200   
    lgd_cluster_bsas.py  gso_${step}.out

    for file in *; 
    do  
        mv -v \${file} ${meta}__\${file}; 
    done
    """
}

process Make_swarm_dir {
    input:
        path test

    output:
        path "*", emit: swarm_dir

    """
    for file in \$(ls ); do
        if [[ -f "\$file" ]]; then
            # Extract the prefix before "__"
            prefix="\${file%%__*}"
            
            # Create the directory if it doesn't exist
            mkdir -p "\$prefix"
            
            # Move the file into the corresponding directory
            mv "\$file" "\$prefix/\${file/\${prefix}__/''}"
        fi
    done
    """
}

process Rank_LighDock{
    publishDir params.outdir

    input:
        path swarm_dir
        val step

    output:
        path "rank_by_scoring.list" , emit: rank_by_scoring

    """
    s=`ls -d ./swarm_* | wc -l`;

    swarms=\$((s-1));	     
    lgd_rank_swarm.py \$s $step;
    lgd_rank.py \$s $step;
    """
}


process Extract_best_200_docking {
    publishDir params.outdir
    
    input:
        path rank
        path swarm_dir

    output:
        path "TOP200_PDBs/*", emit: top200
        
 
    """
    #Here we just want to extract the swarm directory with the better score
    mkdir TOP200_PDBs
    touch TOP200.dat
    head -n 201 rank_by_scoring.list >> TOP200.dat 

    cat TOP200.dat| tail -n +2 >> TOP200_without_head.dat 

    echo "Swarm Glowworm Name LD_Score" > Master
    cat TOP200.dat| awk '{print \$1,\$2,\$16,\$18}' | tail -n +3 >> Master
             
    while read line ; do
        swarm=\$(echo \$line | awk '{print \$1}')
        id=\$(echo \$line | awk '{print \$2}')
        echo "cp swarm_\${swarm}/lightdock_\${id}.pdb TOP200_PDBs/Lightdock_swarm_\${swarm}_\${id}.pdb"
        cp swarm_\${swarm}/lightdock_\${id}.pdb TOP200_PDBs/Lightdock_swarm_\${swarm}_\${id}.pdb

    done < TOP200_without_head.dat
    """
}
	

    //export PATH="/data3/rmarin/for_gilberto/pipeline/DockQ/:$PATH"
process Dock_Q_scores {
    //publishDir params.outdir
    tag "${meta}"

    input:
        tuple val(meta), path(pdb)
        path ref
        
    output:
        path "DockQ.dat"
        
 
    """
    DockQ.py $pdb $ref -useCA -perm1 -perm2 -verbose > DockQ_.dat
    cat DockQ_.dat | tail -5 > DockQ_summary.dat
    
    #..extract data to summary file
    Fnat=\$(cat DockQ_summary.dat | awk '{print \$2}' | head -1 )
    int_RMSD=\$(cat DockQ_summary.dat | awk '{print \$2}' | head -3 | tail -1 )
    Ligand_RMSD=\$(cat DockQ_summary.dat | awk '{print \$2}' | head -4 | tail -1 )
    DockQ=\$(cat DockQ_summary.dat | awk '{print \$2}' | tail -1 )

    new_field=''
    if (( \$(echo "(\$Fnat >= 0.5 && \$int_RMSD <= 1.09)" | bc -l) )) || (( \$(echo "(\$DockQ >= 0.80 )" | bc -l) )); then
        new_field+="High"
    elif (( \$(echo "(\$Fnat >= 0.3 && \$Fnat <= 0.5 && \$int_RMSD <= 2.09)" | bc -l) )) || (( \$(echo "(\$Fnat >= 0.5 && \$int_RMSD > 1.09 && \$int_RMSD <= 2.09)" | bc -l) )) || (( \$(echo "(\$DockQ >= 0.49 && \$DockQ < 0.80 ) " | bc -l) )); then
        new_field+="Medium"
    elif (( \$(echo "(\$Fnat >= 0.1 && \$Fnat <= 0.3 && \$int_RMSD <= 4.09)" | bc -l) )) || (( \$(echo "(\$Fnat >= 0.3 && \$int_RMSD > 2.09 && \$int_RMSD <= 4.09)" | bc -l) )) || (( \$(echo "(\$DockQ >= 0.23 && \$DockQ < 0.49) " | bc -l) )); then
        new_field+="Acceptable"
    elif (( \$(echo "(\$Fnat < 0.1 || \$int_RMSD > 4.09)  " | bc -l) )) && (( \$(echo "(\$DockQ < 0.23 )" | bc -l) )); then
        new_field+="Incorrect"
    else
        new_field+="Undefined"
    fi

    echo  "Name Fnat int-RMSD Ligand-RMSD DockQ-Score Rank" > DockQ.dat
    echo -e "$meta \$Fnat \$int_RMSD \$Ligand_RMSD \$DockQ \$new_field" >> DockQ.dat
    """
}

process Voromqa_scores{
    //publishDir params.outdir
    tag "${meta}"

    input:
        tuple val(meta), path(pdb)
        
    output:
        path "VORO_data.dat"
        
    """
    #!/bin/bash
    
    voronota-voromqa -i $pdb --score-inter-chain | tail +1 > VORO_data.dat
    """
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete" 
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}. Please contact xxx.mail"
}

workflow {
    Setup_calculation(params.ligase , params.target , params.restreint, params.glowworms, params.simulationstep)
    
    Get_conformation_and_cluster( params.ligase , params.target ,Setup_calculation.out.lightdock_files,params.simulationstep, Setup_calculation.out.swarm_directories.flatten().map { [it.toString().split('/')[-2], it] }.groupTuple() )
    
    Make_swarm_dir( Get_conformation_and_cluster.out.collect())

    Rank_LighDock( Make_swarm_dir.out.swarm_dir, params.simulationstep)

    Extract_best_200_docking( Rank_LighDock.out.rank_by_scoring, Make_swarm_dir.out.swarm_dir )

    Dock_Q_scores( Extract_best_200_docking.out.top200.flatten().map { [it.toString().split('/')[-1], it] }, params.refpdb)
        .concat( Channel.of( "Name Fnat int-RMSD Ligand-RMSD DockQ-Score Rank") )
        .collectFile(name: 'resulat_Dock_Q_scores.txt', skip: 1)
        .subscribe { f -> 
			f.copyTo("${params.outdir}")
        }
    
    Voromqa_scores( Extract_best_200_docking.out.top200.flatten().map { [it.toString().split('/')[-1], it] })
        .concat( Channel.of( "Name Fnat int-RMSD Ligand-RMSD DockQ-Score Rank") )
        .collectFile(name: 'resulat_Voromqa_scores.txt', skip: 1)
        .subscribe { f -> 
			f.copyTo("${params.outdir}")
        }

}