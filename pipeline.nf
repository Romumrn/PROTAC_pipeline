#!/usr/bin/env nextflow
 
params.ligase = "$baseDir/data/6hax/ligase.pdb"
params.target = "$baseDir/data/6hax/target.pdb"
params.refpdb = "$baseDir/data/6hax/ref.pdb"
params.restreint = "$baseDir/data/6hax/rest.dat"
params.outdir = "$baseDir/result"

nextflow.enable.dsl=2


log.info """\
            PIPELINE 
================================
ligase = $params.ligase
target = $params.target
reference = $params.refpdb
restreint = $params.restreint
outdir = $params.outdir
"""


// ./nextflow run pipeline.nf -with-docker pipeline_gilberto:v0

process Setup_calculation {
    input : 
        path ligase
        path target
        path restreint

    output:
        path "lightdock_*" , emit: lightdock_files
        path "swarm*/*" , emit: swarm_directories
        
    """
    lightdock3_setup.py $ligase $target -g 20 --noxt --now --noh -rst $restreint --verbose_parser 
    
    lightdock3.py setup.json 10 -s fastdfire -c 20
    """
}

process Get_conformation_and_cluster {
    //publishDir params.outdir
    tag "${meta}"

    input:
        path ligase
        path target
        path lighdock_file
        tuple val(meta), path(filepath)

    output:
        path "*" 
        
    """
    lgd_generate_conformations.py $ligase $target gso_10.out 200   
    lgd_cluster_bsas.py  gso_10.out

    for file in *; do mv -v \${file} ${meta}__\${file}; done
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

process Rank_LD{
    publishDir params.outdir

    input:
        path swarm_dir

    output:
        path "rank_by_scoring.list" , emit: rank_by_scoring

    """
    s=`ls -d ./swarm_* | wc -l`;

    swarms=\$((s-1));	     
    lgd_rank_swarm.py \$s 10;
    lgd_rank.py \$s 10;
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
process Dock_Q {
    publishDir params.outdir

    input:
        path TOP200
        path ref
        
    output:
        path "DockQ.dat"
        
 
    """
    PDBS=\$(ls Lightdock*_swarm* | sed 's/.pdb//g')
    echo "Name Fnat int-RMSD Ligand-RMSD DockQ-Score" > DockQ.dat

    for kk in \$PDBS ; do
        echo \$kk
        DockQ.py \${kk}.pdb $ref -useCA -perm1 -perm2 -verbose > DockQ_\${kk}.dat
        cat DockQ_\${kk}.dat | tail -5 > DockQ_\${kk}_summary.dat
        
        #..extract data to summary file
        Fnat=\$(cat DockQ_\${kk}_summary.dat | awk '{print \$2}' | head -1 )
        iRMS=\$(cat DockQ_\${kk}_summary.dat | awk '{print \$2}' | head -3 | tail -1 )
        LRMSD=\$(cat DockQ_\${kk}_summary.dat | awk '{print \$2}' | head -4 | tail -1 )
        DockQ=\$(cat DockQ_\${kk}_summary.dat | awk '{print \$2}' | tail -1 )
        echo "\${kk:10} \${Fnat} \${iRMS} \${LRMSD} \${DockQ}" >> DockQ.dat
    done
    """
}

//cut into 3 ctagerories and sort, and generate different top file ?

process Process_scores{
    publishDir params.outdir
    input:
        path score
        
    output:
       path "processed_score.dat"
        
    """
    #!/usr/bin/python3
newfile = ''
with open( "$score", "r") as score_file:
    for raw_line in score_file.readlines():
        line = raw_line.strip("\\n")
        if line.startswith("Name"): 
            newfile += line + " Rank\\n"
        else:
            Name, Fnat, int_RMSD, Ligand_RMSD, DockQ_Score = line.split(" ") 

            if (float(Fnat) >= 0.5 and float(int_RMSD) <= 1.09) or float(DockQ_Score) >= 0.80:
                newfile += line + ' High\\n' 
            elif (float(Fnat) >= 0.3 and float(Fnat) <= 0.5 and float(int_RMSD) <= 2.09) or (float(Fnat) >= 0.5 and float(int_RMSD) > 1.09 and float(int_RMSD) <= 2.09) or (float(DockQ_Score) >= 0.49 and float(DockQ_Score) < 0.80):
                newfile += line + ' Medium\\n'
            elif (float(Fnat) >= 0.1 and float(Fnat) <= 0.3 and float(int_RMSD) <= 4.09) or (float(Fnat) >= 0.3 and float(int_RMSD) > 2.09 and float(int_RMSD) <= 4.09) or (float(DockQ_Score) >= 0.23 and float(DockQ_Score) < 0.49):
                newfile += line + ' Acceptable\\n'
            elif (float(Fnat) < 0.1 or float(int_RMSD) > 4.09) and float(DockQ_Score) < 0.23:
                newfile += line + ' Incorrect\\n'
            else:
                newfile += line + ' Undefined\\n'

with open( "processed_score.dat", "w") as out:
    out.write(newfile)
    """
}

// process Voromqa_scores{
//     publishDir params.outdir
//     input:
//         path TOP200
        
//     output:
//         path "VORO_data.dat"
        
//     """
//     #!/bin/bash
 
//     PDBS=\$(ls Lightdock*_swarm* | sed 's/.pdb//g')
//     echo "Name Fnat int-RMSD Ligand-RMSD DockQ-Score" > DockQ.dat

//     for kk in \$PDBS ; do
//         voronota-voromqa -i \${kk}.pdb --score-inter-chain > Voro_\${kk}.dat
//         cat Voro_\${kk}.dat >> VORO_data.dat
//     done
//     sed -i 's/_/-/g' VORO_data.dat
//     """
// }

// process Voromqa_scores{
//     //publishDir params.outdir
//     tag "${meta}"

//     input:
//         tuple val(meta), path(pdb)
        
//     output:
//         path "VORO_data.dat"
        
//     """
//     #!/bin/bash
//     #echo "Name Fnat int-RMSD Ligand-RMSD DockQ-Score" > DockQ.dat

//     voronota-voromqa -i $pdb --score-inter-chain > VORO_data.dat
//     """
// }

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete" 
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}. Please contact xxx.mail"
}

workflow {
    Setup_calculation(params.ligase , params.target , params.restreint)
    
    Get_conformation_and_cluster( params.ligase , params.target ,Setup_calculation.out.lightdock_files, Setup_calculation.out.swarm_directories.flatten().map { [it.toString().split('/')[-2], it] }.groupTuple() )
    
    Make_swarm_dir( Get_conformation_and_cluster.out.collect())

    Rank_LD( Make_swarm_dir.out.swarm_dir)

    Extract_best_200_docking( Rank_LD.out.rank_by_scoring, Make_swarm_dir.out.swarm_dir )

    Dock_Q( Extract_best_200_docking.out.top200, params.refpdb)
    
    Process_scores( Dock_Q.out)

    //Voromqa_scores( Extract_best_200_docking.out.top200.flatten().map { [it.toString().split('/')[-2], it] }.groupTuple())
        // .collectFile(name: 'sample.txt', newLine: true)
        // .subscribe { f -> 
		// 	f.copyTo("${params.result}")
        // }
    
}