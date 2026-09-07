#!/bin/sh

FOLDX="/home/exouser/Installers/foldx/foldx"
ROTABASE="/home/exouser/Installers/foldx/rotabase.txt"

PDB_DIR=""
MUTFILE="$PDB_DIR/individual_list.txt"

RUNS=5

#########################################

for pdb in "$PDB_DIR"/*.pdb
do
    base=$(basename "$pdb" .pdb)

    session="foldx_${base}"
    logfile="${PDB_DIR}/${base}.foldx.log"

    echo "Launching $session  →  $logfile"

    tmux new-session -d -s "$session" "
        cd '$PDB_DIR' &&

        {
            echo '===== FoldX job: ${base} ====='
            date

            # Always keep FoldX-required files local
            cp -n '$ROTABASE' ./rotabase.txt
            cp '$MUTFILE' ./individual_list.txt

            repaired='${base}.pdb'

            if [ ! -f \"\$repaired\" ]; then
                echo 'ERROR: repaired PDB not found: '\$repaired
                exit 1
            fi

            echo 'Running BuildModel'
            '$FOLDX' --command=BuildModel \
                     --pdb=\$repaired \
                     --mutant-file=individual_list.txt \
                     --numberOfRuns=$RUNS

            echo 'DONE ${base}'
            date
        } > '$logfile' 2>&1
    "
done

echo "All tmux jobs started."
echo "Logs saved as: *.foldx.log"
echo "Check sessions with: tmux ls"

