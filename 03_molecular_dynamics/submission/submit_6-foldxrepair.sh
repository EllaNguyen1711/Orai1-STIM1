#!/bin/sh

#########################################
# CONFIG
#########################################

FOLDX="/home/exouser/Installers/foldx/foldx"
ROTABASE="/home/exouser/Installers/foldx/rotabase.txt"

#Give PATH to PDB storing directory

PDB_DIR=""
MUTFILE="$PDB_DIR/individual_list.txt"

RUNS=5

#########################################

for pdb in "$PDB_DIR"/*.pdb
do
    base=$(basename "$pdb" .pdb)

    # Skip already repaired files (POSIX shell compatible)
    case "$base" in
        *_Repair) continue ;;
    esac

    session="foldx_${base}"

    echo "Launching $session"

    tmux new-session -d -s "$session" "
        cd '$PDB_DIR' &&

        cp -n '$ROTABASE' . &&

        if [ -f '${base}_Repair.pdb' ]; then
            echo 'Repair exists -> skipping RepairPDB'
            repaired='${base}_Repair.pdb'
        else
            echo 'Running RepairPDB for ${base}'
            '$FOLDX' --command=RepairPDB --pdb='${base}.pdb'
            repaired='${base}_Repair.pdb'
        fi &&

        echo 'DONE ${base}'
    "
done

echo "All tmux jobs started."
echo "Check sessions with: tmux ls"

