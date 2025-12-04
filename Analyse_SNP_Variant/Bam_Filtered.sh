#!/bin/bash

echo "=== TRIER ET INDEXER LES FICHIERS BAM ==="
echo ""

# Vérifier que samtools est installé
if ! command -v samtools &> /dev/null; then
    echo "Erreur: samtools n'est pas installé"
    exit 1
fi

# ===== SÉLECTION DES FICHIERS BAM =====
echo "=== SÉLECTION DES FICHIERS BAM ==="
BAM_FILES=(*.bam)

if [ ${#BAM_FILES[@]} -eq 0 ] || [ ! -f "${BAM_FILES[0]}" ]; then
    echo "Erreur: Aucun fichier BAM trouvé"
    exit 1
fi

echo "Fichiers BAM disponibles:"
for i in "${!BAM_FILES[@]}"; do
    echo "  [$i] ${BAM_FILES[$i]}"
done
echo ""

read -p "Traiter tous les fichiers BAM? (o/n): " process_all

if [[ $process_all == "o" || $process_all == "O" ]]; then
    SELECTED_BAM=("${BAM_FILES[@]}")
else
    echo "Entrez les numéros des fichiers (séparés par des espaces):"
    read -p "Exemple: 0 2 3 : " selected_numbers
    
    SELECTED_BAM=()
    for num in $selected_numbers; do
        if [[ $num =~ ^[0-9]+$ ]] && [ $num -ge 0 ] && [ $num -lt ${#BAM_FILES[@]} ]; then
            SELECTED_BAM+=("${BAM_FILES[$num]}")
        fi
    done
    
    if [ ${#SELECTED_BAM[@]} -eq 0 ]; then
        echo "Erreur: Aucun fichier sélectionné"
        exit 1
    fi
fi

echo "Fichiers sélectionnés: ${#SELECTED_BAM[@]}"
echo ""

read -p "Lancer le tri et l'indexation? (o/n): " confirm
if [[ $confirm != "o" && $confirm != "O" ]]; then
    echo "Opération annulée"
    exit 0
fi
echo ""

# ===== TRAITEMENT =====
echo "=== TRAITEMENT EN COURS ==="
SUCCESS=0

for BAM_INPUT in "${SELECTED_BAM[@]}"; do
    NUM=$(echo "$BAM_INPUT" | grep -o '[0-9]\+')
    BAM_SORTED="aln-se-sort.${NUM}.bam"
    
    echo "Traitement: $BAM_INPUT"
    
    # Trier le fichier BAM
    samtools sort "$BAM_INPUT" -o "$BAM_SORTED" 2>/dev/null
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Tri: $BAM_SORTED"
        
        # Indexer le fichier trié
        samtools index "$BAM_SORTED" 2>/dev/null
        
        if [ $? -eq 0 ]; then
            echo "  ✓ Index: ${BAM_SORTED}.bai"
            SUCCESS=$((SUCCESS+1))
        else
            echo "  ✗ Erreur lors de l'indexation"
        fi
    else
        echo "  ✗ Erreur lors du tri"
    fi
done

echo ""
echo "=== TERMINÉ ==="
echo "$SUCCESS fichier(s) BAM trié(s) et indexé(s) avec succès"
