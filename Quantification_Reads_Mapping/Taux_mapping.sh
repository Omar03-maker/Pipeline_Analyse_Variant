#!/bin/bash

# Répertoire courant
WORK_DIR=$(pwd)
OUTPUT_DIR="kallisto_output"

echo "=== SCRIPT KALLISTO AUTOMATIQUE ==="
echo ""

# ===== SÉLECTION DU FICHIER FASTA =====
echo "=== SÉLECTION DU FICHIER FASTA ==="
FASTA_FILES=(*.fasta)

if [ ${#FASTA_FILES[@]} -eq 0 ] || [ ! -f "${FASTA_FILES[0]}" ]; then
    echo "Erreur: Aucun fichier FASTA trouvé dans le répertoire"
    exit 1
elif [ ${#FASTA_FILES[@]} -eq 1 ]; then
    echo "Fichier FASTA trouvé: ${FASTA_FILES[0]}"
    read -p "Voulez-vous utiliser ce fichier? (o/n): " confirm
    if [[ $confirm != "o" && $confirm != "O" ]]; then
        echo "Opération annulée"
        exit 0
    fi
    FASTA_FILE="${FASTA_FILES[0]}"
else
    echo "Plusieurs fichiers FASTA trouvés:"
    for i in "${!FASTA_FILES[@]}"; do
        echo "  [$i] ${FASTA_FILES[$i]}"
    done
    read -p "Entrez le numéro du fichier à utiliser: " fasta_choice
    if [[ $fasta_choice =~ ^[0-9]+$ ]] && [ $fasta_choice -ge 0 ] && [ $fasta_choice -lt ${#FASTA_FILES[@]} ]; then
        FASTA_FILE="${FASTA_FILES[$fasta_choice]}"
        echo "Fichier sélectionné: $FASTA_FILE"
    else
        echo "Choix invalide"
        exit 1
    fi
fi
echo ""

# ===== SÉLECTION/CRÉATION DE L'INDEX =====
echo "=== SÉLECTION DE L'INDEX ==="
INDEX_FILES=(*.index)

if [ ${#INDEX_FILES[@]} -eq 0 ] || [ ! -f "${INDEX_FILES[0]}" ]; then
    echo "Aucun fichier index trouvé"
    BASE_NAME=$(basename "$FASTA_FILE" .fasta)
    INDEX_FILE="${BASE_NAME}.index"
    echo "Création d'un nouvel index: $INDEX_FILE"
    read -p "Continuer? (o/n): " confirm
    if [[ $confirm != "o" && $confirm != "O" ]]; then
        echo "Opération annulée"
        exit 0
    fi
elif [ ${#INDEX_FILES[@]} -eq 1 ]; then
    echo "Fichier index trouvé: ${INDEX_FILES[0]}"
    echo "Options:"
    echo "  [1] Utiliser cet index existant"
    echo "  [2] Régénérer l'index"
    read -p "Votre choix (1/2): " index_choice
    
    if [ "$index_choice" == "1" ]; then
        INDEX_FILE="${INDEX_FILES[0]}"
        REGENERATE_INDEX=false
    elif [ "$index_choice" == "2" ]; then
        INDEX_FILE="${INDEX_FILES[0]}"
        REGENERATE_INDEX=true
    else
        echo "Choix invalide"
        exit 1
    fi
else
    echo "Plusieurs fichiers index trouvés:"
    for i in "${!INDEX_FILES[@]}"; do
        echo "  [$i] ${INDEX_FILES[$i]}"
    done
    echo "  [n] Créer un nouvel index"
    read -p "Entrez votre choix: " index_choice
    
    if [[ $index_choice =~ ^[0-9]+$ ]] && [ $index_choice -ge 0 ] && [ $index_choice -lt ${#INDEX_FILES[@]} ]; then
        INDEX_FILE="${INDEX_FILES[$index_choice]}"
        echo "Index sélectionné: $INDEX_FILE"
        read -p "Voulez-vous régénérer cet index? (o/n): " regen
        if [[ $regen == "o" || $regen == "O" ]]; then
            REGENERATE_INDEX=true
        else
            REGENERATE_INDEX=false
        fi
    elif [[ $index_choice == "n" || $index_choice == "N" ]]; then
        BASE_NAME=$(basename "$FASTA_FILE" .fasta)
        INDEX_FILE="${BASE_NAME}.index"
        REGENERATE_INDEX=true
        echo "Création d'un nouvel index: $INDEX_FILE"
    else
        echo "Choix invalide"
        exit 1
    fi
fi
echo ""

# ===== CRÉATION/RÉGÉNÉRATION DE L'INDEX =====
if [ ! -f "$INDEX_FILE" ] || [ "$REGENERATE_INDEX" == true ]; then
    echo "=== GÉNÉRATION DE L'INDEX ==="
    kallisto index -i "$INDEX_FILE" "$FASTA_FILE"
    
    if [ $? -eq 0 ]; then
        echo "✓ Index créé/régénéré avec succès"
    else
        echo "✗ Erreur lors de la génération de l'index"
        exit 1
    fi
    echo ""
fi

# ===== SÉLECTION DES FICHIERS FASTQ =====
echo "=== SÉLECTION DES FICHIERS FASTQ ==="
FASTQ_FILES=(*.fastq)

if [ ${#FASTQ_FILES[@]} -eq 0 ] || [ ! -f "${FASTQ_FILES[0]}" ]; then
    echo "Erreur: Aucun fichier FASTQ trouvé dans le répertoire"
    exit 1
fi

echo "Fichiers FASTQ disponibles:"
for i in "${!FASTQ_FILES[@]}"; do
    echo "  [$i] ${FASTQ_FILES[$i]}"
done
echo ""

read -p "Voulez-vous traiter tous les fichiers FASTQ? (o/n): " process_all

if [[ $process_all == "o" || $process_all == "O" ]]; then
    SELECTED_FASTQ=("${FASTQ_FILES[@]}")
    echo "Tous les fichiers seront traités (${#SELECTED_FASTQ[@]} fichiers)"
else
    echo "Entrez les numéros des fichiers à traiter (séparés par des espaces):"
    read -p "Exemple: 0 2 3 : " selected_numbers
    
    SELECTED_FASTQ=()
    for num in $selected_numbers; do
        if [[ $num =~ ^[0-9]+$ ]] && [ $num -ge 0 ] && [ $num -lt ${#FASTQ_FILES[@]} ]; then
            SELECTED_FASTQ+=("${FASTQ_FILES[$num]}")
        else
            echo "Attention: Numéro invalide ignoré: $num"
        fi
    done
    
    if [ ${#SELECTED_FASTQ[@]} -eq 0 ]; then
        echo "Erreur: Aucun fichier valide sélectionné"
        exit 1
    fi
    
    echo "Fichiers sélectionnés:"
    for file in "${SELECTED_FASTQ[@]}"; do
        echo "  - $file"
    done
fi
echo ""

read -p "Confirmer le lancement de l'analyse? (o/n): " final_confirm
if [[ $final_confirm != "o" && $final_confirm != "O" ]]; then
    echo "Opération annulée"
    exit 0
fi
echo ""

# ===== CRÉATION DU RÉPERTOIRE DE SORTIE =====
mkdir -p "$OUTPUT_DIR"

# ===== TRAITEMENT DES FICHIERS =====
echo "=== TRAITEMENT EN COURS ==="
echo "Traitement de ${#SELECTED_FASTQ[@]} fichier(s) FASTQ..."
echo ""

for FASTQ_FILE in "${SELECTED_FASTQ[@]}"
do
    BASE_NAME=$(basename "$FASTQ_FILE" .fastq)
    
    kallisto quant \
        -i "$INDEX_FILE" \
        -o "$OUTPUT_DIR/${BASE_NAME}_output" \
        --single \
        -l 200 \
        -s 20 \
        "$FASTQ_FILE" > /dev/null 2>&1
done

# ===== AFFICHAGE DES RÉSULTATS =====
echo "=== RÉSULTATS D'ALIGNEMENT ==="
echo ""

for FASTQ_FILE in "${SELECTED_FASTQ[@]}"
do
    BASE_NAME=$(basename "$FASTQ_FILE" .fastq)
    JSON_FILE="$OUTPUT_DIR/${BASE_NAME}_output/run_info.json"
    
    if [ -f "$JSON_FILE" ]; then
        if command -v jq &> /dev/null; then
            PROCESSED=$(jq -r '.n_processed' "$JSON_FILE")
            ALIGNED=$(jq -r '.n_pseudoaligned' "$JSON_FILE")
            PERCENTAGE=$(echo "scale=2; 100 * $ALIGNED / $PROCESSED" | bc)
            
            echo "Fichier: $FASTQ_FILE"
            echo "  Lectures traitées       : $PROCESSED"
            echo "  Lectures pseudoalignées : $ALIGNED"
            echo "  Pourcentage d'alignement: ${PERCENTAGE}%"
            echo ""
        else
            PROCESSED=$(grep -oP '"n_processed":\s*\K[0-9]+' "$JSON_FILE")
            ALIGNED=$(grep -oP '"n_pseudoaligned":\s*\K[0-9]+' "$JSON_FILE")
            PERCENTAGE=$(echo "scale=2; 100 * $ALIGNED / $PROCESSED" | bc)
            
            echo "Fichier: $FASTQ_FILE"
            echo "  Lectures traitées       : $PROCESSED"
            echo "  Lectures pseudoalignées : $ALIGNED"
            echo "  Pourcentage d'alignement: ${PERCENTAGE}%"
            echo ""
        fi
    else
        echo "Fichier: $FASTQ_FILE - Erreur: pas de résultats"
        echo ""
    fi
done

echo "=== TRAITEMENT TERMINÉ ==="
echo "Les résultats sont dans le dossier: $OUTPUT_DIR"
