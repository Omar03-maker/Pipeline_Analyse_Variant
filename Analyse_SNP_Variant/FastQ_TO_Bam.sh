#!/bin/bash

echo "=== CONVERSION FASTQ -> BAM AVEC BWA ET SAMTOOLS ==="
echo ""

# Vérifier que BWA et Samtools sont installés
if ! command -v bwa &> /dev/null; then
    echo "Erreur: bwa n'est pas installé."
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "Erreur: samtools n'est pas installé."
    exit 1
fi

# ===== SÉLECTION DU FICHIER FASTA =====
echo "=== SÉLECTION DU FICHIER FASTA DE RÉFÉRENCE ==="
FASTA_FILES=(*.fasta)

if [ ${#FASTA_FILES[@]} -eq 0 ] || [ ! -f "${FASTA_FILES[0]}" ]; then
    echo "Erreur: Aucun fichier FASTA trouvé"
    exit 1
elif [ ${#FASTA_FILES[@]} -eq 1 ]; then
    echo "Fichier FASTA trouvé: ${FASTA_FILES[0]}"
    read -p "Utiliser ce fichier? (o/n): " confirm
    if [[ $confirm != "o" && $confirm != "O" ]]; then
        echo "Opération annulée"
        exit 0
    fi
    FASTA_REF="${FASTA_FILES[0]}"
else
    echo "Fichiers FASTA disponibles:"
    for i in "${!FASTA_FILES[@]}"; do
        echo "  [$i] ${FASTA_FILES[$i]}"
    done
    read -p "Numéro du fichier à utiliser: " choice
    if [[ $choice =~ ^[0-9]+$ ]] && [ $choice -ge 0 ] && [ $choice -lt ${#FASTA_FILES[@]} ]; then
        FASTA_REF="${FASTA_FILES[$choice]}"
    else
        echo "Choix invalide"
        exit 1
    fi
fi
echo ""

# ===== VÉRIFICATION/CRÉATION DE L'INDEX BWA =====
echo "=== VÉRIFICATION DE L'INDEX BWA ==="
if [ ! -f "${FASTA_REF}.bwt" ]; then
    echo "Index BWA non trouvé. Création de l'index..."
    bwa index "$FASTA_REF"
    if [ $? -eq 0 ]; then
        echo "✓ Index créé avec succès"
    else
        echo "✗ Erreur lors de la création de l'index"
        exit 1
    fi
else
    echo "Index BWA existant trouvé"
    read -p "Régénérer l'index? (o/n): " regen
    if [[ $regen == "o" || $regen == "O" ]]; then
        echo "Régénération de l'index..."
        bwa index "$FASTA_REF"
        if [ $? -eq 0 ]; then
            echo "✓ Index régénéré avec succès"
        else
            echo "✗ Erreur lors de la régénération"
            exit 1
        fi
    fi
fi
echo ""

# ===== SÉLECTION DES FICHIERS FASTQ =====
echo "=== SÉLECTION DES FICHIERS FASTQ ==="
FASTQ_FILES=(*.fastq)

if [ ${#FASTQ_FILES[@]} -eq 0 ] || [ ! -f "${FASTQ_FILES[0]}" ]; then
    echo "Erreur: Aucun fichier FASTQ trouvé"
    exit 1
fi

echo "Fichiers FASTQ disponibles:"
for i in "${!FASTQ_FILES[@]}"; do
    echo "  [$i] ${FASTQ_FILES[$i]}"
done
echo ""

read -p "Traiter tous les fichiers FASTQ? (o/n): " process_all

if [[ $process_all == "o" || $process_all == "O" ]]; then
    SELECTED_FASTQ=("${FASTQ_FILES[@]}")
else
    echo "Entrez les numéros des fichiers (séparés par des espaces):"
    read -p "Exemple: 0 2 3 : " selected_numbers
    
    SELECTED_FASTQ=()
    for num in $selected_numbers; do
        if [[ $num =~ ^[0-9]+$ ]] && [ $num -ge 0 ] && [ $num -lt ${#FASTQ_FILES[@]} ]; then
            SELECTED_FASTQ+=("${FASTQ_FILES[$num]}")
        fi
    done
    
    if [ ${#SELECTED_FASTQ[@]} -eq 0 ]; then
        echo "Erreur: Aucun fichier sélectionné"
        exit 1
    fi
fi

echo "Fichiers sélectionnés: ${#SELECTED_FASTQ[@]}"
echo ""

read -p "Lancer la conversion? (o/n): " confirm
if [[ $confirm != "o" && $confirm != "O" ]]; then
    echo "Opération annulée"
    exit 0
fi
echo ""

# ===== TRAITEMENT =====
echo "=== CONVERSION EN COURS ==="
SUCCESS=0

for FASTQ_FILE in "${SELECTED_FASTQ[@]}"; do
    NUM=$(echo "$FASTQ_FILE" | grep -o '[0-9]\+')
    SAM_OUT="aln-se.${NUM}.sam"
    BAM_OUT="aln-se.${NUM}.bam"
    
    echo "Traitement: $FASTQ_FILE"
    
    # Alignement avec BWA MEM
    bwa mem "$FASTA_REF" "$FASTQ_FILE" > "$SAM_OUT" 2>/dev/null
    
    # Conversion SAM en BAM
    samtools view -S -b "$SAM_OUT" > "$BAM_OUT" 2>/dev/null
    
    if [ $? -eq 0 ]; then
        echo "  ✓ BAM généré: $BAM_OUT"
        SUCCESS=$((SUCCESS+1))
        rm "$SAM_OUT"  # Supprimer le SAM intermédiaire
    else
        echo "  ✗ Erreur pour $FASTQ_FILE"
    fi
done

echo ""
echo "=== TERMINÉ ==="
echo "$SUCCESS fichier(s) BAM généré(s) avec succès"
