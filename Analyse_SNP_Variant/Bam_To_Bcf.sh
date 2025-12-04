#!/bin/bash

echo "=== GÉNÉRATION DE FICHIERS BCF À PARTIR DE BAM TRIÉS ==="
echo ""

# Vérifier que bcftools et samtools sont installés
for tool in bcftools samtools; do
    if ! command -v $tool &> /dev/null; then
        echo "Erreur: $tool n'est pas installé"
        exit 1
    fi
done

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

# ===== SÉLECTION DES FICHIERS BAM =====
echo "=== SÉLECTION DES FICHIERS BAM ==="
echo "Types de fichiers BAM disponibles:"
echo "  [1] Fichiers BAM triés (aln-se-sort.*.bam) - RECOMMANDÉ"
echo "  [2] Tous les fichiers BAM"
read -p "Votre choix (1/2): " bam_type

if [ "$bam_type" == "1" ]; then
    BAM_FILES=(aln-se-sort.*.bam)
elif [ "$bam_type" == "2" ]; then
    BAM_FILES=(*.bam)
else
    echo "Choix invalide"
    exit 1
fi

if [ ${#BAM_FILES[@]} -eq 0 ] || [ ! -f "${BAM_FILES[0]}" ]; then
    echo "Erreur: Aucun fichier BAM trouvé"
    exit 1
fi

echo ""
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

read -p "Lancer la génération des fichiers BCF? (o/n): " confirm
if [[ $confirm != "o" && $confirm != "O" ]]; then
    echo "Opération annulée"
    exit 0
fi
echo ""

# ===== TRAITEMENT =====
echo "=== GÉNÉRATION DES FICHIERS BCF ==="
SUCCESS=0

for BAM_SORTED in "${SELECTED_BAM[@]}"; do
    NUM=$(echo "$BAM_SORTED" | grep -o '[0-9]\+')
    BCF_OUT="aln-se-sort.${NUM}.bcf"
    
    echo "Traitement: $BAM_SORTED"
    
    # Vérifier que le BAM est indexé
    if [ ! -f "${BAM_SORTED}.bai" ]; then
        echo "  Indexation du fichier BAM..."
        samtools index "$BAM_SORTED" 2>/dev/null
    fi
    
    # Générer le fichier BCF
    bcftools mpileup -f "$FASTA_REF" "$BAM_SORTED" -o "$BCF_OUT" 2>/dev/null
    
    if [ $? -eq 0 ]; then
        echo "  ✓ BCF généré: $BCF_OUT"
        SUCCESS=$((SUCCESS+1))
    else
        echo "  ✗ Erreur lors de la génération du BCF"
    fi
done

echo ""
echo "=== TERMINÉ ==="
echo "$SUCCESS fichier(s) BCF généré(s) avec succès"
