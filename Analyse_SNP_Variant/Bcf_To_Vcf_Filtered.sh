#!/bin/bash

echo "=== GÉNÉRATION DE FICHIERS VCF FILTRÉS À PARTIR DE BCF ==="
echo ""

# Vérifier que bcftools et java sont installés
for tool in bcftools java; do
    if ! command -v $tool &> /dev/null; then
        echo "Erreur: $tool n'est pas installé"
        exit 1
    fi
done

# Vérifier que SnpSift.jar existe
if [ ! -f "Snpsift.jar" ]; then
    echo "Erreur: Le fichier 'Snpsift.jar' est introuvable dans le répertoire courant"
    exit 1
fi

# ===== SÉLECTION DES FICHIERS BCF =====
echo "=== SÉLECTION DES FICHIERS BCF ==="
BCF_FILES=(*.bcf)

if [ ${#BCF_FILES[@]} -eq 0 ] || [ ! -f "${BCF_FILES[0]}" ]; then
    echo "Erreur: Aucun fichier BCF trouvé"
    exit 1
fi

echo "Fichiers BCF disponibles:"
for i in "${!BCF_FILES[@]}"; do
    echo "  [$i] ${BCF_FILES[$i]}"
done
echo ""

read -p "Traiter tous les fichiers BCF? (o/n): " process_all

if [[ $process_all == "o" || $process_all == "O" ]]; then
    SELECTED_BCF=("${BCF_FILES[@]}")
else
    echo "Entrez les numéros des fichiers (séparés par des espaces):"
    read -p "Exemple: 0 2 3 : " selected_numbers
    
    SELECTED_BCF=()
    for num in $selected_numbers; do
        if [[ $num =~ ^[0-9]+$ ]] && [ $num -ge 0 ] && [ $num -lt ${#BCF_FILES[@]} ]; then
            SELECTED_BCF+=("${BCF_FILES[$num]}")
        fi
    done
    
    if [ ${#SELECTED_BCF[@]} -eq 0 ]; then
        echo "Erreur: Aucun fichier sélectionné"
        exit 1
    fi
fi

echo "Fichiers sélectionnés: ${#SELECTED_BCF[@]}"
echo ""

# ===== PARAMÈTRES DE FILTRAGE =====
echo "=== PARAMÈTRES DE FILTRAGE ==="
read -p "Profondeur minimale de lecture (DP) [défaut: 3]: " min_dp
min_dp=${min_dp:-3}

echo "Filtre sélectionné: DP >= $min_dp"
echo ""

read -p "Lancer la conversion et le filtrage? (o/n): " confirm
if [[ $confirm != "o" && $confirm != "O" ]]; then
    echo "Opération annulée"
    exit 0
fi
echo ""

# ===== TRAITEMENT =====
echo "=== CONVERSION ET FILTRAGE EN COURS ==="
SUCCESS=0

for BCF_FILE in "${SELECTED_BCF[@]}"; do
    NUM=$(echo "$BCF_FILE" | grep -o '[0-9]\+')
    VCF_FILE="aln-se-sort.${NUM}.vcf"
    VCF_FILTERED="aln-se-sort-filter.${NUM}.vcf"
    
    echo "Traitement: $BCF_FILE"
    
    # Générer le fichier VCF
    bcftools call -v -m "$BCF_FILE" > "$VCF_FILE" 2>/dev/null
    
    if [ $? -eq 0 ]; then
        echo "  ✓ VCF généré: $VCF_FILE"
        
        # Filtrer les SNPs avec profondeur >= min_dp
        java -jar Snpsift.jar filter "(DP >= $min_dp)" -f "$VCF_FILE" > "$VCF_FILTERED" 2>/dev/null
        
        if [ $? -eq 0 ]; then
            echo "  ✓ VCF filtré généré: $VCF_FILTERED"
            SUCCESS=$((SUCCESS+1))
        else
            echo "  ✗ Erreur lors du filtrage VCF"
        fi
    else
        echo "  ✗ Erreur lors de la génération du VCF"
    fi
done

echo ""
echo "=== TERMINÉ ==="
echo "$SUCCESS fichier(s) VCF filtré(s) généré(s) avec succès"
