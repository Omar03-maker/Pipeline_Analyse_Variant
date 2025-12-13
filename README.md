# PIPELINE D'ANALYSE DE VARIANTS 
Il s'agit dune pipeline automatisé pour l'analyse et la visualisation de variants à partir de données de sequencage génomique (NGS). Dans notre exemple nous avons utilisé la sequence de reference de l'espece *Mycobacterium tuberculosis*. Cette pipeline permet notamment de : 
- Quantifier les reads qui s'alignent sur le génome de *M.uberculosis* (Reference : AL123456)
- Identifier les variants et les SNPs présent
- Faire le filtrage et la visualisation des variants avec IGV

## 📝 Notes importantes
- **Pour votre analyse remplacer les fichiers fastq et fasta par vos propres fichiers** : Ce sont des données d'exemples pour la pipeline
- **Installer tous les outils et fichiers dans les repertoires correspondant pour une meilleure efficacité** 
- **Tous les scripts sont interactifs** : Cette pipeline vous guide étape par étape jusqu'à la visualisation de vos variants
- **Le filtrage DP ≥ 3** est le plus recommandé pour réduire les faux positifs
- **L'ordre d'exécution** des scripts est important (respectez la séquence indiquée)


## Outils nécessaires
- BWA pour l'alignement de séquences
- Samtools pour la manipulation de fichiers BAM/SAM
- Bcftools pour la génération de fichiers BCF/VCF
- Kallisto pour la quantification des reads
- SnpsefF.jar pour le filtrage des SNPs


## 📦 Installation des outils
### 1. Installation BWA, Samtools et Bcftools
#### Methode 1 : En ligne de commande 
sudo apt-get update

sudo apt-get install bwa

sudo apt-get install samtools

sudo apt-get install bcftools
#### Methode 2 : Via le dépôt github 
Uitliser le lien pour recuperer les fichiers correspondant : https://github.com/COMBINE-lab/salmon/releases

### 2. Installation Kallisto
#### Méthode 1 : 
sudo apt-get install kallisto
#### Méthode 2 : Compilation depuis GitHub
cd ~

git clone https://github.com/pachterlab/kallisto

cd kallisto

mkdir build

cd build

cmake ..

make

sudo make install
### 3. Télécharger snpEff
- wget https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
- unzip snpEff_latest_core.zip``


## Utilisation de la pipeline
### Partie A : Quantification avec Kallisto_tools
- Aller au niveau du repertoire kallisto_tools en utilisant la commande suivante : cd kallisto_tools/
- Puis executer le fichier sh suivant : ./Taux_Mapping.sh 
#### Le script vous guidera dans votre demarche :
- Sélection du fichier FASTA de référence
- Création de l'index Kallisto si nécessaire
- Choix des fichiers FASTQ à analyser
- Génération automatiquement des statistiques d'alignement
#### Résultats attendus :
- Nombre total de lectures traitées
- Nombre de lectures pseudoalignées
- Taux de mapping (pourcentage d'alignement)
 
Les résultats seront ensuite transférer dans le dossier kallisto_output/

### Partie B : Analyse des SNPs (Détection de variants)
#### Étape 1 : Conversion FASTQ → BAM
- Aller au niveau du repertoire sam_tools en utilisant la commande suivante : cd salmon_tools/
- Puis executer le fichier sh suivant : ./FastQ_TO_Bam.sh
#### Étape 2 : Tri et indexation des BAM
- Toujours dans le meme repertoire executer le fichier sh suivant :./Bam_Filtered.sh
Cette étape :
- Trie les fichiers BAM par coordonnées génomiques
- Indexe les fichiers BAM triés
##### Visualiser des fichiers
- Pour les fichiers BAM non trié effectuer la commande suivante : samtools view aln-se.1.bam | head
- Pour les fichiers BAM non trié effectuer la commande suivante : samtools view aln-se-sort.1.bam | head
#### Étape 3 : Génération des fichiers BCF
- Toujours dans le meme repertoire executer le fichier sh suivant : ./Bam_To_Bcf.sh
#### Étape 4 : Génération et filtrage des VCF
- Toujours dans le meme repertoire executer le fichier sh suivant : ./Bcf_To_Vcf_Filtered.sh
Cette étape :
- Convertit BCF → VCF
- Filtre les SNPs avec une profondeur de lecture ≥ 3 (paramètre modifiable)


## 🔍 Visualisation des résultats
### Visualisation avec IGV (Integrative Genomics Viewer)
#### 1. Installer IGV
- Téléchargez IGV depuis le site officiel :
  
https://software.broadinstitute.org/software/igv/download

- Ou installation via ligne de commande linus :
  
wget https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_Linux_2.16.2_WithJava.zip

unzip IGV_Linux_2.16.2_WithJava.zip

cd IGV_Linux_2.16.2

./igv.sh

#### 2. Charger les données dans IGV
##### a) Charger le génome de référence :
- Menu : Genomes → Load Genome from File
- Sélectionner : AL123456.3.fasta
##### b) Charger les fichiers d'alignement et de variants :
- Menu : File → Load from File
- Sélectionner : aln-se-sort.1.bam (fichier BAM trié)
- Menu : File → Load from File
- Sélectionner : aln-se-sort-filter.1.vcf (variants filtrés)

#### 3. Explorer les variants
- La piste **BAM** montre la couverture (coverage) et l'alignement des reads
- La piste **VCF** affiche uniquement les SNPs avec DP ≥ 3
- Zoomez sur une région pour voir les variants en détail
- Les SNPs présents sur moins de 3 reads ne sont pas affichés (filtrage)


## 📊 Résultats attendus
### Détection de variants
Les fichiers VCF filtrés contiennent les SNPs détectés avec :
- Position chromosomique
- Nucléotide de référence
- Nucléotide alternatif (variant)
- Profondeur de lecture (DP)
- Qualité du variant

# Auteur
- El Hadji Omar Dia
- GitHub : @Omar03-maker
- Mail : elhadjiomardia@esp.sn

# Si vous trouver ce projet intéressant, n'hésitez pas à lui donner une étoile ⭐ !
