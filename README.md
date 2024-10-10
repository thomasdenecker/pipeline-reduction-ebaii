# pipeline-reduction-ebaii
Pipeline de réduction de données pour la formation EBAII 

## Objectif
Ce pipeline snakemake prend des fichiers FASTQ.GZ pour un organisme donné (à préciser dans les paramètres) et arrive à des fichiers FASTQ.GZ de taille réduite. 

La réduction des fichiers FASTQ permet un gain du temps de calcul pour que, dans le cadre d'une formation, les apprenants puissent se concentrer sur la structuration du pipeline de traitement des données (cas d'usage : traitement de données RNAseq dont le résultat est une table de comptage) plutôt que d'attendre l'exécution du traitement de chaque étape. Un jeu de donnée rapide à traiter permet une démarche par essai-erreur dynamique. Les résultats obtenus à partir des données réduites ne doivent servir qu'à l'entrainement et ne peuvent être utilisés par la suite.

## Inputs

Plusieurs fichiers :
- Des FASTQ.GZ
- Génome de référence (séquences au format fasta uniquement)

## Output

- fastq.gz réduits


## Protocole

La réduction pour chaque fichier FASTQ.GZ est rélisée avec 2 sélections successives de reads (en + des lignes du header) : 
- les "n" premiers reads ayant mappés sur l'un des chromosomes (l'identifiant du chr est à définir dans les paramètres ainsi que la valeur n, codée `nbAlignChr` dans le fichier de configuration du snakemake, cf. `data.yml`)
- les "m" premiers reads du FASTQ (autres que les reads déjà sélectionnés ci-dessus, m est codée `nbAlginRandom`) afin de conserver une part de la variabilité des reads.

Etapes suivies :
- décompression des fichiers FATSQ.GZ (avec gunzip, temporaire*)
- contrôle qualité des fichiers FASTQ (avec fastqc)
- création de l'index pour le logiciel de mapping (avec hisat2)
- mapping des reads de chaque fichier FASTQ.GZ (avec hisat2, temporaire*)
- selection des "n" reads mappés sur le chromosome choisi (avec samtools)
- selection des "m" premiers reads (avec samtools + bash, remarque**)
- création des fichiers FASTQ réduits (avec samtools)
- compression des fichiers FASTQ réduits (avec gzip)

*temporaire : les fichiers résultants ne sont pas conservés une fois qu'ils ne servent plus afin d'éviter la saturation de l'espace disque

**remarque : le nombre de reads résultants n'est pas exact à "n" + nombre de reads mappés sur le chromosome choisi (dépend du nombre de reads mappés sur le chromosome choisi inclus dans les "n" premiers reads ainsi que, en cas de séquençage paired-end, du nombre de reads mappés sans sa paire)

## Paramétrages

Les paramétrages pour le lancement sont réunis dans un fichier yaml (ex. `data.yml`) qu'il faut éditer afin de modifier les valeurs définies par défaut pour le jeu de données test :
- chemin d'accès à la séquence fna
- chemin d'accès au répertoire des fastq.gz à réduire
- nom du répertoire des fastq.gz réduits
- identifiant du chromosome d'intérêt
- nombres "n" et "m" de premiers reads à conserver

## Jeu de données test

Le jeu de données choisi pour tester le bon fonctionement du pipeline doit être téléchargé :
- RNAseq paired-end (attention à renommer les "_1." en "_R1.") SRR4308679_R1.fastq.gz et SRR4308679_R2.fastq.gz : ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR430/009/SRR4308679/SRR4308679_1.fastq.gz ; ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR430/009/SRR4308679/SRR4308679_2.fastq.gz 
- genome GCF_000214025.2 : [genome.fna](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/015/GCF_000214015.3_version_140606/GCF_000214015.3_version_140606_genomic.fna.gz)

## Utilisation du pipeline 

### Préparation des fichiers

+ Créer le répertoire de données ("DataPOC" pour le jeu de données test et y rapatrier la séquence du genome et les RNASeq, voir ci-dessus "Jeu de données test")
+ Adapter l'arborescence des fichiers dans le fichier de paramétrage du pipeline (éditer le fichier `data.yml` ; pas d'adaptation dans le cas des données "DataPOC")

### sur le [cluster de l'IFB](https://www.france-bioinformatique.fr/clusters-ifb/)

+ Charger les modules nécessaires : `module load snakemake fastqc samtools hisat2`
+ S'il n'existe pas déjà un fichier profile "slurm", copier sous `~/.config/snakemake/slurm/config.yaml` le fichier donné en exemple (`ifb_slurm_profile.yaml`) 
+ Se placer dans l'espace projet : `cd /shared/projects/... `
+ Lancer le pipeline : `sbatch snakemake --profile slurm --jobs 4 --cores 4 -s reduction.smk --configfile data.yml `

#### Versions des logiciels

+ snakemake/6.5.0 (testé entre 5.19.2 et 7.7.0)
+ fastqc/v0.11.9 
+ samtools/1.15.1 
+ hisat2/2.2.1

### sur un poste de travail unix + conda:

+ Créer un environnement conda dédié (cf. `ce_RNASeqReduction.yml`) : `conda env create -n RNASeqReduction -f ce_RNASeqReduction.yml`
+ Activer l'environement dédié : `conda activate RNASeqReduction`
+ Lancer le pipeline : `snakemake -s reduction.smk --configfile data.yml --useconda`


## Améliorations

+ mapping d'un des fichiers FASTQ.GZ réduits sur le génome pour vérification (avec hisat2)
+ authoriser les index hisat2 "larges" (génomes de taille > à 4 billion de nucléotides):  les fichiers d'index se terminent par ht2l au lieu de ht2
+ libérer la contrainte "R1/R2.fastq.gz" pour les noms des fichiers pairés : authoriser "1/2" seulement et une variation du suffixe.

