# pipeline-reduction-ebaii
Pipeline de réduction de données pour la formation EBAII 

## Objectif
Proposer un pipeline snakemake capable de prendre des fichiers FASTQ.GZ pour un organisme donné (à préciser dans les paramètres) pour arriver à des fichiers FASTQ.GZ de taille réduite. 

L'intérêt à créer de tels fichiers est la réduction du temps de calcul pour que dans le cadre d'une formation, les apprenants puissent se concentrer sur la structuration du pipeline de traitement des données (cas d'usage initial : traitement de données RNAseq dont le résultat est une table de comptage) et appliquer toutes étapes d'analyses plutôt que d'attendre l'exécution du triatement à chaque progression. Un jeu de donnée rapide à traiter permet une démarche par essai-erreur dynamique. 

## Inputs

Plusieurs fichiers :
- Des FASTQ.GZ
- Génome de référence (séquence uniquement)

## Output

- fastq.gz réduits
- (Tableau de données réduites pouvant aboutir à une table de comptage ?) 


## Protocole

La réduction pour chaque fichier FASTQ.GZ est rélisée avec 2 sélections successive de reads en + des lignes du header : 
- les reads ayant mappés sur l'un des chromosomes (au choix et à définir dans les paramètres)
- les "n" premiers reads du FASTQ (autres que les reads déjà sélectionnés ci-dessus) afin de conserver une part de la variabilité des reads.

Etapes suivies :
- contrôle qualité des fichiers FASTQ.GZ (avec fastqc)
- création de l'index pour le logiciel de mapping (avec hisat2)
- mapping des reads de chaque fichier FASTQ.GZ (avec hisat2)
- selection des reads mappés sur le chromosome choisi (avec samtools)
- selection des "n" premiers reads (avec samtools + bash)
- création des fichiers FASTQ réduits (avec samtools)
- compression des fichiers FASTQ réduits (avec gzip)
- TODO mapping d'un des fichiers FASTQ.GZ réduits sur le génome pour vérification (avec hisat2)

Remarque : le nombre de reads résultants n'est pas exact à "n" + nombre de reads mappés sur le chromosome choisi (dépend du nombre de reads mappés sur le chromosome choisi inclus dans les "n" premiers reads ainsi que, en cas de séquençage paired-end, du nombre de reads mappés sans sa paire)

## Paramétrages

Les paramétrages pour le lancement sont réunis dans un fichier yaml. 
Editer afin de modifier les valeurs définies par défaut pour le jeu de données test :
- chemin d'accès à la séquence fna
- chemin d'accès au répertoire des fastq.gz à réduire
- nom du chromosome d'intérêt
- nombe "n" de premiers reads

## Jeu de données test

Voici un jeu de données choisi pour tester le bon fonctionement du pipeline :
- RNAseq paired-end : SRR4308679_R1.fastq.gz et SRR4308679_R2.fastq.gz (ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR430/009/SRR4308679/SRR4308679_1.fastq.gz, attention à renommer les "_1." en "_R1.")
- genome : GCF_000214025.2 [genome.fna](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/015/GCF_000214015.3_version_140606/GCF_000214015.3_version_140606_genomic.fna.gz)

## Lieu d'utilisation de ce pipeline 

Utilisé sur le [cluster de l'IFB](https://www.france-bioinformatique.fr/clusters-ifb/)

### sur un poste de travail unix + conda:

Créer un environnement conda dédié (cf. fichier yml) :
`conda env create -n RNASeqReduction -f ce_RNASeqReduction.yml`
Activer cet environement dédié :
`conda activate RNASeqReduction`
Créer le répertoire de données ("DataPOC" pour le jeu de donénes test) et y rapatrier la séquence du genome et les RNASeq
Adapter à votre arborescence de fichiers le fichier de paramétrage du pipeline (éditer data.yml)
Lancer le pipeline
`snakemake -s reduction.smk --configfile data.yml`

(à compléter)
