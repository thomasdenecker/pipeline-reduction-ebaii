# pipeline-reduction-ebaii
Pipeline de réduction de données pour la formation EBAII 

## Objectif
Proposer un pipeline capable de prendre des fichiers FASTQ.GZ pour un organisme donné (à préciser dans les paramètres) pour arriver à des fichiers FASTQ.GZ réduits. 

L'objectif de créer des fichiers FASTQ.GZ réduits est de réduire les temps de calcul dans le cadre d'une formation afin que les apprenants puissent se concentrer sur la structuration du pipeline de traitement des données (cas d'usage traitement de données RNAseq dont le résultat est une table de comptage) et appliquer toutes étapes d'analyses qui seront présentées lors de la formation. 

## Inputs

Plusieurs fichiers :
- Des FASTQ.GZ
- Génome de référence (séquence uniquement)

## Output

- fastq réduits
- (Tableau de données réduites pouvant aboutir à une table de comptage ?) 


## Protocole

La réduction pour chaque fichier FASTQ.GZ est rélisée avec 2 sélections successive de reads : 
- les reads ayant mappés sur l'un des chromosomes (au choix et à définir dans les paramètres)
- les "n" premiers reads du FASTQ (autres que les reads déjà sélectionnés ci-dessus) afin de conserver une part de la variabilité des reads.

Etapes suivies :
- contrôle qualité des fichiers FASTQ.GZ (avec fastqc)
- si nécessaire, création de l'index pour le logiciel de mapping (avec hisat2)
- mapping des reads de chaque fichier FASTQ.GZ (avec hisat2)
- selection des reads mappés sur le chromosome choisi (avec samtools)
- selection des "n" premiers reads (avec bash)
- création des fichiers FASTQ.GZ réduits (avec samtools)
- mapping d'un des fichiers FASTQ.GZ réduits sur le génome pour vérification (avec hisat2)

Remarque : en cas de séquençage paired-end, le nombre de reads résultants ne sera pas exact à "n".

## Paramétrages

un fichier yaml a éditer afin de modifier les valeurs définies par défaut pour le jeu de données test

## Jeu de données test

Le jeu test suivant est choisi pour tester le bon fonctionement du pipeline.
Il est téléchargé automatiquement par le pipeline.
RNAseq paired-end : SRR4308679
genome : GCF_000214025.2

## Lieu d'utilisation de ce pipeline 

Ce pipeline sera utilisé sur le [cluster de l'IFB](https://www.france-bioinformatique.fr/clusters-ifb/)

implémentation du pipeline en snakemake, paramétrages fixés un fichier yaml.

(à compléter)
