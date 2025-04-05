#Creat environment
conda create -n diamond 

#Install diamond
conda install -c bioconda diamond 

#Activate
conda activate diamond
diamond makedb --in honeybee.faa -d Apis_mellifera.dmnd --threads 8 #Apis mellifera database
diamond makedb --in Drosophila_melanogaster.faa -d Drosophila_melanogaster.dmnd --threads 8 #Drosophila melanogaster database

#Start homology comparison
diamond blastp --db Drosophila_melanogaster.dmnd --query honeybee.faa -o Dros_AM_protein.txt --threads 8

