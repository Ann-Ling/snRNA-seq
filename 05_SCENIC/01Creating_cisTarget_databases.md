we created the cisTarget database as described in [Create cisTarget databases](https://github.com/aertslab/create_cisTarget_databases), which included creating a conda environment, installing necessary packages, and creating cisTarget motif databases. 
# Environment configuration
Please refer to this workflow:[Create cisTarget databases](https://github.com/aertslab/create_cisTarget_databases)

# Document preparation
```shell
# Set the GTF file paths
gtf=/data2/liuhuiling/SCENIC_preperation/Creating_cisTarget_motif_databases/Genome_assembly_Amel_HAv3.1.gtf

# Step 1: Extract the chromosome location and gene name of the gene and save it to tss.bed
awk 'BEGIN{FS=OFS="\t"}$1!~/^#/ && $3=="gene"{
    end=$4+1;
    # Extract the gene name, the gene name is in the gene field
    if (match($9, /gene "([^"]+)"/, m)) {
        gene_name = m[1];
    }
    # Output chromosome, start and end position, gene name, direction, gene type
    print $1, $4, end, gene_name, $6, $7
}' $gtf | sed 's/"//g' > tss.bed

# Step 2: Expand the gene region according to the TSS coordinates and save to gene_region.bed
awk 'BEGIN{FS=OFS="\t"}{
    if ($2 >= 1999) {
        # If the starting position is greater than or equal to 1999, then expand the front and back areas
        start = $2 - 1999;
        end = $2 + 4999;
        print $1, start, end, $4
    } else {
        # If the starting position is less than 1999, then expand from 0
        end = $2 + 4999;
        print $1, "0", end, $4
    }
}' tss.bed > gene_region.bed
```

# Start building database
```shell
create_cistarget_databases_dir=/data2/liuhuiling/create_cisTarget_databases
# FASTA file with sequences per region IDs / gene IDs.
fasta_filename=/data2/liuhuiling/SCENIC_preperation/Creating_cisTarget_motif_databases/gene_region.fasta
# Directory with motifs in Cluster-Buster format.
motifs_dir=/data2/liuhuiling/SCENIC_preperation/Creating_cisTarget_motif_databases/motif_cb
# File with motif IDs (base name of motif file in ${motifs_dir}).
motifs_list_filename=/data2/liuhuiling/SCENIC_preperation/Creating_cisTarget_motif_databases/motif_cb_name_list.txt
# cisTarget motif database output prefix.
db_prefix=motifs
nbr_threads=15

"${create_cistarget_databases_dir}/create_cistarget_motif_databases.py" \
    -f "${fasta_filename}" \
    -M "${motifs_dir}" \
    -m "${motifs_list_filename}" \
    -o "${db_prefix}" \
    -t "${nbr_threads}"
```
