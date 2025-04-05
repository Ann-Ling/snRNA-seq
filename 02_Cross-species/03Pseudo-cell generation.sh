# Combine the UMI matrix of KC between Mpha and Amel
python universal.combine.KC.umi.py Dros.midgut.umi.geneMatrix.txt AM.umi.geneMatrix.txt Dros_AM.ort.genes Dros AM > Dros_AM.midgut.umi.geneMatrix.txt

cat Dros.pheno.cell.txt AM.pheno.cell.txt | grep -v 'Sample_ID'|sed '1iSample_ID\tStudy_ID\tCelltype' > Dros_AM.midgut.pheno.cell.txt

# Generate pseudo cell umi Matrix
## we will get：merged.pheno.cell.txt, merged.exp.txt, and merged.umi.txt
perl random.combine.pl Dros_AM.midgut.umi.geneMatrix.txt Dros_AM.midgut.pheno.cell.txt
