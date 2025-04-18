# step1 grn 
cd /data2/liuhuiling/gut_snRNA/10SCENIC/02ileum
pyscenic grn --num_workers 10 --output grn.tsv --method grnboost2 ileum_scenic.loom TF_list.txt

nohup pyscenic grn --num_workers 6 --output grn.tsv --method grnboost2 ileum_scenic.loom TF_list.txt > pyscenic_output.log 2>&1 &

# step2 ctx
pyscenic ctx \
  grn.tsv motifs.regions_vs_motifs.rankings.feather \
  --annotations_fname motifs-v9-nr.bee-m0.001-o0.1.tbl \
  --expression_mtx_fname ileum_scenic.loom \
  --mode "dask_multiprocessing" \
  --output ctx_2.csv \
  --num_workers 10 \
  --mask_dropouts \
  --rank_threshold 1000 \
  --auc_threshold 0.02 \
  --nes_threshold 2.0

# step3 AUCell  
pyscenic aucell \
  ileum_scenic.loom \
  ctx_2.csv \
  --output aucell_2.loom \
  --num_workers 15
  
