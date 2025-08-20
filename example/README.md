# Script usage (detail)
## ST_analysis.py

```
python3 ST_analysis.py -i /mnt/VisiumHD_Colon/binned_outputs/square_008um \
                -o /mnt/output/write_first \
                -p /mnt/VisiumHD_Colon/8um_squares_annotation_pathologist.csv \
                -t Neoplasm \
                --auto-resolution \
                --cpus 4 \
                --annotation \
                --gsea \
                --survival \
                --tcga /stDrug/data/TCGA \
                --id TCGA-COAD \
                --not_treated
```
## Sub-cluster analysis

```
python3 ST_analysis.py -i /mnt/output/write_first/scanpyobj.h5ad \
                -f h5ad \
                -o /mnt/output/write_sub \
                -r 0.6 \
                --cpus 4 \
                -c '0,2' \
                --annotation \
                --gsea \
                --survival \
                --tcga /stDrug/data/TCGA \
                --not_treated
```

## drug response prediction.py
```
# PRISM
python3 drug_response_prediction.py \
    -i /mnt/output/write_sub/scanpyobj.sub.h5ad \
    -o /mnt/output/

# GDSC
python3 drug_response_prediction.py \
    -i /mnt/output/write_sub/scanpyobj.sub.h5ad \
    -o /mnt/output/ \
    -m GDSC
```