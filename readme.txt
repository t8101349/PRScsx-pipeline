note: no need to copy scripts to target folder
start from PRSCSX_for_target.sh
example:
PRSCSX_for_target.sh  -b /SNParray/SourceShare/20240321_50w_Imputation/step12-pgen2bed/Axiom_imputed_r2.MAF  -g finaloutput_fixed_input.txt  -v v2  -n 654461 -c all  -o meta_prs/oralcancer_fixed_META_PRS  --target_list /home/Weber/Cancer/oralcancer/20250709/oralcancer_split.ind_target_set.tsv --val_list /home/Weber/Cancer/oralcancer/20250709/oralcancer_split.ind_validate_set.tsv

PRSCS_sumstat_adjust.py: output sumstats.prscs.txt
PRS_Distribution_Plot_csx.py: plot graph(for csx)
calc_prs_metrics.py: paint ROC curve, Precision-Recall curve