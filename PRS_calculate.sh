# ===============================================================
# Describe: automatic PRS calculation pipeline
# Author: Rogen
# Date: 2022.05.23
# Parameters: Base summary statistics, Target bfile, pheno file, 
# 			score normalize: {1:std | 0:avg(default)}, Output
# None: if the MAF command is specificity, the base 
# 		summary statistics file must contain MAF column
# ===============================================================
#! /usr/bin/bash
# source activate R3.5

Help()
{
	# Display Help
	echo
	echo "***Pipeline for PRS v1.1***"
	echo "Add description of the script functions here:"
	echo
	echo "  -a | --train       : PRS train prefix {it is specified only with --validate for run validate dataste}."
	echo "  -g | --gwas        : Summary statisitcs file of GWAS."
	echo "  -t | --target      : Target set that you wanna calculate PRS (bfile)."
	echo "  -v | --validate    : Validate set that you wanna calculate PRS (bfile)."
	echo "  -k | --keep        : Remain what sample you wanna analyza (Optional)."
	echo "  -r | --remove      : Delete what sample you don't wanna analyza (Optional)."
	echo "  -l | --regression  : Associated with logit or linear regression method [ logit {default} | linear ] (Optional)."
	echo "  -e | --extract     : Pick up what SNP you wanna analyze (Optional)."
	echo "  -x | --exclude     : Drop out what SNP you don't wanna analyze (Optional)."
	echo "  -p | --pheno       : Filepath of phenotype table (including covariants)."
	echo "  -z | --ignore_fid  : for PLINK2.0 phenotable (First column is IID, default: 0) (Optional)."
	echo "  -s | --score       : PRS calcuated manna [ avg (default) | sum | std ] (Optional)."
	echo "  -m | --trait       : Column name of trait in phenotype table."
	echo "  -i | --snp_id      : Column name of SNP of input summary statistics file {default: SNP} (Optional)."
	echo "  -c | --covar       : Colunm name of covariant in phenotype {Comma delimiter. eg: Sex,Age,@PC[1-2]} (Optional)."
	echo "  -f | --maf_filter  : Fiter variant with MAF {default: MAF:0.01}"
	echo "  -o | --out         : Output name (default: prsice as prifix of ouput name)."
	echo
}


re='^(-+help|-h)$';
if [[ $1 =~ $re ]]; then
	Help
else
	while [ "$#" -gt 0 ]; do
		case "$1" in
			-a|--train) train="$2"; shift 2;;
			-g|--gwas) gwas="$2"; shift 2;;
			-v|--validate) validate="$2"; shift 2;;
			-t|--target) target="$2"; shift 2;;
			-k|--keep) keep="$2"; shift 2;;
			-r|--remove) remove="$2"; shift 2;;
			-l|--regression) regression="$2"; shift 2;;
			-e|--extract) extract="$2"; shift 2;;
			-x|--exclude) exclude="$2"; shift 2;;
			-f|--maf_filter) maf_filter="$2"; shift 2;;
			-i|--snp_id) snp_id="$2"; shift 2;;
			-p|--pheno) pheno="$2"; shift 2;;
			-z|--ignore_fid) ignore_fid="$2"; shift 2;;
			-m|--trait) trait="$2"; shift 2;;
			-s|--score) score="$2"; shift 2;;
			-c|--covar) covar="$2"; shift 2;;
			-o|--out) out="$2"; shift 2;;
			*) echo "unknown option: $1" >&2; exit 1;;
			*) handle_argument "$1"; shift 1;;
		esac
	done
	

	train_cmd="bash ./PRS_for_target.sh \
			-g $gwas \
			-t $target \
			-p $pheno "

	val_cmd="bash ./PRS_for_validate.sh \
			-g $gwas \
			-v $validate \
			-p $pheno "


    # ====== Set association method ======
    if [ -z $regression ] || [ $regression == 'logit' ];
    then
        regression='logit'

        cmd+="-l logit "
    elif [ $regression == 'linear' ];
    then
        cmd+="-l $regression "
    fi

	# ====== Set keep and remove ======
	if [[ ! -z $keep ]] && [[ ! -z $remove ]]; then
		echo -e '--keep is mutually exclusive with --remove...'
		exit

	elif [[ ! -z $keep ]]; then
		cmd+="-k $keep "

	elif [[ ! -z $remove ]]; then
		cmd+="-r $remove "

	else
		cmd=$cmd
	fi

	# ====== Set extract and exclude ======
	if [[ ! -z $extract ]] && [[ ! -z $exclude ]]; then
		echo -e '--extract is mutually exclusive with --exclude...'
		exit

	elif [[ ! -z $extract ]]; then
		cmd+="-e $extract "

	elif [[ ! -z $exclude ]]; then
		cmd+="-x $exclude "

	else
		cmd=$cmd
	fi
	
	# ====== Set MAF filter ======
	if [[ -z $maf_filter ]];
	then
		cmd+="-f MAF:0.01 "
	else
		cmd+="-f $maf_filter "
	fi

	# ====== Set how type PRS score would be calculated ======
	if [ -z $score ];
	then
		score="avg"
	fi

	# ====== Set covariant ======
	if [[ ! -z $covar ]];
	then
		cmd+="-c $covar "
	fi

	# ====== Set SNP ======
	if [[ ! -z $snp_id ]];
	then
		cmd+="-i $snp_id "
	fi

	# ===== Ignore FID for PRSice2 =====
	if [[ ! -z $ignore_fid ]] && [[ $ignore_fid == 1 ]];
	then
		cmd+="-z 1 "
	fi

	# ====== Set output name ======
	if [ -z $out ];
	then
		out="prsice"
	fi

	cmd+="-s $score \
		-m $trait "

	train_cmd+=$cmd
	train_cmd+="-o $out.target"

	val_cmd+=$cmd
	val_cmd+="-a $out.target \
			-o $out.validate"


	# ====== Show progressing ======
	echo -e '\n\nStarting to calcuate PRS...'
	printf -- '=%.0s' {1..100}
	echo -e ""
	echo $train_cmd
	printf -- '=%.0s' {1..100}

	eval $train_cmd

	echo -e "\nTraining have done...\n\n"

	echo -e "Starting valdate..."
	printf -- '=%.0s' {1..100}
	echo -e ""
	echo $val_cmd
	printf -- '=%.0s' {1..100}

	eval $val_cmd

	echo -e "\nValdate have done...\n\n"
fi


if [ -f '.RData' ]; then
	rm .RData
fi


<< command
=========== Example ===========
./PRS_for_target.sh \
	-g /DATA6/Rogen/Auto_Run_PLINK_GWAS/Gout_Demo/Gout_TPMI_imputed_adjGWAS.Gout.glm.logistic \
	-t /DATA6/SNParray/30k_qc/impute/all_0214/pgen/bfile/CMUH_imputed_r2.MAF.HWE \
	-k /DATA6/Rogen/Auto_Run_PLINK_GWAS/Gout_Demo/Gout_Demo_Pheno.txt \
	-p /DATA6/Rogen/Auto_Run_PLINK_GWAS/Gout_Demo/Gout_Demo_Pheno.txt \
	-s std \
	-i ID \
	-c UA,Sex,Age \
	-o prs \
	-m Gout \
	-z 1 \
	-f MAF:0.01

./PRS_for_validate.sh \
	-a prs.train \
	-g /DATA6/Rogen/Auto_Run_PLINK_GWAS/Gout_Demo/Gout_TPMI_imputed_adjGWAS.Gout.glm.logistic \
	-v /DATA6/SNParray/30k_qc/impute/all_0214/pgen/bfile/CMUH_imputed_r2.MAF.HWE \
	-r /DATA6/Rogen/Auto_Run_PLINK_GWAS/Gout_Demo/Gout_Demo_Pheno.txt \
	-p /DATA6/Rogen/Auto_Run_PLINK_GWAS/Gout_Demo/Gout_Demo_Pheno.txt \
	-s std \
	-i ID \
	-c UA,Sex,Age \
	-z 1 \
	-o prs
command