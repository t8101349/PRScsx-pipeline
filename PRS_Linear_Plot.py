# =======================================
# Describe: PRS plot in linear trait
# Author: Rogen 
# Build_Date: 2024.05.29
# Last_update: 2024.06.21
# Update: 
# 2024.09.05 - fix CV method into KFold because StratifiedKFold can't split continous target
# 2024.06.14 - change pearson corr to pearson square
# 2024.06.21 - revise few bugs
# =======================================
#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
from scipy import stats
from sklearn.linear_model import Lars, LarsCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import r2_score
from argparse import ArgumentParser


# custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks")

def argparser():
    parser = ArgumentParser()
    parser.add_argument('--phenotable','-i',type=str,required=True,help='Input pheno file (with tab separate).')
    parser.add_argument('--fid_absent','-f',type=bool,required=False,default=False,help='fid column is absent or not.')
    parser.add_argument('--ignore_nan_pheno','-g',type=bool,required=False,default=False,help='Ignore phenotype which is none or -9.')
    parser.add_argument('--phename','-m',type=str,required=True,help='Input pheno colunm name.')
    parser.add_argument('--prs','-p',type=str,required=False,default='PRS',help='Input PRS column name.')
    parser.add_argument('--ratio','-r',type=float,required=False,default=0.8,help='How much ratio separate as train dataset (default: 0.8).')
    parser.add_argument('--metrics','-s',type=str,required=False,default='pearsonr',help='linear correlation score [ rscore | pearsonr(default) ] .')
    parser.add_argument('--covname','-n',type=str,required=False,default='',help='Types what column names as covariance (eg. Sex,Age,PC1-10).')
    parser.add_argument('--out','-o',type=str,required=False,default='linear_plot',help='Output name (default: linear_plot).')
    
    args = parser.parse_args()

    return(args)


def sample_split(df,ratio,covname,seed=1):
    base = pd.DataFrame()
    target = pd.DataFrame()
    df_g = df.groupby(by=covname)
    condictioner = list(df_g.groups.keys())
    
    for c in condictioner:
        sub_df = df_g.get_group(c)

        if sub_df.shape[0] < 2:
            base = base.append(sub_df)
        else:
            sub_df_base, sub_df_target = train_test_split(sub_df,stratify=sub_df[covname],shuffle=True,train_size=ratio,random_state=seed)
            base = base.append(sub_df_base)
            target = target.append(sub_df_target)
            
    return(base,target)


def plot(data,metrics,name='relplot'):
    g = sns.relplot(data=data,y='Real',x='Predict',hue='Diff>0.1',s=20, alpha=0.8,linewidth=0)
    plt.axline((data.Predict.min(), data.Predict.min()), (data.Predict.max(), data.Predict.max()), linewidth=1, linestyle='--', color='black')
  
    ax = g.axes[0,0]
    y_text_text_loc = int(data.Real.min()+(data.Real.max()-data.Real.min())/2)
    x_axis_text_loc = int(data.Predict.min()+(data.Predict.max()-data.Predict.min())/2)
#     ax.plot([data.Predict.min(), data.Predict.min()], [data.Predict.max(), data.Predict.max()], linewidth=1, linestyle='--', color='black')
    
    if metrics == 'pearsonr':
        val_corr, val_pvale = stats.pearsonr(data.Real,data.Predict)
        ax.text(x_axis_text_loc,y_text_text_loc,'Pearson\'s square: %.4f, p.val: %.4f' % (round(np.square(val_corr),4),round(val_pvale,4)))
    elif metrics == 'rscore':
        r2_test = r2_score(data.Real, data.Predict)
        ax.text(x_axis_text_loc,y_text_text_loc,'R-square: %.4f' % round(r2_test,4))

    g.fig.set_figwidth(8)
    g.fig.set_figheight(5)
    g.fig.savefig("%s.png"%name,dpi=400)

    data.to_csv('%s_predval.tsv'%name,sep='\t',index=False)


if __name__ == '__main__':
    args = argparser()
    print('-'*110)
    print(args)
    print('-'*110)

    phenotable = args.phenotable
    phenoname = args.phename
    prs = args.prs
    ratio = args.ratio
    covname = args.covname
    out = args.out
    metrics = args.metrics
    fid_absent = args.fid_absent
    ignore_nan_pheno = args.ignore_nan_pheno


    df = pd.read_table(phenotable,low_memory=False)
    # df = pd.read_excel(phenotable)

    if ignore_nan_pheno:
        df = df[(~df[phenoname].isna()) & (df[phenoname] != -9)]

    if fid_absent and 'FID' in df.columns:
        df.drop(['FID'],axis=1,inplace=True)


    train, test = sample_split(df,ratio,['Sex','Age'],seed=0)

    sp_covname = covname.split(',')
    for i in sp_covname:
        match = re.match(r'PC(\d+)\-(\d+)',i)
        if match:
            number = match.groups()
            PCs = list(map(lambda x: 'PC'+str(x) ,range(int(number[0]),int(number[1])+1)))
            sp_covname.remove(i)
            sp_covname = sp_covname + PCs

    x_train = train[[prs]+sp_covname]
    x_test = test[[prs]+sp_covname]
    y_train = train[phenoname]
    y_test = test[phenoname]

    # print(x_train.shape, y_train.shape, x_test.shape, y_test.shape)
    feature = x_train.columns
    
    scaler = StandardScaler()
    x_train = scaler.fit_transform(x_train)
    x_test = scaler.fit_transform(x_test)

    lars = LarsCV(verbose=True,cv=KFold(n_splits=10,shuffle=True,random_state=2),n_jobs=-1).fit(x_train, y_train)

    cept = lars.intercept_
    fi = feature[np.argsort(abs(lars.coef_))[::-1]].tolist()
    fs = sorted(abs(lars.coef_),reverse=True)
    fi = ['Intercept'] + fi
    fs = [cept] + fs
    fs_df = pd.DataFrame([fi,fs]).T
    fs_df.columns = ['Feature','Score']
    fs_df.to_csv('{}_feature_importance.tsv'.format(out),sep='\t',index=False)
    
    print('Feature importance:')
    print(fs_df)
    
    self_predict = lars.predict(x_train)
    val_predict = lars.predict(x_test)

    train_diff_01 = (abs(y_train - self_predict) / y_train) > 0.1
    train_data = pd.DataFrame([y_train.tolist(),self_predict,train_diff_01],index=['Real','Predict','Diff>0.1']).T

    test_diff_01 = (abs(y_test - val_predict) / y_test) > 0.1
    test_data = pd.DataFrame([y_test.tolist(),val_predict,test_diff_01],index=['Real','Predict','Diff>0.1']).T

    train_data = train_data.astype(np.float_)
    train_data['Diff>0.1'] = train_data['Diff>0.1'].astype(np.bool_)
    train_data['IID'] = train.IID.astype(np.str_)

    test_data = test_data.astype(np.float_)
    test_data['Diff>0.1'] = test_data['Diff>0.1'].astype(np.bool_)
    test_data['IID'] = test.IID.astype(np.str_)

    plot(train_data,metrics=metrics,name='{}_train'.format(out))
    plot(test_data,metrics=metrics,name='{}_test'.format(out))

    val_corr, val_pvalue = stats.pearsonr(df[phenoname],df[prs])
    g = sns.relplot(data=df,x=df[prs],y=df[phenoname],s=20,alpha=0.8,linewidth=0)
    g.fig.set_figwidth(8)
    g.fig.set_figheight(5)
    ax = g.axes[0,0]
    ax.text(0,df[phenoname].max()/2,'Pearson\'s square: %.4f' % round(np.square(val_corr),4))
    g.fig.savefig("{}_realdata_vs_prs.png".format(out),dpi=400)
