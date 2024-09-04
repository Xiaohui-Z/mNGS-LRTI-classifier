# import model
import numpy as np
import pandas as pd
from sklearn.svm import SVC
import seaborn as sns
from sklearn import preprocessing

# set pdf style
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# import expression and label data
expres_m = pd.read_csv("DEG_expression.csv",sep=",",index_col=0)


# splite features and label
expres = expres_m.iloc[:,0:-1]
label_s = expres_m.label

#select best knearl
Kernel = ["linear","poly","rbf","sigmoid"]
for kernel in Kernel:
    clf= SVC(kernel = kernel
             , gamma="auto"
             ,degree=1
             , cache_size=5000
            )
    val = cross_val_score(clf
                        ,expres,label_s,cv=5)
    print("The mean accuracy under kernel %s is %f" % (kernel,val.mean()),val.std().round(4)*100)
    
    
# adjust C parameter
score = []
C_range = np.linspace(0.1,3,50) 
for i in C_range:
    clf = SVC(kernel="linear",C=i,cache_size=5000)
    val = cross_val_score(clf,x_wrapper,label_s,cv=10)
    
    score.append(val.mean())
    
print(max(score), C_range[score.index(max(score))])
plt.plot(C_range,score)
plt.show()

# get best parameter
clf = SVC(kernel="linear",C=1.8,cache_size=5000)
val = cross_val_score(clf,expres,label_s,cv=10)
val.mean(),val.std()

#  get the correlation coefficients for each variable in the mode
trained_clf = clf.fit(x_wrapper,label_s)
gene_coef = trained_clf.coef_
gene_coef = pd.DataFrame(gene_coef).T
gene_coef.loc[:,"ens"] = expres.loc[:,status_x].columns

merged_df = pd.merge(genename,gene_coef, on='ens', how='inner')

# heatmap the feature correlation
df = merged_df.iloc[:,2:7]
df.set_index('gene', inplace=True)


sns.clustermap(df, cmap="coolwarm", standard_scale=1, row_cluster=True, 
               col_cluster=False, annot=None, figsize=(5, 12), 
               dendrogram_ratio=(.2, .05), cbar_pos=(-0.1, 0.2, 0.03, 0.4))