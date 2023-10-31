#!/usr/bin/env python
# coding: utf-8


# import data processing model
import numpy as np
import pandas as pd
from sklearn import preprocessing

# set pdf font fromat
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42



# import all expression data of all DEGs and label
expres = pd.read_csv("all_deg_gene_expres.csv",sep=",",index_col=0)
group_data = pd.read_table("group.csv",sep=",",index_col = 0)


# encode LRTIs as 1 and non-LRTI as 0
label_encoder = {"Noninfection":0,"Infection":1}
label_s = group_data.infection.apply(lambda x:label_encoder[x])



# split the discovery cohort to training and testing data with ratio of 7:3
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV

Xtrain, Xtest, Ytrain, Ytest = train_test_split(expres,label_s,test_size=0.3,
                                                random_state=0)


# feature selection using RFE
rfc_s = LR(penalty="l2",solver="liblinear",C=0.1, max_iter=300)

from sklearn.feature_selection import RFE

score = []
for i in range(5,20,1):
    x_wrapper = RFE(rfc_s,n_features_to_select=i, step=10).fit_transform(expres,label_s)
    once = cross_val_score(rfc_s,x_wrapper,label_s,cv=10).mean()
    score.append(once)
    
    print("shape of i={}:".format(i),x_wrapper.shape[1],
          "with score: ",once)



# when 14 features to selected, model showed best corss_validation 
selector_x= RFE(rfc_s,n_features_to_select=14, step=10).fit(expres,label_s)
x_wrapper = selector_x.transform(expres)


# check which 14 features were selected
expres.loc[:,status_x].columns

