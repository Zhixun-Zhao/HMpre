#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 07/12/17 16:13 PM
# @Author  : ZHIXUN ZHAO
# @File    : HMpre.py
# @Software: PyCharm Community Edition

import xgboost as xgb
import numpy as np
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import math
from sklearn.externals import joblib
import matplotlib.pyplot as plt


def performace_est(test_label, predict_label):
    tn = 0
    fn = 0
    tp = 0
    fp = 0
    label = np.column_stack((test_label, predict_label))
    for i in range(len(label)):
        if label[i][0] == 0:
            if label[i][1] == 0:
                tn = tn + 1
            else:
                fp = fp + 1
        else:
            if label[i][1] == 0:
                fn = fn + 1
            else:
                tp = tp + 1
    Sn = tp / (tp + fn)
    Sp = tn / (tn + fp)
    ACC = (tp + tn) / (tp + tn + fp + fn)
    MCC = ((tp * tn) - (fp * fn)) / math.sqrt((tp + fn) * (tp + fp) * (tn + fn) * (tn + fp))
    Precision = tp / (tp + fp)
    Recall = tp / (tp + fn)
    F1 = Precision * Recall * 2 / (Precision + Recall)
    return Sn, Sp, ACC, MCC, Precision, Recall, F1


def results(true_label, predictscore, h):
    test_lab = true_label
    y_hat = predictscore

    fpr, tpr, thresholds = roc_curve(test_lab, y_hat)
    roc_auc = auc(fpr, tpr)
    roc_auc = auc(fpr, tpr)
    precision, recall, threshold = precision_recall_curve(test_lab, y_hat)
    auprc = auc(recall, precision)
    y_hat[y_hat > h] = 1
    y_hat[~(y_hat > h)] = 0
    Sn, Sp, ACC, MCC, Precision, Recall, F1 = performace_est(test_lab, y_hat)
    return Sn, Sp, ACC, MCC, Precision, Recall, F1, roc_auc, auprc


train_data = np.loadtxt('train_data.txt', dtype=float)
train_label = np.loadtxt('train_label.txt', dtype=float)
test_data = np.loadtxt('test_data.txt', dtype=float)
test_label = np.loadtxt('test_label.txt')

data_train = xgb.DMatrix(train_data, label=train_label)

param = {'lambda': 700, 'max_depth': 6, 'eta': 0.1, 'silent': 1, 'objective': 'binary:logistic',
         'booster': 'gbtree','scale_pos_weight': 6, 'eval_metric': 'auc'}

bst = xgb.train(param, data_train, num_boost_round=400)
pred_prob = bst.predict(test_data)

Sn, Sp, ACC, MCC, Precision, Recall, F1, roc_auc, auprc = results(test_label, pred_prob, 0.5)
print('HMpre prediction results')
print('Sn=', Sn, 'Sp=', Sp, 'Precision=', Precision, 'Recall=', Recall, 'F1=', F1, 'ACC=', ACC, 'MCC=', MCC, 'ROC=',
        roc_auc, 'AUC=', auprc)
