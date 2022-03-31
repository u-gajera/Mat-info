#!/usr/bin/env python
# coding: utf-8
# Author: Uday Gajera


import numpy as np
import pandas as pd
import re
from collections import Counter
import matplotlib.pylab as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import os
import formula_reader

def cleaning_the_dataset(dataset_material):
    for ii in dataset_material.keys():
        if (ii.find('Unnamed') != -1):
            del dataset_material[ii]
        elif list(set(pd.isnull(dataset_material[ii].values)))[0]==True:
            del dataset_material[ii]
        else:
            pass


dataset_material = pd.read_excel('Data_Ghiringhelli.xlsx', 'Material Data')
dataset_atomic = pd.read_excel('Data_Ghiringhelli.xlsx',   'Atomic Data')
dataset_atomic.columns = ['Z', 'A', 'IP', 'EA', 'HOMO', 'LUMO', 'rs', 'rp', 'rd', '1D']
dataset_test = pd.read_csv('dataset_test.csv')
dataset_test_new = pd.read_csv('dataset_test_new.csv')


cleaning_the_dataset(dataset_atomic)
cleaning_the_dataset(dataset_material)

column_name = formula_reader.columns_name(dataset_atomic_columns=dataset_atomic.keys()[2:-1], 
                           num_of_element=['A','B'])

equation_list = ['(rp_B**3 - exp(rs_B))/((rp_A**2))', 
                 '(sqrt(fabs(EA_A))+sqrt(fabs(IP_B)))/((rp_A**2))']


df_form = pd.read_excel('Final wrap up.xlsx', sheet_name='ZincBlend 1D')

all_num,all_equations,all_avg_mse = [],[],[]
for str_lst in df_form['OAD_gen_1 RMSE']:
    try :
        lindex = str_lst.index('(')
        rindex = str_lst.rindex(')')
        equation=str_lst[lindex:rindex+1]
        all_equations.append(equation)
    except:
        pass

new_equation_list = []
old_rmse, old_r2, old_cls   =[],[],[]
new_rmse, new_r2, new_cls   =[],[],[]
coeff_list                  =[]
min_slops, min_intercepts   =[],[]
r2_scores, cls_acss         =[],[]
regression = 'On'
rounding   = 12
target_property = 'r2_score' # 'RMSE' or 'classification_accuracy'

actute_species = ['RS' if i < 0 else 'ZB' for i in dataset_test.DE]
def target_func(pred, actu):
    rmse = np.round(np.sqrt(mean_squared_error(pred, actu)),rounding)
    r2 = r2_score(pred, dataset_test.DE)
    pred_species = ['RS' if i < 0 else 'ZB' for i in pred]
    quant_res    = ([i for i, j in zip(actute_species, pred_species) if i == j])
    cls_acs      = (len(quant_res)*100/82)
    return rmse, r2, cls_acs
        
for equation in all_equations[:10]:#+all_equations[4:7]+all_equations[8:9]:
    slops, intercepts = [],[]
    equation = equation.replace('HOMOKS','HOMO')
    equation = equation.replace('LUMOKS','LUMO')
    eq_obj = formula_reader.equation_reader(equation_string=equation, 
                                            dataset=dataset_test)
    feature_array = np.array(eq_obj.run_equation_reader())
    feature_array = np.reshape(feature_array,(-1,1))
    a,b,c = eq_obj.get_individual_element()
    formula = eq_obj.get_formula()
    new_equation_list.append(formula)
    actu = dataset_test['DE'].values

    regressor = LinearRegression()
    regressor.fit(feature_array, actu)
    
    pred = regressor.predict(np.array(feature_array))
    rmse, r2, cls_acs = target_func(pred=pred, actu=actu)
    old_rmse.append(target_func(pred=pred, actu=actu)[0])
    old_r2.append(target_func(pred=pred, actu=actu)[1])
    old_cls.append(target_func(pred=pred, actu=actu)[2])
    print('rmse: ',old_rmse[-1], 'r2 score: ', old_r2[-1], 'class accuracy: ',old_cls[-1])
    coeff = regressor.coef_
    intercept = regressor.intercept_
    
    temp_rmse, value_lst, temp_r2, temp_cls_acs = [],[],[],[]
    for ii in np.arange(-1.1,1.1,0.01):
        for jj in np.arange(-1.1,1.1,0.01):
                if np.round(ii,4)==0. or np.round(jj,4)==0.:
                    pass
                else:
                    x = (np.round(ii,4)*a + np.round(jj,4)*b)/(c)
                    if regression == 'On':
                        regressor = LinearRegression()
                        regressor.fit(np.array(x).reshape(-1,1), dataset_test.DE.values)
                        pred = regressor.predict(np.array(x).reshape(-1,1))
                        rmse, r2, cls_acs = target_func(pred=pred, actu=actu)
                        temp_rmse.append(rmse)
                        temp_r2.append(r2)
                        temp_cls_acs.append(cls_acs)
                        value_lst.append([np.round(ii,4), np.round(jj,4)])
                        slops.append(regressor.coef_)
                        intercepts.append(regressor.intercept_)
                    else:
                        pred = (coeff*x)+intercept
                        rmse, r2, cls_acs = target_func(pred=pred, actu=actu)
                        temp_rmse.append(rmse)
                        temp_r2.append(r2)
                        temp_cls_acs.append(cls_acs)
                        value_lst.append([np.round(ii,4), np.round(jj,4)])
    
    if target_property == 'RMSE':
        args = np.argmin(temp_rmse)
    
    elif target_property == 'classification_accuracy':
        args = np.argmax(temp_cls_acs)
        
    elif target_property == 'r2_score':
        args = np.argmax(temp_r2)
    
    print(value_lst[args], temp_rmse[args], temp_r2[args], temp_cls_acs[args])
    diff=  [0.1,0.01,0.001,0.0001]
    interval = [0.01,0.001,0.0001,0.00001]
    for ll in range(len(diff)):
        range_lst = value_lst[args]
        for ii in (np.arange(range_lst[0]-diff[ll],
                             range_lst[0]+diff[ll],interval[ll])):
            for jj in (np.arange(range_lst[1]-diff[ll], 
                                 range_lst[1]+diff[ll],interval[ll])):
                    if np.round(ii,4)==0. or np.round(jj,4)==0.:
                        pass
                    else:
                        x = (np.round(ii,4)*a + np.round(jj,4)*b)/(c)
                        if regression == 'On':
                            regressor = LinearRegression()
                            regressor.fit(np.array(x).reshape(-1,1), dataset_test.DE.values)
                            pred = regressor.predict(np.array(x).reshape(-1,1))
                            rmse, r2, cls_acs = target_func(pred=pred, actu=actu)
                            temp_rmse.append(rmse)
                            temp_r2.append(r2)
                            temp_cls_acs.append(cls_acs)
                            value_lst.append([np.round(ii,4), np.round(jj,4)])
                            slops.append(regressor.coef_)
                            intercepts.append(regressor.intercept_)
                        else:                            
                            pred = (coeff*x)+intercept
                            rmse, r2, cls_acs = target_func(pred=pred, actu=actu)
                            temp_rmse.append(rmse)
                            temp_r2.append(r2)
                            temp_cls_acs.append(cls_acs)
                            value_lst.append([np.round(ii,4), np.round(jj,4)])
        print(value_lst[args], temp_rmse[args], temp_r2[args], temp_cls_acs[args])    
    
    print('------------------------------------',len(temp_rmse))
    
    new_rmse.append(temp_rmse[args])
    new_r2.append(temp_r2[args])
    new_cls.append(temp_cls_acs[args])
    coeff_list.append(value_lst[args])
    min_slops.append(slops[args])
    min_intercepts.append(intercepts[args])



df = pd.DataFrame()
df['Equations'] = new_equation_list
df['old_rmse_for_full_dataset']=old_rmse
df['new_rmse']=new_rmse
df['old_r2'] = old_r2
df['new_r2'] = new_r2
df['old_class_accuracy']= old_cls
df['new_class_accuracy']= new_cls
df['individual_coeff'] = coeff_list
df['slops'] = [round(num[0], 3) for num in min_slops]
df['intecepts'] = min_intercepts

df.to_excel('optimization_using_r2_score.xlsx')

test_func = 0.633*(-0.85*dataset_test['rp_B'].values**2+1.09*np.exp(dataset_test['rs_B'].values))/dataset_test['rp_A'].values**2  -0.341319

pred_species = ['RS' if i < 0 else 'ZB' for i in test_func]
quant_res    = ([i for i, j in zip(actute_species, pred_species) if i == j])




