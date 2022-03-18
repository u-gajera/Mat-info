#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import re


# In[2]:


print('read the flow using formula_reader.standard_flow()')


# In[3]:


def standard_flow():
    
    print("------Calling and cleaning the datasets------")
    print("dataset_material_Fenad = pd.read_excel('FENAD.xlsx', 'list_higher_Pauling_electronega')")
    print("dataset_atomic_Fenad  = pd.read_excel('FENAD.xlsx', 'NAD AtomicData')")
    print("dataset_atomic_Fenad  = dataset_atomic_Fenad.apply(pd.to_numeric,errors='ignore')")
    print("dataset_material_Fenad['B']=dataset_material_Fenad['B'].str.strip()")
    print("dataset_material_Fenad['C']=dataset_material_Fenad['C'].str.strip()")
    print("\n")
    print("------Cleaning the dataset from inside------")
    print("cleaning_the_dataset(dataset_atomic_Fenad)")
    print("cleaning_the_dataset(dataset_material_Fenad)")
    print("\n")
    print("------Assinging the column name------")
    print("column_name = columns_name(dataset_atomic_columns=dataset_atomic_Fenad.keys()[1:], num_of_element=['A','B','C'])")
    print("\n")
    print("------Combining the datasets------")
    print("dataset_test=combining_the_dataframe(dataset_material=dataset_material_Fenad, dataset_atomic=dataset_atomic_Fenad, \n                                 column_name_list=column_name)")
    
    print("\n")
    print("------Final equation reading------")
    print("equation = '(Mr_B - sqrt(fabs(Mr_C)))/(LUMO_B**2 + exp(HOMO_A))")
    print("object = equation_reader(equation_string=equation,dataset=dataset_test)")
    print("feature_array = np.array(check.run_equation_reader())")
    print("feature_array = np.reshape(feature_array,(-1,1))")
    print("a,b,c,d = check.get_individual_element()")
    print("formula = check.get_formula()")


# In[4]:


def cleaning_the_dataset(dataset_material):
    for ii in dataset_material.keys():
        if (ii.find('Unnamed') != -1):
            del dataset_material[ii]
        elif list(set(pd.isnull(dataset_material[ii].values)))[0]==True:
            del dataset_material[ii]
        else:
            pass


# In[5]:


def columns_name(dataset_atomic_columns, num_of_element=['A','B']):
    column_name = []
    for ii in num_of_element:
        for jj in dataset_atomic_columns:
            column_name.append(jj+'_'+ii)
    return column_name


# In[6]:


def combining_the_dataframe(dataset_material, dataset_atomic, column_name_list):
    df_test = pd.DataFrame(columns=column_name_list)
    for ii,value in enumerate(dataset_material.iterrows()):
        df_test.loc[ii]=(list(dataset_atomic[dataset_atomic['A']==value[1][1]].values[0][1:]) +                         list(dataset_atomic[dataset_atomic['A']==value[1][2]].values[0][1:]) +                         list(dataset_atomic[dataset_atomic['A']==value[1][3]].values[0][1:])  )
    dataset_material = pd.concat([dataset_material, df_test], axis=1)
    return dataset_material


# In[7]:


class equation_reader:
    
    def __init__(self, equation_string, dataset):
        self.equation_string=equation_string
        self.dataset        = dataset
        self.operation_list = [1,1,1,1]
        self.sign_list      = ['+','+','+','+']
        
    def first_split(self):
        equation = self.equation_string.split('/')
        temp_list = []
        for ii,value in enumerate(equation):
            if '+' in equation[ii]:
                nom1,nom2 = equation[ii].split('+')
                temp_list.append(nom1)
                temp_list.append(nom2)
            elif '-' in equation[ii]:
                self.sign_list[1] = '-'
                nom1,nom2 = equation[ii].split('-')
                temp_list.append(nom1)
                temp_list.append(nom2)
            else:
                temp_list.append(equation[ii])
        return temp_list
            
    def clear_temp_list(self):
        temp_list = self.first_split()
        temp_list2 = []
        characters_to_remove = ["sqrt","fabs","exp"]
        for jj,value in enumerate(temp_list):
            result = re.sub('[*()]', ' ', value).lstrip().rstrip()
            new_string = result
            
            if 'sqrt' in new_string:
                self.operation_list[jj] = 'sqrt'
            elif 'exp' in new_string:
                self.operation_list[jj] = 'exp'
            for character in characters_to_remove:
                new_string = new_string.replace(character, "")
            if len(re.findall(r'\d+',result)) !=0:
                self.operation_list[jj] = int(re.findall(r'\d+',result)[0])
                new_string = re.sub(r'[0-9]+', '', new_string)

            temp_list2.append(new_string.lstrip().rstrip())
        return temp_list2
    
    def collecting_information_from_dataset(self):
        main_element_list = self.clear_temp_list()
        feature_array_temp= []
        for ii in main_element_list:
            try:
                feature_array_temp.append(self.dataset[ii])
            except:
                ValueError
        return feature_array_temp
    
    def operation_with_operators(self):
        feature_array = []
        feature_array_temp = self.collecting_information_from_dataset()
    
        for kk,value in enumerate(self.operation_list):
            try:
                if type(value) == int:
                    feature_array.append(feature_array_temp[kk]**value)
                elif value == 'exp':
                    feature_array.append(np.exp(feature_array_temp[kk]))
                elif value == 'sqrt':
                    feature_array.append(np.sqrt(np.abs(feature_array_temp[kk])))
                else:
                    print ('Operator are not the number exp or sqrt')
            except:
                pass
        return feature_array
    
    def operation_with_sign(self):
        feature_array= self.operation_with_operators()
        for ll, value in enumerate(self.sign_list):
            if value == '-':
                feature_array[ll] = feature_array[ll]*(-1) 
            else:
                pass
        return feature_array
    
    def final_combined(self):
        feature_array = self.operation_with_sign()
        try:
            return (feature_array[0]+feature_array[1])/(feature_array[2]+feature_array[3])
        except:
            return (feature_array[0]+feature_array[1])/(feature_array[2])

    def run_equation_reader(self):
        fnct= self.final_combined()
        print('origional eqution: ',self.equation_string)
        print('operation list   : ',self.operation_list)
        print('sign list        : ',self.sign_list)
        print('feature list     : ',self.clear_temp_list())
        return fnct
    
    def get_individual_element(self):
        return self.operation_with_operators()
    
    def get_formula(self):
        return self.equation_string

