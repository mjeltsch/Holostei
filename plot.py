#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
import json
import pandas as pd
import matplotlib.pyplot as plt

def zero_div(x, y):
    if x == None:
        x = 'NaN'
    if y == None:
        y = 'NaN'
    try:
        return float(x) / float(y)
    except ZeroDivisionError:
        return 'NaN'

with open('holostei_vegfcd_mRNA.nxs.SLAC.json') as f:
    data = json.load(f)
    # Iterating through the json list
    #for i in data['MLE']:
    #    print(i)
    dN_div_dS_list = []
    p_value_list = []
    for i in data['MLE']['content']['0']['by-site']['AVERAGED']:
        print(zero_div(i[6], i[5]))
        dN_div_dS_list.append(zero_div(i[6], i[5]))
        if i[9] > 0.1:
            p_value_list.append(3)
        elif i[9] > 0.05: 
            p_value_list.append(2)
        elif i[9] > 0.01:
            p_value_list.append(1)
        else:
            p_value_list.append(0)

#create DataFrame
df = pd.DataFrame({'x': range(1,len(dN_div_dS_list)+1),
                   'y': dN_div_dS_list,
                   'z': p_value_list})

# Create scatterplot
plt.rcParams["figure.figsize"] = [15.00, 7.0]
plt.rcParams["figure.autolayout"] = True
plt.scatter(df.x, df.y, s=25, c=df.z, cmap='plasma', marker ='x')
plt.savefig("myimg.svg")

#fig, ax = plt.subplots(figsize = (15, 7))
#ax.scatter(df.x, df.y, s=25, c=df.z, cmap='copper')
#ax.set_yscale("log")
#plt.savefig("myimg.svg")