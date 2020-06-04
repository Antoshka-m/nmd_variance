#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculation of variance in non-overlapping windows
"""


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os.path
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
from scipy.stats import linregress

do_flattening = True
flat_window = 60 #len of window for flattening in seconds
var_window = 10 # window for calculating variance
x_name = 'Il'
y_name = 'Iz'
# if normalize over sum, use the following:
x_norm_name = 'Il norm'
y_norm_name = 'Iz norm'

def get_filenames():
    """
    get filenames with graphical interface
        
    Returns
    -------
    file_list : list
        list of chosen filenames
    """
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    file_list=askopenfilenames()
    return file_list

def get_exp_labels(file_list):
    print('By default, the labels of experiments will be taken from filenames:\n')
    exp_labels = [os.path.basename(file)[:-4] for file in file_list]
    for i, exp in enumerate(exp_labels):
        print(str(i+1)+'. ' + exp)
    print('If you would like to you the default ones, just press enter.')
    print('Otherwise enter all the labels manually (dont use the same ones twice!)')
    label=input('Enter the label of the first file %s\n' %exp_labels[0])
    if label=='':
        return exp_labels
    else:
        exp_labels[0]=label
        for i in range(1, len(exp_labels)):
            label=input('Enter the labels of the file %s\n' %exp_labels[i])
            while (label in exp_labels[:i]):
                label=input('Label %s already exists! Enter the other label of the file %s\n' %(label, exp_labels[i]))
            exp_labels[i]=label
    return exp_labels

def flattening(data, fps, chunk_l):
    """
    Make flattening of raw data using windows (time chunks) of given size
    
    Parameters
    -------
    
    data: list or np.array
        timeseries data (deflection magnitude)
    fps: float
        frames/sec of recording (or acquisition frequency)
    chunk_l: float
        length of the chunk window (in seconds)
        
    Returns
    -------
    data_fit: np.array
        fitting line
    data_flat: np.array
        flattened data
    """
    # get chunk length in datapoints instead of seconds
    # chunk_l = int(chunk_l//dt )
    chunk_l = int(fps*chunk_l)
    # get size of output flattened data (cut out remaining part after division)
    flat_data_size = len(data)//chunk_l
    data_fit = np.empty((flat_data_size, chunk_l))
    for chunk_i in range(flat_data_size):
        start_i = chunk_i*chunk_l
        end_i = (chunk_i+1)*chunk_l
        y = data[start_i:end_i]
        # make linear regression for each chunk
        k_coef, b_coef, _, _, _ = linregress([i for i in range(len(y))], y)
        yfit = k_coef * np.array([i for i in range(len(y))]) + b_coef
        data_fit[chunk_i] = yfit
    data_flat = data[:(flat_data_size*chunk_l)]-data_fit.flatten()
    return data_fit.flatten(), data_flat


def var_calc(data, fps, chunk_l):
    """Calculate variance of data in chosen non-overlapping windows
    
    Parameters
    -------
    
    data: list or np.array
        timeseries data (deflection magnitude)
    fps: float
        frames/sec of recording (or acquisition frequency)
    chunk_l: float
        length of the chunk window (in seconds)
        
    Returns
    -------
    var_data: np.array
        variance
    """
    # get chunk length in datapoints instead of seconds
    chunk_l = int(fps*chunk_l)
    # get size of data which is taken for calculation (cut out remaining part after division)
    data_size = len(data)//chunk_l
    var_data = np.empty(data_size)
    for chunk_i in range(data_size):
        start_i = chunk_i*chunk_l
        end_i = (chunk_i+1)*chunk_l
        x = data[start_i:end_i]
        var_data[chunk_i] = np.var(x)
    return var_data

def plot_multiple_sigs(t, sig, data, title=None):
    plt.figure(figsize=(8, 6))
    sns.set()
    plt.title(title)
    fig=sns.lineplot(x=t, y=sig, data=data, hue='experiment', ci=None)
    fig.set_xlabel('t, s', fontsize=22)
    fig.set_ylabel(sig+', a.u.', fontsize=22)
    fig.tick_params(labelsize=18)
    plt.legend(loc='upper right')
    plt.setp(fig.get_legend().get_texts(), fontsize='14') # for legend text
    plt.setp(fig.get_legend().get_title(), fontsize='18') # for legend title
    plt.tight_layout()


filenames = get_filenames()
print('following files will be imported:')
for file in filenames:
    print(file)
exp_labels=get_exp_labels(filenames)
data = pd.DataFrame()
data_flat = pd.DataFrame()
var_data_flat = pd.DataFrame()
i=0
for file in filenames:
    print('Analyzing %s...' % os.path.basename(file))
    print('Importing file...')
    df=pd.read_csv(file)
    df=pd.concat([df, pd.DataFrame({'experiment': np.full(len(df), exp_labels[i])})],
                 axis=1)
    data=pd.concat([data,
                    df[['Iz', 'Il', 'sum', 'experiment']]],
                   axis=0,
                   ignore_index=True)
    if 'fps' in df.columns:
        fps=df['fps'].iloc[0]
    else:
        fps=input('Enter fps\n')
        fps=float(fps)
    y_fit, y_flat = flattening(df[y_name], fps, chunk_l=flat_window)
    x_fit, x_flat = flattening(df[x_name], fps, chunk_l=flat_window)
    # y_norm_fit, y_norm_flat = flattening(df[y_norm_name], fps, chunk_l=flat_window)
    # x_norm_fit, x_norm_flat = flattening(df[x_norm_name], fps, chunk_l=flat_window)
    var_x = var_calc(data=x_flat, fps=fps, chunk_l=var_window)
    var_y = var_calc(data=y_flat, fps=fps, chunk_l=var_window)
    print(var_x.mean())
    print(var_y.mean())
    # var_x_norm = var_calc(data=x_norm_flat, fps=fps, chunk_l=var_window)
    # var_y_norm = var_calc(data=y_norm_flat, fps=fps, chunk_l=var_window)
    data_flat = pd.concat([data_flat,
                           pd.DataFrame({'Iz': y_flat,
                                          'Il': x_flat,
                                          'experiment': np.full(len(y_flat), exp_labels[i])})],
                          axis=0)
    var_data_flat = pd.concat([var_data_flat,
                               pd.DataFrame({'var_x': var_x,
                                             'var_y': var_y,
                                             'experiment': np.full(len(var_y), exp_labels[i])})]).reset_index()
    i+=1
data_flat =pd.DataFrame({'t': np.arange(0, len(data_flat)/fps, 1/fps),
                         'Iz': data_flat['Iz'],
                         'Il': data_flat['Il'],
                         'experiment': data_flat['experiment']})
var_data_flat =pd.DataFrame({'t': np.arange(0, len(var_data_flat)*var_window, var_window),
                         'var_x': var_data_flat['var_x'],
                         'var_y': var_data_flat['var_y'],
                         'experiment': var_data_flat['experiment']}).reset_index()


plot_multiple_sigs('t', 'Iz', data=data_flat, title='flattened vertical deflection' )
plot_multiple_sigs('t', 'Il', data=data_flat, title='flattened horizontal deflection')
plot_multiple_sigs('t', 'var_x', data=var_data_flat, title='flattened vertical variance' )
plot_multiple_sigs('t', 'var_y', data=var_data_flat, title='flattened horizontal variance')



    
    