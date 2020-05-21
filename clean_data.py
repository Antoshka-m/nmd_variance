#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cleaning timeseries data:
flattening, filtering (low-pass), removing spikes if needed.
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os.path
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
from scipy.stats import linregress

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


def flattening(data, fps, chunk_l):
    """
    Make flattening of raw data using windows (time chunks) of given size
    
    Parameters
    -------
    
    data: list or np.array
        timeseries data (deflection magnitude)
    dt: float
        time resolution of data
    chunk_l: float
        length of the chunk (in seconds)
        
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

def plot_flattening(t, x, x_flat, x_fit):
    plt.figure(figsize=(8, 6))
    # plt.grid()
    plt.plot(t, x, label='raw timeseries data', color='blue')
    plt.plot(t[:len(x_fit)], x_fit, color='red')
    plt.plot(t[:len(x_flat)], x_flat, label='flattened', color='orange')
    plt.xlabel('t, s', fontsize=24)
    plt.ylabel('Displacement, a.u.', fontsize=24)
    plt.legend(loc='upper right')
    plt.tight_layout()
    # plt.show()

def plot_signal(t, x, title=None):
     plt.figure(figsize=(8, 6))
     plt.title(title)
     plt.grid()
     plt.plot(t, x)
     plt.xlabel('t, s', fontsize=24)
     plt.ylabel('Displacement, a.u.', fontsize=24)
# import data

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
    fps=df['fps'].iloc[0]
    y_norm_fit, y_norm_flat = flattening(df['Iz']/df['sum'], fps, chunk_l=60)
    x_norm_fit, x_norm_flat = flattening(df['Il']/df['sum'], fps, chunk_l=60)
    data_flat = pd.concat([data_flat,
                           pd.DataFrame({'Iz norm': y_norm_flat,
                                          'Il norm': x_norm_flat,
                                          'experiment': np.full(len(y_norm_flat), exp_labels[i])})],
                          axis=0)
    plot_flattening(df['t'], df['Iz']/df['sum'], y_norm_flat, y_norm_fit)
    plot_signal(df['t'][:len(y_norm_flat)], y_norm_flat, 'vertical deflection')
    plot_flattening(df['t'], df['Il']/df['sum'], x_norm_flat, x_norm_fit)
    plot_signal(df['t'][:len(x_norm_flat)], x_norm_flat, 'horizontal deflection')
    i+=1

data_flat =pd.DataFrame({'t': np.arange(0, len(data_flat)/fps, 1/fps),
                         'Iz norm': data_flat['Iz norm'],
                         'Il norm': data_flat['Il norm'],
                         'experiment': data_flat['experiment']})
plot_multiple_sigs('t', 'Iz norm', data=data_flat, title='flattened vertical deflection (over sum)')
plot_multiple_sigs('t', 'Il norm', data=data_flat, title='flattened horizontal deflection (over sum)')
plt.show()
# plot_signal(np.arange(0, len(data_flat)/fps, 1/fps), data_flat['Iz norm'])
# plot_signal(np.arange(0, len(data_flat)/fps, 1/fps), data_flat['Il norm'])