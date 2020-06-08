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
from scipy.stats import median_absolute_deviation
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

def get_export_dir(file):
    """
    create export directory (if it doesn't exist yet) in the same folder as import folder
        
    Returns
    -------
    exp_dir_path_list : list
        list of path of each export folder
    """
    # exp_dir_name = "export_files"
    exp_dir_name ='var_plots'
    exp_dir_path = os.path.join(os.path.dirname(file), exp_dir_name) # export directory for everything
    #exp_dir_path_all_types = []
    plots_dir_list = ['timeser_plots', 'barplots', 'boxplots']
    exp_dir_path_list = [] # path for folder with particular plots
    for plots_dir in plots_dir_list:
        if not os.path.isdir(os.path.join(exp_dir_path, plots_dir)): #create directory if it doesnt exist yet
            os.makedirs(os.path.join(exp_dir_path, plots_dir))
            print("\nExport directory %s was created. %s plots will be saved there" % (os.path.join(exp_dir_path, plots_dir), plots_dir))
        else:
            print("\nExport directory %s already exists. %s plots will be saved there" % (os.path.join(exp_dir_path, plots_dir), plots_dir))
        exp_dir_path_list.append(os.path.join(exp_dir_path, plots_dir))
    return exp_dir_path_list


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
    plt.figure(figsize=(8, 6), num=title)
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
    if title is not None:
        f_name = title.replace(' ', '_')
    plt.savefig(fname=os.path.join(export_dirs[0],
                             f_name + '.png'))
    
    
def plot_var_bars_mean(data, var_label, title=None):
    plt.figure(num=title)
    plt.title(title)
    sns.set()
    sns.barplot(x='experiment', y=var_label, data=data, ci='sd')
    if title is not None:
        f_name = title.replace(' ', '_')
    plt.savefig(fname=os.path.join(export_dirs[1],
                             f_name + '.png'))

def plot_var_bars_med(data, var_label, exp_labels, title=None):
    mad=np.empty(len(exp_labels))
    for i, exp in enumerate(exp_labels):
        y = data[data['experiment']==exp][var_label]
        mad[i]=median_absolute_deviation(y)
    plt.figure(num=title)
    plt.title(title)
    sns.set()
    sns.barplot(x='experiment', y=var_label, data=data, ci=None, estimator=np.median, yerr=mad)
    if title is not None:
        f_name = title.replace(' ', '_')
    plt.savefig(fname=os.path.join(export_dirs[1],
                             f_name + '.png'))
        

def plot_var_boxplot(data, var_label, title=None):
    plt.figure(num=title)
    plt.title(title)
    sns.set()
    sns.boxplot(x='experiment', y=var_label, data=data, showfliers=False)
    if title is not None:
        f_name = title.replace(' ', '_')
    plt.savefig(fname=os.path.join(export_dirs[2],
                             f_name + '.png'))
    
def rem_spikes_one_mean(x, n=3):
    mean = np.mean(x)
    std = np.std(x)
    x_clean = np.copy(x)
    # df = pd.DataFrame({'x': x, 't': t})
    # x_clean = df[abs(df['x'])>abs(mean+2*std)]['x'].apply(
            # lambda x: np.random.uniform(mean-2*std, mean+2*std))
    # mask = (df['x'])
    # x_clean = df['x'].apply(lambda x: 'np.random.uniform(mean-2*std, mean+2*std)' if abs(x)>abs(mean+2*std))
    # x_clean[abs(x_clean)>abs(mean+n*std)]=np.random.uniform(mean-n*std, mean+n*std, len(x_clean[abs(x_clean)>abs(mean+n*std)]))
    x_clean[(x_clean>mean+n*std)|(x_clean<mean-n*std)]=np.random.uniform(
        mean-n*std,
        mean+n*std,
        len(x_clean[(x_clean>mean+n*std)|(x_clean<mean-n*std)])) 
    return x_clean, (mean-n*std, mean+n*std)


def rem_spikes_one_median(x, n=3):
    med = np.median(x)
    mad = median_absolute_deviation(x)
    x_clean = np.copy(x)
    # df = pd.DataFrame({'x': x, 't': t})
    # x_clean = df[abs(df['x'])>abs(mean+2*std)]['x'].apply(
            # lambda x: np.random.uniform(mean-2*std, mean+2*std))
    # mask = (df['x'])
    # x_clean = df['x'].apply(lambda x: 'np.random.uniform(mean-2*std, mean+2*std)' if abs(x)>abs(mean+2*std))
    # x_clean[abs(x_clean)>abs(mean+n*std)]=np.random.uniform(med-n*mad, med+n*mad, len(x_clean[abs(x_clean)>abs(med+n*mad)]))
    x_clean[(x_clean>med+n*mad)|(x_clean<med-n*mad)]=np.random.uniform(
        med-n*mad,
        med+n*mad,
        len(x_clean[(x_clean>med+n*mad)|(x_clean<med-n*mad)]))
    return x_clean, (med-n*mad, med+n*mad)


def plot_flattening(t, x, x_flat, x_fit, ch_l, title=None):
    plt.figure(figsize=(8, 6))
    plt.title(title)
    # plt.grid()
    plt.plot(t, x, label='raw timeseries data', color='blue')
    for ch_i in range(len(x_flat)//ch_l):
        plt.plot(t[ch_i*ch_l:(ch_i+1)*ch_l], x_fit[ch_i*ch_l:(ch_i+1)*ch_l], color='red')
    #plt.plot(t[:len(x_fit)], x_fit, color='red')
    plt.plot(t[:len(x_flat)], x_flat, label='flattened', color='orange')
    plt.xlabel('t, s', fontsize=24)
    plt.ylabel('Displacement, a.u.', fontsize=24)
    plt.legend(loc='upper right')
    plt.tight_layout()
    # plt.show()
    

def calc_plot_var(data,
                  exp_labels,
                  filenames,
                  do_flat=True,
                  rem_spkies_median=False,
                  rem_spikes_mean=False,):
    """Doing calculation of the variance and plotting all necessary plots
    for all possible cases and types of corrections, that can be switched off
    or on using the boolean input parameters of the function"""
    
    
    
    
    
filenames = get_filenames()
print('following files will be imported:')
for file in filenames:
    print(file)
export_dirs=get_export_dir(filenames[0])
exp_labels=get_exp_labels(filenames)
data = pd.DataFrame()
data_flat = pd.DataFrame()
data_flat_3sd_cleaned=pd.DataFrame()
data_flat_mad_cleaned=pd.DataFrame()
var_data_flat = pd.DataFrame()
var_data_flat_3sd_cleaned = pd.DataFrame()
var_data_flat_mad_cleaned = pd.DataFrame()
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
    # do flattening within given non-overlapping flat_windows (in seconds)
    y_fit, y_flat = flattening(df[y_name], fps, chunk_l=flat_window)
    x_fit, x_flat = flattening(df[x_name], fps, chunk_l=flat_window)
    # plot_flattening(df['t'], df[x_name], x_flat, x_fit, ch_l=flat_window, title='x signal flattening')
    # plot_flattening(df['t'], df[y_name], y_flat, y_fit, ch_l=flat_window, title='y signal flattening')
    # now remove spikes of the timeseries data: those that are outside of 3SD from mean
    x_flat_3sd_cleaned, _ = rem_spikes_one_mean(x_flat)
    y_flat_3sd_cleaned, _ = rem_spikes_one_mean(y_flat)
    x_flat_mad_cleaned, _ = rem_spikes_one_median(x_flat)
    y_flat_mad_cleaned, _ = rem_spikes_one_median(y_flat)
    # y_norm_fit, y_norm_flat = flattening(df[y_norm_name], fps, chunk_l=flat_window)
    # x_norm_fit, x_norm_flat = flattening(df[x_norm_name], fps, chunk_l=flat_window)
    var_x = var_calc(data=x_flat, fps=fps, chunk_l=var_window)
    var_y = var_calc(data=y_flat, fps=fps, chunk_l=var_window)
    # also calculate variance for the cleaned data (without spikes)
    var_x_3sd_cleaned = var_calc(data=x_flat_3sd_cleaned, fps=fps, chunk_l=var_window)
    var_y_3sd_cleaned = var_calc(data=y_flat_3sd_cleaned, fps=fps, chunk_l=var_window)
    var_x_mad_cleaned = var_calc(data=x_flat_mad_cleaned, fps=fps, chunk_l=var_window)
    var_y_mad_cleaned = var_calc(data=y_flat_mad_cleaned, fps=fps, chunk_l=var_window)
    print(var_x.mean())
    print(var_y.mean())
    # var_x_norm = var_calc(data=x_norm_flat, fps=fps, chunk_l=var_window)
    # var_y_norm = var_calc(data=y_norm_flat, fps=fps, chunk_l=var_window)
    data_flat = pd.concat([data_flat,
                           pd.DataFrame({'Iz': y_flat,
                                          'Il': x_flat,
                                          'experiment': np.full(len(y_flat), exp_labels[i])})],
                          axis=0)
    data_flat_3sd_cleaned = pd.concat([data_flat_3sd_cleaned,
                           pd.DataFrame({'Iz': y_flat_3sd_cleaned,
                                          'Il': x_flat_3sd_cleaned,
                                          'experiment': np.full(len(y_flat_3sd_cleaned), exp_labels[i])})],
                          axis=0)
    data_flat_mad_cleaned = pd.concat([data_flat_mad_cleaned,
                           pd.DataFrame({'Iz': y_flat_mad_cleaned,
                                          'Il': x_flat_mad_cleaned,
                                          'experiment': np.full(len(y_flat_mad_cleaned), exp_labels[i])})],
                          axis=0)
    var_data_flat = pd.concat([var_data_flat,
                               pd.DataFrame({'var_x': var_x,
                                             'var_y': var_y,
                                             'experiment': np.full(len(var_y), exp_labels[i])})]).reset_index()
    var_data_flat_3sd_cleaned = pd.concat([var_data_flat_3sd_cleaned,
                               pd.DataFrame({'var_x': var_x_3sd_cleaned,
                                             'var_y': var_y_3sd_cleaned,
                                             'experiment': np.full(len(var_y_3sd_cleaned), exp_labels[i])})]).reset_index()
    var_data_flat_mad_cleaned = pd.concat([var_data_flat_mad_cleaned,
                               pd.DataFrame({'var_x': var_x_mad_cleaned,
                                             'var_y': var_y_mad_cleaned,
                                             'experiment': np.full(len(var_y_mad_cleaned), exp_labels[i])})]).reset_index()
    i+=1
data_flat =pd.DataFrame({'t': np.arange(0, len(data_flat)/fps, 1/fps),
                         'Iz': data_flat['Iz'],
                         'Il': data_flat['Il'],
                         'experiment': data_flat['experiment']})
data_flat_3sd_cleaned =pd.DataFrame({'t': np.arange(0, len(data_flat_3sd_cleaned)/fps, 1/fps),
                         'Iz': data_flat_3sd_cleaned['Iz'],
                         'Il': data_flat_3sd_cleaned['Il'],
                         'experiment': data_flat_3sd_cleaned['experiment']})
data_flat_mad_cleaned =pd.DataFrame({'t': np.arange(0, len(data_flat_mad_cleaned)/fps, 1/fps),
                         'Iz': data_flat_mad_cleaned['Iz'],
                         'Il': data_flat_mad_cleaned['Il'],
                         'experiment': data_flat_mad_cleaned['experiment']})
var_data_flat =pd.DataFrame({'t': np.arange(0, len(var_data_flat)*var_window, var_window),
                         'var_x': var_data_flat['var_x'],
                         'var_y': var_data_flat['var_y'],
                         'var_sum': var_data_flat['var_x']+var_data_flat['var_y'],
                         'experiment': var_data_flat['experiment']}).reset_index()
var_data_flat_3sd_cleaned =pd.DataFrame({'t': np.arange(0, len(var_data_flat_3sd_cleaned)*var_window, var_window),
                         'var_x': var_data_flat_3sd_cleaned['var_x'],
                         'var_y': var_data_flat_3sd_cleaned['var_y'],
                         'var_sum': var_data_flat_3sd_cleaned['var_x']+var_data_flat_3sd_cleaned['var_y'],
                         'experiment': var_data_flat_3sd_cleaned['experiment']}).reset_index()
var_data_flat_mad_cleaned =pd.DataFrame({'t': np.arange(0, len(var_data_flat_mad_cleaned)*var_window, var_window),
                         'var_x': var_data_flat_mad_cleaned['var_x'],
                         'var_y': var_data_flat_mad_cleaned['var_y'],
                         'var_sum': var_data_flat_mad_cleaned['var_x']+var_data_flat_3sd_cleaned['var_y'],
                         'experiment': var_data_flat_mad_cleaned['experiment']}).reset_index()


plot_multiple_sigs('t', 'Iz', data=data_flat, title='flattened y deflection')
plot_multiple_sigs('t', 'Il', data=data_flat, title='flattened x deflection')
plot_multiple_sigs('t', 'Iz', data=data_flat_3sd_cleaned, title='flattened y deflection 3SD cleaned' )
plot_multiple_sigs('t', 'Il', data=data_flat_3sd_cleaned, title='flattened x deflection 3SD cleaned')
plot_multiple_sigs('t', 'Iz', data=data_flat_mad_cleaned, title='flattened y deflection mad cleaned' )
plot_multiple_sigs('t', 'Il', data=data_flat_mad_cleaned, title='flattened x deflection mad cleaned')
plot_multiple_sigs('t', 'var_x', data=var_data_flat, title='variance of flattened x deflection')
plot_multiple_sigs('t', 'var_y', data=var_data_flat, title='variance of flattened y deflection')
plot_multiple_sigs('t', 'var_sum', data=var_data_flat, title='sum variance of flattened signals')
plot_multiple_sigs('t', 'var_x', data=var_data_flat_3sd_cleaned, title='variance of x flattened 3SD cleaned')
plot_multiple_sigs('t', 'var_y', data=var_data_flat_3sd_cleaned, title='variance of y flattened 3SD cleaned')
plot_multiple_sigs('t', 'var_sum', data=var_data_flat_3sd_cleaned, title='sum variance of signals flattened 3SD cleaned')
plot_multiple_sigs('t', 'var_x', data=var_data_flat_mad_cleaned, title='variance of x flattened mad cleaned')
plot_multiple_sigs('t', 'var_y', data=var_data_flat_mad_cleaned, title='variance of y flattened mad cleaned')
plot_multiple_sigs('t', 'var_sum', data=var_data_flat_mad_cleaned, title='sum variance of signals flattened mad cleaned')
var_data_norm = var_data_flat.copy()
var_x_mean = var_data_norm[var_data_norm['experiment']==var_data_norm['experiment'].iloc[0]]['var_x'].mean()
var_data_norm['var_x']=var_data_norm['var_x']/var_x_mean
var_y_mean = var_data_norm[var_data_norm['experiment']==var_data_norm['experiment'].iloc[0]]['var_y'].mean()
var_data_norm['var_y']=var_data_norm['var_y']/var_y_mean
# plot_var_bars(data=var_data_norm, var_label='var_x')
# plot_var_bars(data=var_data_norm, var_label='var_y')
plot_var_bars_mean(data=var_data_flat, var_label='var_x', title='barplot variance x flattened')
plot_var_bars_mean(data=var_data_flat, var_label='var_y', title='barplot variance y flattened')
plot_var_bars_mean(data=var_data_flat, var_label='var_sum', title='barplot sum of variance of flattened signals')
plot_var_bars_med(data=var_data_flat, var_label='var_x', exp_labels=exp_labels, title='barplot median variance x flattened')
plot_var_bars_med(data=var_data_flat, var_label='var_y', exp_labels=exp_labels, title='barplot median variance y flattened')
plot_var_bars_med(data=var_data_flat, var_label='var_sum', exp_labels=exp_labels, title='barplot median sum of variance of flattened signals')
plot_var_boxplot(data=var_data_flat, var_label='var_x', title='variance x flattened (boxplot)')
plot_var_boxplot(data=var_data_flat, var_label='var_y', title='variance y flattened (boxplot)')
plot_var_boxplot(data=var_data_flat, var_label='var_sum', title='sum of variance of flattened signals (boxplot)')

plot_var_bars_mean(data=var_data_flat_3sd_cleaned, var_label='var_x', title='barplot variance x, flattened, 3SD cleaned ')
plot_var_bars_mean(data=var_data_flat_3sd_cleaned, var_label='var_y', title='barplot variance y flattened, 3SD cleaned ')
plot_var_bars_mean(data=var_data_flat_3sd_cleaned, var_label='var_sum', title='barplot sum of variance of flattened, 3SD cleaned signals')
plot_var_bars_med(data=var_data_flat_3sd_cleaned, var_label='var_x', exp_labels=exp_labels, title='barplot median variance x, flattened, 3SD cleaned ')
plot_var_bars_med(data=var_data_flat_3sd_cleaned, var_label='var_y', exp_labels=exp_labels, title='barplot median variance y flattened, 3SD cleaned ')
plot_var_bars_med(data=var_data_flat_3sd_cleaned, var_label='var_sum', exp_labels=exp_labels, title='barplot median sum of variance of flattened, 3SD cleaned signals')
plot_var_boxplot(data=var_data_flat_3sd_cleaned, var_label='var_x', title='boxplot variance x flattened, 3SD cleaned (boxplot)')
plot_var_boxplot(data=var_data_flat_3sd_cleaned, var_label='var_y', title='boxplot variance y flattened, 3SD cleaned (boxplot)')
plot_var_boxplot(data=var_data_flat_3sd_cleaned, var_label='var_sum', title='boxplot sum of variance of flattened, 3SD cleaned signals (boxplot)')

plot_var_bars_med(data=var_data_flat_mad_cleaned, var_label='var_x', exp_labels=exp_labels, title='barplot median variance x, flattened, mad cleaned ')
plot_var_bars_med(data=var_data_flat_mad_cleaned, var_label='var_y', exp_labels=exp_labels, title='barplot median variance y flattened, mad cleaned ')
plot_var_bars_med(data=var_data_flat_mad_cleaned, var_label='var_sum', exp_labels=exp_labels, title='barplot median sum of variance of flattened, mad cleaned signals')
plot_var_boxplot(data=var_data_flat_mad_cleaned, var_label='var_x', title='boxplot variance x flattened, mad cleaned')
plot_var_boxplot(data=var_data_flat_mad_cleaned, var_label='var_y', title='boxplot variance y flattened, mad cleaned')
plot_var_boxplot(data=var_data_flat_mad_cleaned, var_label='var_sum', title='boxplot sum of variance of flattened, mad cleaned signals')




    
    