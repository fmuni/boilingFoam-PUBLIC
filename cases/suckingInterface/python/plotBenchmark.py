# Plot the content of a csv file; file based on H. Scheufler's library https://github.com/DLR-RY/TwoPhaseFlow
import scipy
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
import csv
import seaborn as sns
import pandas as pd

rcParams.update({'figure.autolayout': True})
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True

def plotdiameter():
    # Initialize empty arrays

    t1 = []
    x1 = []
    t2 = []
    x2 = []



    with open('python/data.dat','r') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        next(data)
        for row in data:
            t1.append(float(row[0]))
            x1.append(float(row[2]))

    with open('bubbleInfo.csv','r') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        for row in data:
            t2.append(float(row[0]))
            x2.append((float(row[1])-5e-12)/1e-8)

    # Plot
    plt.figure(figsize=(250 /25.4, 200 / 25.4))
    plt.plot(t1,x1, 'b--', label='Analytical solution',linewidth=3)
    plt.plot(t2,x2, 'r-', label='Simulation',linewidth=3)
    plt.xlabel('Time [s]')
    plt.ylabel('Position [m]')
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #plt.title('Interesting Graph\nCheck it out')
    plt.legend(frameon=True,loc='upper left')
    #plt.show()
    plt.draw()
    # Save in pdf format
    plt.savefig('suckingInterface.pdf')


plotdiameter()
plt.show()
