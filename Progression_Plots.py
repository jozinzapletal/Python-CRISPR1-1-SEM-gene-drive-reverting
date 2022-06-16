import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler

generations = 60;
developmentTime = 12;

data1 = pd.read_excel('CRISPR1_1 Results.xlsx')

gammaList = data1['gamma'].unique()
alphaList = data1['alpha'].unique()
xPlot = data1['time'].unique()

ft = 6

for g in gammaList:

    fig1 = plt.figure(num=None, figsize=(2.25, 1.83), dpi=300)

    for a in alphaList:
        ax = plt.subplot(111)

        yPlot = data1[(data1['gamma'] == g) & (data1['alpha'] == a)]
        yPlot = yPlot.drop(['gamma', 'alpha', 'total_population'], axis=1)

        colors = ['firebrick', 'red', 'coral', 'lightcoral', 'darkviolet', 'orchid']
        plt.rc('axes', prop_cycle=(cycler('color', colors) + cycler('linestyle', ['-', '-', '-', '-', '-', '-'])))

        plt.plot(xPlot, yPlot, label=u'α = ' + str(a))
        plt.ylim(0, 1)
        plt.xlim(0, xPlot.max())
        plt.ylabel('Proportion of Wildtype', fontsize=ft)
        plt.xlabel('Time (generations)', fontsize=ft)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xticks([0, 120, 240, 360, 480, 600, 720])
        ax.set_xticklabels(['0', '10', '20', '30', '40', '50', '60'])

        plt.tick_params(axis='both', which='major', labelsize=ft)
        ax.tick_params(direction="in")

        plt.legend(fontsize=ft)
        plt.legend(loc='right', prop={'size': 6}, frameon=False, bbox_to_anchor=(1.4, 0.7))

        plt.text(generations / 2.5, 1.05, 'γ = ' + str(g), color='black', fontsize=ft)

        plt.show()