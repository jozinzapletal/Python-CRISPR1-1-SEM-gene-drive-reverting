import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler


def generate_plots():

    generations = 60

    data1 = pd.read_excel('CRISPR1_1 Results.xlsx')
    male_gammaList = data1['male_gamma'].unique()
    male_alphaList = data1['male_alpha'].unique()
    female_gammaList = data1['female_gamma'].unique()
    female_alphaList = data1['female_alpha'].unique()
    xPlot = data1['time'].unique()

    ft = 6
    plot_count = 1

    alleleList = ['w', 'v', 'u', 'r', 'g', 's']

    for a1 in male_alphaList:
        for a2 in female_alphaList:
            for g1 in male_gammaList:
                for g2 in female_gammaList:

                    fig1 = plt.figure(num=None, figsize=(3, 2), dpi=300)

                    for b in alleleList:

                        yPlot = data1[(data1['male_gamma'] == g1) & (data1['male_alpha'] == a1) &
                                      (data1['female_gamma'] == g2) & (data1['female_alpha'] == a2)][b]

                        colors = ['steelblue', 'deepskyblue', 'darkviolet', 'orchid', 'firebrick', 'salmon']
                        plt.rc('axes', prop_cycle=(cycler('color', colors) + cycler('linestyle', ['-', '-', '-', '-', '-', '-'])))

                        plt.ylim(0, 1)
                        plt.xlim(0, xPlot.max())
                        plt.ylabel('Allele Frequency', fontsize=ft)
                        plt.xlabel('Time (generations)', fontsize=ft)

                        ax = plt.subplot(111)
                        ax.spines['right'].set_visible(False)
                        ax.spines['top'].set_visible(False)
                        ax.yaxis.set_ticks_position('left')
                        ax.xaxis.set_ticks_position('bottom')
                        ax.set_xticks([0, 10, 20, 30, 40, 50, 60])
                        ax.set_xticklabels(['0', '10', '20', '30', '40', '50', '60'])

                        plt.tick_params(axis='both', which='major', labelsize=ft)
                        ax.tick_params(direction="in")

                        plt.plot(xPlot, yPlot.T, label=b)

                        plt.legend(fontsize=ft)
                        plt.legend(loc='right', prop={'size': 6}, frameon=False, bbox_to_anchor=(1.05, 0.8))

                        plt.text(generations / 4, 1.05, 'male (α = ' + str(a1) + ', γ = ' + str(g1) + ')\n' +
                                                        'female (α = ' + str(a2) + ', γ = ' + str(g2) + ')',
                                 color='black', fontsize=ft)

                        plt.subplots_adjust(left=0.14, bottom=0.16, right=0.90, top=0.80, wspace=0, hspace=0)

                    plt.savefig(str(plot_count) + '.png')
                    plt.show()
                    plot_count += 1
