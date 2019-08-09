"""
Programmer: Sijun Shen
Date: Aug 8th, 2019
The Comparator is used to compare the performance of Dimension_Transformer and PCA_Transformer
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Dimension_Transformer as DT
from sklearn.decomposition import PCA
import os

class Comparator:
    def __init__(self, dataSet='motion_data.csv'):
        self.df = pd.read_csv(dataSet)
        self.df = self.df[self.df['time'].isin([11, 12, 14, 15, 16, 19, 20, 21, 22, 24, 25, 26])]

    def plot_figure(self, column_name):
        duplicated = self.df.duplicated(subset = 'time',keep='first')
        x_axis_ix = self.df.index[duplicated != True]
        x_axis_labels = self.df.time[duplicated != True]

        corr1 = "Corr with vehicle_trajectory: " + str(
            round(self.df[column_name].corr(self.df['vehicle_trajectory']), 2))
        corr2 = "Corr with vehicle_orthogonal: " + str(
            round(self.df[column_name].corr(self.df['vehicle_orthogonal']), 2))

        plt.figure()
        plt.plot(self.df.index, self.df[column_name], 'y-', label=column_name)
        plt.plot(self.df.index, self.df['vehicle_trajectory'], 'r--', label='vehicle trajectory')
        plt.plot(self.df.index, self.df['vehicle_orthogonal'], 'g-.', label='vehicle orthogonal')
        plt.legend(loc='lower right')
        plt.ylabel("G-force", size=10)
        plt.ylim([-1, 1])
        plt.xticks(x_axis_ix, x_axis_labels, fontsize=8)
        plt.xlabel("time")
        plt.axhline(y=0.3, color='black', linestyle='--')
        plt.axhline(y=-0.3, color = 'black', linestyle='--')
        plt.annotate('N', (x_axis_ix[x_axis_labels == 12][0], -0.4))
        plt.annotate('E', (x_axis_ix[x_axis_labels == 15][0], -0.4))
        plt.annotate('W', (x_axis_ix[x_axis_labels == 21][0], -0.4))
        plt.annotate('S', (x_axis_ix[x_axis_labels == 25][0], -0.4))
        plt.annotate(corr1, (x_axis_ix[x_axis_labels == 12][0], -.6))
        plt.annotate(corr2, (x_axis_ix[x_axis_labels == 12][0], -.7))
        plt.title(column_name, size=14)
        plt.savefig(os.path.join('plots', column_name))
        plt.show()
        print("Corr with vehicle_trajectory: " + str(round(self.df[column_name].corr(self.df['vehicle_trajectory']), 2)))

if __name__ == "__main__":
    comparator = Comparator("motion_data.csv")
    comparator.plot_figure("PCA_x")
    comparator.plot_figure("PCA_y")

    comparator = Comparator("motion_data_PLS.csv")
    comparator.plot_figure("PLS_x")
    comparator.plot_figure("PLS_y")






