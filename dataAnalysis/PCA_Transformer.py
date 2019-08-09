"""
Programmer: Sijun Shen
Date: Aug 7th, 2019
The PCA_Transformer class is used to convert the phone axes to vehicle axes using
principal component analysis algorithm
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Dimension_Transformer as DT
from sklearn.decomposition import PCA
import os

class PCA_Transformer():
    def __init__(self):
        self.transformer = DT.Dimension_Transformer("ZDEVICEMOTION_SUNDAY06162019.csv", "ZLOCATIONS_SUNDAY06162019.csv")
        self.transformer.phone_frame_to_reference_frame()
        self.transformer.reference_frame_to_vehicle_frame()
        self.transformer.motion_data['PCA_x'] = 0
        self.transformer.motion_data['PCA_y'] = 0
        self.df = self.transformer.motion_data # to save the space of code

    def transform(self, time_window = 20):
        cur_index = 0
        next_index = cur_index + time_window
        features = ['user_original_x', 'user_original_y']
        while next_index <= len(self.df.index):
            data_temp = self.df.iloc[cur_index : next_index, :]
            res = self.transform_helper(data_temp.loc[:, features])
            self.df.loc[cur_index : (next_index - 1), 'PCA_x'] = res['PC_1']
            self.df.loc[cur_index : (next_index - 1), 'PCA_y'] = res['PC_2']
            cur_index = next_index
            next_index = cur_index + time_window

    def transform_helper(self, data):
        data_temp = data - data.mean(axis=0) # make the mean of columns equal to zero
        # data_temp = data_temp.to_numpy()
        pca = PCA(n_components=2)
        pca.fit(data_temp)
        # access and vectors
        vectors = pca.components_.T
        res = data.dot(vectors)
        res.columns = ['PC_1', 'PC_2']
        return res

    def plot_figure(self, col_name, title):
        duplicated = self.df.duplicated(subset='time', keep='first')
        x_axis_ix = self.df.index[duplicated != True]
        x_axis_labels = self.df.time[duplicated != True]

        periods = [[11, 13], [14, 17], [19, 23], [24, 27]]
        locs = ["N", "E", "W", "S"]
        periods_label = []

        for period in periods:
            periods_label.append(
                [x_axis_ix[x_axis_labels == period[0]].values, x_axis_ix[x_axis_labels == period[1]].values])

        plt.figure()
        plt.plot(self.df.index, self.df[col_name], '-o', linewidth=0.5, markersize=0.5)
        plt.xticks(x_axis_ix, x_axis_labels, fontsize=8)
        plt.ylabel(title, size=10)
        plt.ylim([-1, 1])
        plt.xlabel("Time", size=10)
        plt.title(title, size=14)
        plt.axhline(y=0.3, color='r', linestyle='--')
        plt.axhline(y=-0.3, color='r', linestyle='--')

        for i in range(len(periods_label)):
            plt.axvline(x=periods_label[i][0], color = "y", linestyle="-")
            plt.axvline(x=periods_label[i][1], color="y", linestyle="-")
            x_loc = periods_label[i][0] / 2.0 + periods_label[i][1] / 2.0 - 20
            plt.annotate(locs[i], (x_loc[0], 0.8))

        plt.savefig(os.path.join('plots', title))
        plt.show()

if __name__ == "__main__":
    pca_transformer = PCA_Transformer()
    pca_transformer.transform()
    pca_transformer.plot_figure('PCA_x', "PCA on X")
    pca_transformer.plot_figure('PCA_y', "PCA on Y")
    pca_transformer.df.to_csv(r'C:\Users\Owner\Dropbox\Sijun_Research\BTW analysis\dataAnalysis\motion_data.csv', index = None, header=True)








