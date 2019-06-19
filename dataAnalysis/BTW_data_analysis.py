"""
Programmer: Sijun Shen
Date: June 19, 2019
The dimension_transformer class is used to convert the phone axes to vehicle axes
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyquaternion import Quaternion
import math

class dimension_transformer:
    def __init__(self, motion_file, location_file):
        self.motion_data = pd.read_csv(motion_file)
        self.motion_data['time'] = self.motion_data.ZTIMESTAMPEST.str[28:30]
        self.motion_data['time'] =  pd.to_numeric(self.motion_data['time'], errors='coerce')
        self.motion_data = self.motion_data[self.motion_data['time'] < 33]
        self.location_data = pd.read_csv(location_file)

    def phone_frame_to_reference_frame_helper(self, w, x, y, z, tran_x, tran_y, tran_z):
        my_quaternion = Quaternion(q1=w, q2=x, q3=y, q4=z)
        v = np.array([tran_x, tran_y, tran_z])
        v_prime = my_quaternion.rotate(v)
        return v_prime

    # function to convert angle to radian
    def convert_to_radian(self, angle):
        return angle * 3.1415926 / 180

    # rotate 2D points around the origin clockwise
    def rotate_point(self, x_origin, y_origin, angle):
        radian = self.convert_to_radian(angle)
        x_new = x_origin * math.cos(radian) - y_origin * math.sin(radian)
        y_new = x_origin * math.sin(radian) + y_origin * math.cos(radian)
        return x_new, y_new

    # convert phone frame to reference frame
    def phone_frame_to_reference_frame(self):
        outcome_temp = self.motion_data.apply(lambda x: self.phone_frame_to_reference_frame_helper(x.ZQUATW, x.ZQUATX, x.ZQUATY, x.ZQUATZ, x.ZGRAVITYX, x.ZGRAVITYY, x.ZGRAVITYZ), axis = 1)
        outcome_m = np.matrix(outcome_temp.tolist())
        df = pd.DataFrame(data=outcome_m)
        df.columns = ['gravity_original_x', "gravity_original_y", "gravity_original_z"]
        self.motion_data = pd.concat([self.motion_data, df], axis=1)

        outcome_temp = self.motion_data.apply(lambda x: self.phone_frame_to_reference_frame_helper(x.ZQUATW, x.ZQUATX, x.ZQUATY, x.ZQUATZ, x.ZUSERACCELX, x.ZUSERACCELY, x.ZUSERACCELZ), axis = 1)
        outcome_m = np.matrix(outcome_temp.tolist())
        df = pd.DataFrame(data=outcome_m)
        df.columns = ['user_original_x', "user_original_y", "user_original_z"]
        self.motion_data = pd.concat([self.motion_data, df], axis = 1)

    # convert reference frame to vehicle frame
    def reference_frame_to_vehicle_frame(self):
        self.location_data = self.location_data[self.location_data['ZCOURSE'] > 0]
        interpolated_zcourse = np.interp(self.motion_data.ZTIMESTAMPEPOCH, self.location_data.ZTIMESTAMPEPOCH, self.location_data.ZCOURSE)
        self.motion_data['interpolated_zcourse'] = interpolated_zcourse

        outcome_temp = self.motion_data.apply(lambda x: self.rotate_point(x.user_original_x, x.user_original_y, x.interpolated_zcourse), axis = 1)
        outcome_m = np.matrix(outcome_temp.tolist())
        df = pd.DataFrame(data=outcome_m)
        df.columns = ['vehicle_trajectory', 'vehicle_orthogonal']
        self.motion_data = pd.concat([self.motion_data, df], axis = 1)

    # function used to plot figures
    def plot_figure(self, col_name, title):
        duplicated = self.motion_data.duplicated(subset = 'time',keep='first')
        x_axis_ix = self.motion_data.index[duplicated != True]
        x_axis_labels = self.motion_data.time[duplicated != True]

        periods = [[11, 13], [14, 17], [19, 23], [24, 27]]
        locs = ["N", "E", "W", "S"]
        periods_label = []

        for period in periods:
            periods_label.append([x_axis_ix[x_axis_labels == period[0]].values, x_axis_ix[x_axis_labels == period[1]].values])

        plt.figure()
        plt.plot(self.motion_data.index, self.motion_data[col_name], '-o', linewidth=0.5, markersize=0.5)
        plt.xticks(x_axis_ix, x_axis_labels, fontsize=8)
        plt.ylabel(title, size=10)
        plt.ylim([-1, 1])
        plt.xlabel("Time", size=10)
        plt.title(title, size=14)
        plt.axhline(y=0.3, color='r', linestyle='--')
        plt.axhline(y=-0.3, color='r', linestyle='--')

        for i in range(len(periods_label)):
            plt.axvline(x=periods_label[i][0], color = "y", linestyle="-.")
            plt.axvline(x=periods_label[i][1], color="y", linestyle="-")
            x_loc = periods_label[i][0] / 2.0 + periods_label[i][1] / 2.0 - 20
            plt.annotate(locs[i], (x_loc[0], 0.8))

        plt.savefig(title)
        plt.show()

if __name__ == "__main__":
    transformer = dimension_transformer("ZDEVICEMOTION_SUNDAY06162019.csv", "ZLOCATIONS_SUNDAY06162019.csv")
    transformer.phone_frame_to_reference_frame()
    transformer.reference_frame_to_vehicle_frame()

    transformer.plot_figure("ZGRAVITYX", title="Gravity on Phone Frame x")
    transformer.plot_figure("ZGRAVITYY", title="Gravity on Phone Frame y")
    transformer.plot_figure("ZGRAVITYZ", title="Gravity on Phone Frame z")

    transformer.plot_figure("gravity_original_x", title="Gravity on Reference Frame x (north)")
    transformer.plot_figure("gravity_original_y", title="Gravity on Reference Frame y (west)")
    transformer.plot_figure("gravity_original_z", title="Gravity on Reference Frame z (vertical)")

    transformer.plot_figure("user_original_x", title="User acceleration on Reference Frame x (north)")
    transformer.plot_figure("user_original_y", title="User acceleration on Reference Frame y (west)")
    transformer.plot_figure("user_original_z", title="User acceleration on Reference Frame z (vertical)")

    transformer.plot_figure("vehicle_trajectory", title="Acceleration on vehicle trajectory")
    transformer.plot_figure("vehicle_orthogonal", title="Acceleration orthogonal to vehicle trajectory ")







