"""
Programmer: Sijun Shen
Date: Aug 7th, 2019
The PCA_Transformer class is used to convert the phone axes to vehicle axes using
principal component analysis algorithm
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyquaternion import Quaternion
import math

# class PCA_Trasnformer():
#     def __init__(self, time_window=10):
#         self.transformer =
