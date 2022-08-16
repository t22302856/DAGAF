#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data-adaptive Gaussian average filtering on signal s.
Designed by Yue-Der Lin (ydlin@fcu.edu.tw) and Kai-Chun Liu (t22302856@), Taiwan.

Reference: https://doi.org/10.1016/J.BSPC.2021.103104
     
--------Simple Example----------------     
"""


import numpy as np
import matplotlib.pyplot as plt
from dagaf import dagaf
s = np.random.random(1000)

# Plot the data.
plt.figure()
plt.plot(s)
plt.show()
imf, res = dagaf(s, 4, 1.6, 'd', 'd', 20, 0.001, True)
