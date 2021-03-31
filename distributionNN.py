import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import sys

#hist1=[14212197,1769672,929304,544091,355937,264090,208551,210142,176133,234959]
hist1=[       0,        0,        0,        0,        0,        0,        0,        0,
              0, 11050114,  2809139,  1261838,   797054,   575444,   429672,   320840,
              247362,   198500,   166455,   143776,   126041,   112363,    96280,   135002,
              81964,    98320,   106202,   132027,    16683]
#hist2=[55774,48301,46551,42804,40074,42290,47088,76380,109776,1784001]
hist2=[      0,       0,       0,       0,       0,       0,       0,
             0,       0,   20133,   29602,   26469,   25403,   25087,
         24685,   23448,   22102,   21280,   21582,   22492,   23544,
         25144,   26581,   50852,   39984,   64053,  114364,  542972,
       1143262]

print("0 evt:",sum(hist1))
print("1 evt:",sum(hist2))
bin_edges=[-0.5       , -0.44655172, -0.39310345, -0.33965517, -0.2862069 ,
       -0.23275862, -0.17931034, -0.12586207, -0.07241379, -0.01896552,
        0.03448276,  0.08793103,  0.14137931,  0.19482759,  0.24827586,
        0.30172414,  0.35517241,  0.40862069,  0.46206897,  0.51551724,
        0.56896552,  0.62241379,  0.67586207,  0.72931034,  0.78275862,
        0.8362069 ,  0.88965517,  0.94310345,  0.99655172,  1.05      ]
#[0. , 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 ]
bin_center=[]
for a in range(len(bin_edges)-1):
    bin_center.append((bin_edges[a]+bin_edges[a+1])/2)
ax7 = plt.axes()
ax7.plot(bin_center,hist1,label='notnumuCC')
ax7.plot(bin_center,hist2,label='(a)numuCC')
ax7.plot(bin_center,hist1,'ok')
ax7.plot(bin_center,hist2,'ok')
ax7.legend()
ax7.set_xlabel('NN output')
ax7.set_yscale('log')
plt.show()
