# the goal: look at the ics for the tip, stalk, tip/stalk
# take the percentage of "on" for the stability in the one point cycles
# this will tell us how much of a contribuition does that node play in
# the phenotype adoption

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# ics = [notch, delta, jagged, icd, vegfr, vegf, d_ext, j_ext, n_ext, v_ext] -- bar code definition

tips = pd.read_csv('Tip Intitial Conditions Synchronous Bin', usecols = [1])
stalks = pd.read_csv('Stalk Intitial Conditions Synchronous Bin', usecols = [1])
hybrids = pd.read_csv('Hybrid Intitial Conditions Synchronous Bin', usecols = [1])
tipstalks = pd.read_csv('Tip_Stalk Intitial Conditions Synchronous Bin', usecols = [1])
print(tips)
test_tip = tips.get_values()

tip = tips.to_numpy()
stalk = stalks.to_numpy()
hybrid = hybrids.to_numpy()
tipstalk = tipstalks.to_numpy()
print()
print(test_tip)
print()
print(tip)
print(tip.shape)
d = np.zeros(len(tip))
n = np.zeros(len(tip))
j, icd, vr, v, dext, next, jext, vext = d[:], d[:], d[:], d[:], d[:], d[:], d[:], d[:]
for i in range(len(tip)):
    place_hold = tip[i]
    a = np.fromstring(place_hold[:],dtype=np.int,sep=' ')
    print(a)
    print(place_hold)
    row = np.array(place_hold[0])
    print("np array row:",row)
    print(row)
    print(row)
    n[i], d[i], j[i], icd[i], vr[i], v[i], dext[i], jext[i], next[i], vext[i] = row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9]
print(n)
print()
print(d)
