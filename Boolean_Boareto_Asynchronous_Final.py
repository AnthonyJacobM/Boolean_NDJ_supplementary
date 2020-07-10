# testing new rules for boolean boareto
import boolean2
from boolean2 import Model, util
import matplotlib; matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
import random
from pylab import *

# 2: n *= i or not (n and dext) or not (n and jext) or not (n and j) or not (n and d)
# 2: n *= i or not (n and (d or j))


proteins = ['n', 'd', 'j', 'i', 'vr', 'v']
n_and_d = []
n_not_d = []
d_not_n = []
not_d_not_n = []
tip_bin = []
stalk_bin = []
hybrid_bin = []
protein_bin = []

n_and_d_j = []
n_not_d_j = []
d_not_n_j = []
not_d_not_n_j = []
tip_j_bin = []
stalk_j_bin = []
hybrid_j_bin = []
other_bin = []
other_j_bin = []
protein_j_bin = []


def ics(proteins, ko, oe):
    # input proteins and k/o and o/e arrays
    ics = []
    off = []
    on = []
    for i in range(len(ko)):
        off.append(ko[i])
    for j in range(len(oe)):
        on.append(oe[j])
    text_new = boolean2.modify_states(text=text, turnoff=off, turnon=on)
    return text_new


def inf_norm(x, y):
    assert len(x) == len(y)
    n = len(x)
    max = 0
    z = np.subtract(x, y)
    for i in range(n):
        if np.abs(z[i]) > max:
            max = np.abs(z[i])
    return max


def steady_state(x, y, error):
    distance = inf_norm(x, y)
    if distance <= error:
        return True
    elif distance > error:
        print("Unsteady State")
        return False


def asynch_threshold(protein, threshold=0.30):
    """ asynch determine if on or off """
    if protein >= threshold:
        return 1
    else:
        return 0


def asynch_threshold_list(proteins, threshold=0.25):
    bin = []
    if len(proteins) > 1:
        for i in range(len(proteins)):
            if proteins[i] >= threshold:
                bin.append(1)
            else:
                bin.append(0)
    return bin


def tip_f(icd, jagged, delta, vegfr, notch):
    if delta == 1 and vegfr == 1 and jagged == 0 and icd == 0 and notch == 0:
        return 1
    else:
        return 0


def stalk_f(icd, jagged, delta, vegfr, notch):
    if icd == 1 and jagged == 1 and delta == 0 and vegfr == 0 and notch == 1:
        return 1
    else:
        return 0


def hybrid_f(delta, vegfr, icd, jagged, notch):
    if delta == 1 and vegfr == 1 and icd == 1 and jagged == 1 and notch == 1:
        return 1
    else:
        return 0


def other(tip, hybrid, stalk):
    if tip == 0 and hybrid == 0 and stalk == 0:
        return 1
    else:
        return 0


coll = util.Collector()
flag = True
N = 10
sim = int(30)
t0 = 0
choices_list = ['cis_delta', 'cis_jagged', 'trans_delta', 'trans_jagged']
# while flag == True:

d_ss = []
j_ss = []
i_ss = []
n_ss = []
vr_ss = []
v_ss = []

d_ss_j = []
j_ss_j = []
i_ss_j = []
n_ss_j = []
vr_ss_j = []
v_ss_j = []

notch_ss = []
delta_ss = []
jagged_ss = []
icd_ss = []
vegfr_ss = []
vegf_ss = []

notch_ss2 = []
delta_ss2 = []
jagged_ss2 = []
icd_ss2 = []
vegfr_ss2 = []
vegf_ss2 = []

notch_ss3 = []
delta_ss3 = []
jagged_ss3 = []
icd_ss3 = []
vegfr_ss3 = []
vegf_ss3 = []

d_ics = []
j_ics = []
n_ics = []
i_ics = []
vr_ics = []
v_ics = []
dext_ics = []
jext_ics = []
next_ics = []
vext_ics = []

dext_ics2 = []
jext_ics2 = []
next_ics2 = []
vext_ics2 = []

dext_ics3 = []
jext_ics3 = []
next_ics3 = []
vext_ics3 = []

tip_ss = []
stalk_ss = []
tip_stalk_ss = []
hybrid_ss = []

tip_ss2 = []
stalk_ss2 = []
tip_stalk_ss2 = []
hybrid_ss2 = []

tip_ss3 = []
stalk_ss3 = []
tip_stalk_ss3 = []
hybrid_ss3 = []

n_bin = []
d_bin = []
j_bin = []
i_bin = []
vr_bin = []
v_bin = []

t_bin = []
s_bin = []
t_s_bin = []
h_bin = []

ics_n_j = []
ics_n = []
ics_d_j = []
ics_d = []
ics_j_j = []
ics_j = []
ics_i_j = []
ics_i = []
ics_vr = []
ics_vr_j = []
ics_v_j = []
ics_v = []
ics_vext_j = []
ics_vext = []
ics_next_j = []
ics_next = []
ics_dext_j = []
ics_dext = []
ics_jext_j = []
ics_jext = []

ics_bin = ['n', 'd', 'j', 'vr', 'v', 'dext', 'jext', 'next', 'vext']
R = 2 ** (len(ics_bin))
step = 0
ss_num = 0
ss_num2 = 0
ss_num3 = 0
unsteady_num = 0
choices = ['cis_jagged', 'cis_delta', 'trans_delta', 'trans_jagged']

# change the initial conditions
import pylab

total_tip = 0
total_stalk = 0
total_tip_stalk = 0
total_hybrid = 0
total_tip_j = 0
total_stalk_j = 0
total_tip_stalk_j = 0
total_hybrid_j = 0
total_tip_k = 0
total_stalk_k = 0
total_tip_stalk_k = 0
total_hybrid_k = 0

asynch_bin = np.zeros((20, 6))
jagged_ics = []

for notch in range(0, 2):
    for delta in range(0, 2):
        for jagged in range(2):
            for vegfr in range(0, 2):
                for vegf in range(2):
                    for d_ext in range(2):
                        for j_ext in range(2):
                            for n_ext in range(2):
                                for v_ext in range(2):
                                    for icd in range(0, 2):
                                        step += 1
                                        if jagged == 1:
                                            j0 = "True"
                                        else:
                                            j0 = 'False'
                                        if delta == 1:
                                            d0 = "True"
                                        else:
                                            d0 = "False"
                                        if notch == 1:
                                            n0 = "True"
                                        else:
                                            n0 = "False"
                                        if vegfr == 1:
                                            vr0 = 'True'
                                        else:
                                            vr0 = 'False'
                                        if vegf == 1:
                                            v0 = 'True'
                                        else:
                                            v0 = 'False'
                                        if icd == 1:
                                            i0 = 'True'
                                        else:
                                            i0 = 'False'
                                        if d_ext == 1:
                                            dext0 = 'True'
                                        else:
                                            dext0 = 'False'
                                        if j_ext == 1:
                                            jext0 = 'True'
                                        else:
                                            jext0 = 'False'
                                        if n_ext == 1:
                                            next0 = 'True'
                                        else:
                                            next0 = 'False'
                                        if v_ext == 1:
                                            vext0 = 'True'
                                        else:
                                            vext0 = 'False'

                                        # define the logic gate:
                                        text_ics = """
                                        n = {}
                                        d = {}
                                        j = {}
                                        i = {}
                                        vr = {}
                                        v = {}

                                        tip = vr and d and not j and not i and not n
                                        stalk = i and j and n and not d and not vr
                                        tip_stalk = vr and d and j and n and i
                                        hybrid = not (tip or stalk or tip_stalk)

                                        jagged_oe = False
                                        jagged_ko = True

                                        dext = {}
                                        jext = {}
                                        vext = {}
                                        next = {}

                                        2: i *= n and (dext or jext)
                                        2: n *= i 
                                        2: j *= i and not (j and next) or jagged_oe and jagged_ko
                                        2: v *= vr and vext 
                                        2: d *= v or not i and not (d and next)
                                        2: vr *= v or not i

                                        3: tip *= d and vr and not i and not n and not j
                                        3: stalk *= i and n and j and not vr and not d
                                        3: tip_stalk *= d and vr and i and n and j
                                        3: hybrid *= not (tip or stalk or tip_stalk)
                                        """.format(n0, d0, j0, i0, vr0, v0, dext0, jext0, vext0, next0)
                                        for p in range(10):
                                            model = Model(text=text_ics, mode='async')
                                            model.initialize()
                                            model.iterate(steps=20)
                                            coll.collect(nodes = model.nodes, states = model.states)
                                        avg1 = coll.get_averages()

                                        on = ['jagged_oe', 'j']
                                        text_mod = boolean2.modify_states(text=text_ics, turnon=on)
                                        for r in range(10):
                                            model_jag_oe = Model(text=text_mod, mode='async')
                                            model_jag_oe.initialize()
                                            model_jag_oe.iterate(steps=20)
                                            coll.collect(nodes = model_jag_oe.nodes, states = model_jag_oe.states)
                                        avg2 = coll.get_averages()

                                        off = ['jagged_ko', 'j']
                                        text_mod2 = boolean2.modify_states(text=text_ics, turnoff=off)
                                        for t in range(10):
                                            model_jag_ko = Model(text=text_mod2, mode='async')
                                            model_jag_ko.initialize()
                                            model_jag_ko.iterate(steps=20)
                                            coll.collect(nodes = model_jag_ko.nodes, states = model_jag_ko.states)
                                        avg3 = coll.get_averages()

                                        n = avg1['n']
                                        d = avg1['d']
                                        j = avg1['j']
                                        i = avg1['i']
                                        vr = avg1['vr']
                                        v = avg1['v']
                                        t = avg1['tip']
                                        s = avg1['stalk']
                                        ts = avg1['tip_stalk']
                                        h = avg1['hybrid']

                                        n2 = avg2['n']
                                        d2 = avg2['d']
                                        j2 = avg2['j']
                                        i2 = avg2['i']
                                        vr2 = avg2['vr']
                                        v2 = avg2['v']
                                        t2 = avg2['tip']
                                        s2 = avg2['stalk']
                                        ts2 = avg2['tip_stalk']
                                        h2 = avg2['hybrid']

                                        n3 = avg3['n']
                                        d3 = avg3['d']
                                        j3 = avg3['j']
                                        i3 = avg3['i']
                                        vr3 = avg3['vr']
                                        v3 = avg3['v']
                                        t3 = avg3['tip']
                                        s3 = avg3['stalk']
                                        ts3 = avg3['tip_stalk']
                                        h3 = avg3['hybrid']

                                        """plot(d, 'b')
                                        plot(n, 'g')
                                        plot(j, 'r')
                                        plot(vr, 'm')
                                        plot(i, 'c')
                                        plot(v, 'k')
                                        show()"""

                                        x = [float(n[-1]), float(d[-1]), float(j[-1]), float(i[-1]),
                                             float(vr[-1]), float(v[-1])]
                                        y = [float(n[-2]), float(d[-2]), float(j[-2]), float(i[-2]),
                                             float(vr[-2]), float(v[-2])]
                                        x2 = [float(n2[-1]), float(d2[-1]), float(j2[-1]), float(i2[-1]),
                                              float(vr2[-1]), float(v2[-1])]
                                        y2 = [float(n2[-2]), float(d2[-2]), float(j2[-2]), float(i2[-2]),
                                              float(vr2[-2]), float(v2[-2])]
                                        x3 = [float(n3[-1]), float(d3[-1]), float(j3[-1]), float(i3[-1]),
                                              float(vr3[-1]), float(v3[-1])]
                                        y3 = [float(n3[-2]), float(d3[-2]), float(j3[-2]), float(i3[-2]),
                                              float(vr3[-2]), float(v3[-2])]
                                        ss = steady_state(x, y, 0.05)
                                        if ss:
                                            N = asynch_threshold(n[-1], threshold=0.40)
                                            D = asynch_threshold(d[-1], threshold=0.40)
                                            J = asynch_threshold(j[-1], threshold=0.40)
                                            I = asynch_threshold(i[-1], threshold=0.40)
                                            VR = asynch_threshold(vr[-1], threshold=0.40)
                                            V = asynch_threshold(v[-1], threshold=0.40)
                                            """T = asynch_threshold(tip[-1], threshold = 0.25)
                                            S = asynch_threshold(stalk[-1], threshold = 0.25)
                                            TS = asynch_threshold(ts[-1], threshold = 0.25)"""
                                            T = tip_f(icd=I, jagged=J, delta=D, vegfr=VR, notch=N)
                                            S = stalk_f(icd=I, jagged=J, delta=D, vegfr=VR, notch=N)
                                            TS = hybrid_f(delta=D, vegfr=VR, icd=I, jagged=J, notch=N)
                                            H = other(tip=T, hybrid=TS, stalk=S)

                                            T = tip_f(icd, jagged, delta, vegfr, notch)

                                            notch_ss.append(N)
                                            delta_ss.append(D)
                                            jagged_ss.append(J)
                                            icd_ss.append(I)
                                            vegfr_ss.append(VR)
                                            vegf_ss.append(V)
                                            tip_ss.append(T)
                                            stalk_ss.append(S)
                                            tip_stalk_ss.append(TS)
                                            hybrid_ss.append(H)
                                            dext_ics.append(dext0)
                                            jext_ics.append(jext0)
                                            next_ics.append(next0)
                                            vext_ics.append(vext0)
                                            jagged_ics.append(j0)
                                            ss_num += 1
                                            tip_j_bin.append(T)
                                            stalk_j_bin.append(S)
                                            hybrid_j_bin.append(TS)
                                            other_bin.append(H)
                                            if T == 1:
                                                total_tip += 1
                                            if S == 1:
                                                total_stalk += 1
                                            if H == 1:
                                                total_hybrid += 1
                                            if TS == 1:
                                                total_tip_stalk += 1
                                        else:
                                            print("ss not good")
                                        ss2 = steady_state(x2, y2, 0.05)
                                        if ss2:
                                            N2 = asynch_threshold(n[-1], threshold=0.40)
                                            D2 = asynch_threshold(d[-1], threshold=0.40)
                                            J2 = asynch_threshold(j[-1], threshold=0.40)
                                            I2 = asynch_threshold(i[-1], threshold=0.40)
                                            VR2 = asynch_threshold(vr[-1], threshold=0.40)
                                            V2 = asynch_threshold(v[-1], threshold=0.40)
                                            """T = asynch_threshold(tip[-1], threshold = 0.25)
                                            S = asynch_threshold(stalk[-1], threshold = 0.25)
                                            TS = asynch_threshold(ts[-1], threshold = 0.25)"""
                                            T2 = tip_f(icd=I2, jagged=J2, delta=D2, vegfr=VR2, notch=N2)
                                            S2 = stalk_f(icd=I2, jagged=J2, delta=D2, vegfr=VR2, notch=N2)
                                            TS2 = hybrid_f(delta=D2, vegfr=VR2, icd=I2, jagged=J2, notch=N2)
                                            H2 = other(tip=T2, hybrid=TS2, stalk=S2)

                                            notch_ss2.append(N2)
                                            delta_ss2.append(D2)
                                            jagged_ss2.append(J2)
                                            icd_ss2.append(I2)
                                            vegfr_ss2.append(VR2)
                                            vegf_ss2.append(V2)
                                            tip_ss2.append(T2)
                                            stalk_ss2.append(S2)
                                            tip_stalk_ss2.append(TS2)
                                            hybrid_ss2.append(H2)
                                            dext_ics2.append(dext0)
                                            jext_ics2.append(jext0)
                                            next_ics2.append(next0)
                                            vext_ics2.append(vext0)
                                            ss_num2 += 1
                                            if TS2 == 1:
                                                total_tip_j += 1
                                            if S2 == 1:
                                                total_stalk_j += 1
                                            if H2 == 1:
                                                total_hybrid_j += 1
                                            if TS2 == 1:
                                                total_tip_stalk_j += 1
                                        else:
                                            print("ss not good")
                                        ss3 = steady_state(x3, y3, 0.05)
                                        if ss3:
                                            N3 = asynch_threshold(n3[-1], threshold=0.40)
                                            D3 = asynch_threshold(d3[-1], threshold=0.40)
                                            J3 = asynch_threshold(j3[-1], threshold=0.40)
                                            I3 = asynch_threshold(i3[-1], threshold=0.40)
                                            VR3 = asynch_threshold(vr3[-1], threshold=0.40)
                                            V3 = asynch_threshold(v3[-1], threshold=0.40)
                                            """T = asynch_threshold(tip[-1], threshold = 0.25)
                                            S = asynch_threshold(stalk[-1], threshold = 0.25)
                                            TS = asynch_threshold(ts[-1], threshold = 0.25)"""
                                            T3 = tip_f(icd=I3, jagged=J3, delta=D3, vegfr=VR3, notch=N3)
                                            S3 = stalk_f(icd=I3, jagged=J3, delta=D3, vegfr=VR3, notch=N3)
                                            TS3 = hybrid_f(delta=D3, vegfr=VR3, icd=I3, jagged=J3, notch=N3)
                                            H3 = other(tip=T3, hybrid=TS3, stalk=S3)

                                            notch_ss3.append(N3)
                                            delta_ss3.append(D3)
                                            jagged_ss3.append(J3)
                                            icd_ss3.append(I3)
                                            vegfr_ss3.append(VR3)
                                            vegf_ss3.append(V3)
                                            tip_ss3.append(T3)
                                            stalk_ss3.append(S3)
                                            tip_stalk_ss3.append(TS3)
                                            hybrid_ss3.append(H3)
                                            dext_ics3.append(dext0)
                                            jext_ics3.append(jext0)
                                            next_ics3.append(next0)
                                            vext_ics3.append(vext0)
                                            ss_num2 += 1
                                            if TS3 == 1:
                                                total_tip_k += 1
                                            if S3 == 1:
                                                total_stalk_k += 1
                                            if H3 == 1:
                                                total_hybrid_k += 1
                                            if TS3 == 1:
                                                total_tip_stalk_k += 1
                                        else:
                                            print("ss not good")

print(
ss_num, "steady states achieved normal", ss_num2, "steady jag oe", ss_num3, "steady jagged ko", step, "total runs")
"""fig = plt.figure(figsize = (5, 11))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(n, 'g', label='notch')
ax1.plot(i, 'y', label='nicd')
ax1.plot(v, 'w', label = 'vegf')
plt.legend()
ax2.plot(d, 'b', label='delta')
ax2.plot(j, 'r', label = 'jagged')
ax2.plot(vr, 'm', label = 'vegfr')
plt.legend()
plt.show()"""

"""ics_jagged_ko = np.zeros((len(tip_j_bin), 10))
ics_jagged_oe = np.zeros((len(tip_j_bin), 10))
for k in range(len(tip_j_bin)):
    ics_jagged_oe[k] = [ics_n_j[k], ics_d_j[k],
                        ics_j_j[k], ics_i_j[k],
                        ics_vr_j[k], ics_v_j[k],
                        ics_next_j[k], ics_dext_j[k],
                        ics_jext_j[k], ics_vext_j[k]]

    ics_jagged_ko[k] = [ics_n[k], ics_d[k],
                        ics_j[k], ics_i[k],
                        ics_vr[k], ics_v[k],
                        ics_next[k], ics_dext[k],
                        ics_jext[k], ics_vext[k]]"""

import pandas as pd

data = {'Notch': notch_ss,
        'Delta': delta_ss,
        'Jagged': jagged_ss,
        'Icd': icd_ss,
        'Vegfr': vegfr_ss,
        'Vegf': vegf_ss,
        'Dext': dext_ics,
        'Jext': jext_ics,
        'Next': next_ics,
        'Vext': vext_ics,
        'Tip': tip_ss,
        'Stalk': stalk_ss,
        'Tip/Stalk': tip_stalk_ss,
        'Hybrid': hybrid_ss,
        'Jagged Expression': jagged_ics
        }

data2 = {'Notch': notch_ss2,
         'Delta': delta_ss2,
         'Jagged': jagged_ss2,
         'Icd': icd_ss2,
         'Vegfr': vegfr_ss2,
         'Vegf': vegf_ss2,
         'Dext': dext_ics2,
         'Jext': jext_ics2,
         'Next': next_ics2,
         'Vext': vext_ics2,
         'Tip': tip_ss2,
         'Stalk': stalk_ss2,
         'Tip/Stalk': tip_stalk_ss2,
         'Hybrid': hybrid_ss2}

data3 = {'Notch': notch_ss3,
         'Delta': delta_ss3,
         'Jagged': jagged_ss3,
         'Icd': icd_ss3,
         'Vegfr': vegfr_ss3,
         'Vegf': vegf_ss3,
         'Dext': dext_ics3,
         'Jext': jext_ics3,
         'Next': next_ics3,
         'Vext': vext_ics3,
         'Tip': tip_ss3,
         'Stalk': stalk_ss3,
         'Tip/Stalk': tip_stalk_ss3,
         'Hybrid': hybrid_ss3}

df = pd.DataFrame(data=data)
df2 = pd.DataFrame(data=data2)
df3 = pd.DataFrame(data=data3)

# df.to_csv('ndj.csv', sep = ",")
# dfp.to_csv('ndj_phenotype_transition_mat.csv', sep = ',')
"""dfj.to_csv('ndj_phenotype_jagged.csv', sep = ',')
dfr.to_csv('ndj_transition_rules.csv', sep=',')"""
# df_stalk.to_csv('Tip_no_stalk.csv', sep = ',')
# df_stalk_jagged.to_csv('Tip_w_stalk.csv', sep = ',')


import seaborn as sns;

sns.set()
import matplotlib.pyplot as plt

"""ax = sns.catplot(x = "Delta", y = "Notch", hue = "Jagged Expression", data = df)
ax.set(ylim = [-0.1, 1.1])
plt.show()
ax2 = sns.catplot(x = "Vegfr", y = "Icd", hue = "Jagged Expression", data = df)
ax2.set(ylim = [-0.1, 1.1])
plt.show()"""
"""fig, axs = plt.subplots(1,1)
fig.suptitle('Tip Adoption vs Jagged Mutation')
axs[0,0] = sns.catplot(y = "Tip", x = 'Dext',  hue = 'Vext', data = df)
axs[0,0].set(ylim = [-0.1, 1.1])
plt.show()"""
"""axs[1,0] = sns.catplot(y = "Tip", x = 'Dext',  hue = 'Vext', data = df2)
axs[1,0].set(ylim = [-0.1, 1.1])
axs[2,0] = sns.catplot(y = "Tip", x = 'Dext',  hue = 'Vext', data = df3)
axs[2,0].set(ylim = [-0.1, 1.1])
plt.show()"""


#---------------------------------
ax1 = sns.catplot(x = "Tip", y = "Stalk",  hue = "Jagged Expression", data = df)
ax1.set(ylim = [-0.1, 1.1])
ax1.set(title = 'Asynchronous No Jagged Mutation')
ax1.savefig('Asynch No Jagged Mutation Tip vs Stalk.pdf')
plt.show()
ax2 = sns.catplot(x = "Tip", y = "Stalk", data = df2)
ax2.set(ylim = [-0.1, 1.1])
ax2.set(title = 'Asynchronous Jagged Overexpressed')
ax2.savefig('Asynch Jagged Overexpressed Tip vs Stalk.pdf')
plt.show()
ax3 = sns.catplot(x = "Tip", y = "Stalk", data = df3)
ax3.set(ylim = [-0.1, 1.1])
ax3.set(title = 'Asynchronous Jagged Knockout')
ax3.savefig('Asynch Jagged Knockout Tip vs Stalk.pdf')
plt.show()
#------------------------------------


#------------------------------------
ax10 = sns.catplot(x = 'Dext', y = 'Tip', hue = 'Vext', data = df)
ax10.set(ylim = [-0.1, 1.1])
ax10.set(title = 'Asynchronous No Jagged Mutation')
ax10.savefig('Asynch No Jagged Mutation Tip vs Dext.pdf')
plt.show()
ax13 = sns.catplot(x = 'Dext', y = 'Tip', hue = 'Vext', data = df2)
ax13.set(ylim = [-0.1, 1.1])
ax13.set(title = 'Jagged Overrexpressed')
ax13.savefig('Asynch Jagged Overexpressed Tip vs Dext.pdf')
plt.show()
ax15 = sns.catplot(x = 'Dext', y = 'Tip', hue = 'Vext', data = df3)
ax15.set(ylim = [-0.1, 1.1])
ax15.set(title = 'Asynchronous Jagged Knockout')
ax15.savefig('Asynch Jagged Knockout Tiip vs Dext.pdf')
plt.show()
#-------------------------------------


#---------------------------------------
ax12 = sns.catplot(x = 'Dext', y = 'Tip/Stalk', hue = 'Vext', data = df)
ax12.set(ylim = [-0.1, 1.1])
ax12.set(title = 'Asynchronous No Jagged Mutation')
ax12.savefig('Asynch No Jagged Mutation Tip Stalk vs Dext.pdf')
plt.show()
ax18 = sns.catplot(x = 'Dext', y = 'Tip/Stalk', hue = 'Vext', data = df2)
ax18.set(ylim = [-0.1, 1.1])
ax18.set(title = 'Asynchronous Jagged Overexpressed')
ax18.savefig('Asynch Jagged Overexpressed Tip Stalk vs Dext.pdf')
plt.show()
ax17 = sns.catplot(x = 'Dext', y = 'Tip/Stalk', hue = 'Vext', data = df3)
ax17.set(ylim = [-0.1, 1.1])
ax17.set(title = 'Asynchronous Jagged Knockout')
ax17.savefig('Asynch Jagged Knockout Tip Stalk vs Dext.pdf')
plt.show()
#---------------------------------


#---------------------------------
ax11 = sns.catplot(x = 'Dext', y = 'Stalk', hue = 'Vext', data = df)
ax11.set(ylim = [-0.1, 1.1])
ax11.set(title = 'Asynchronous No Jagged Mutation')
ax11.savefig('Asynch No Jagged Mutation Stalk vs Dext.pdf')
plt.show()
ax14 = sns.catplot(x = 'Dext', y = 'Stalk', hue = 'Vext', data = df2)
ax14.set(ylim = [-0.1, 1.1])
ax14.set(title = 'Asynchronous Jagged Overexpressed')
ax14.savefig('Asynch Jagged Overexpressed Stalk vs Dext.pdf')
plt.show()
ax16 = sns.catplot(x = 'Dext', y = 'Stalk', hue = 'Vext', data = df3)
ax16.set(ylim = [-0.1, 1.1])
ax16.set(title = 'Asynchronous Jagged Knockout')
ax16.savefig('Asynch Jagged Knockout Stalk vs Dext.pdf')
plt.show()
#-------------------------------------


#---------------------------------
ax21 = sns.catplot(x = 'Dext', y = 'Hybrid', hue = 'Vext', data = df)
ax21.set(ylim = [-0.1, 1.1])
ax21.set(title = 'Asynchronous No Jagged Mutation')
ax21.savefig('Asynchn No Jagged Mutation Stalk vs Dext.pdf')
plt.show()
ax24 = sns.catplot(x = 'Dext', y = 'Hybrid', hue = 'Vext', data = df2)
ax24.set(ylim = [-0.1, 1.1])
ax24.set(title = 'Asynchronous Jagged Overexpressed')
ax24.savefig('Asynch Jagged Overexpressed Hybrid vs Dext.pdf')
plt.show()
ax26 = sns.catplot(x = 'Dext', y = 'Hybrid', hue = 'Vext', data = df3)
ax26.set(ylim = [-0.1, 1.1])
ax26.set(title = 'Asynchronous Jagged Knockout')
ax26.savefig('Asynch Jagged Knockout Hybrid vs Dext.pdf')
plt.show()
#-------------------------------------


"""ax14 = sns.catplot(x = 'dext', y = 'notch_ss', hue = "Jagged Initial", data = df)
ax14.set(ylim = [-0.1, 1.1])
ax15 = sns.catplot(x = 'jext', y = 'notch_ss', hue = "Jagged Initial", data = df)
ax15.set(ylim = [-0.1, 1.1])
ax16 = sns.catplot(x = 'next', y = 'notch_ss', hue = "Jagged Initial", data = df)
ax16.set(ylim = [-0.1, 1.1])
ax17 = sns.catplot(x = 'vext', y = 'notch_ss', hue = "Jagged Initial", data = df)
ax17.set(ylim = [-0.1, 1.1])

ax18 = sns.catplot(x = 'dext', y = 'jagged_ss', hue = "Jagged Initial", data = df)
ax18.set(ylim = [-0.1, 1.1])
ax19 = sns.catplot(x = 'jext', y = 'jagged_ss', hue = "Jagged Initial", data = df)
ax19.set(ylim = [-0.1, 1.1])
ax20 = sns.catplot(x = 'next', y = 'jagged_ss', hue = "Jagged Initial", data = df)
ax20.set(ylim = [-0.1, 1.1])
ax22 = sns.catplot(x = 'vext', y = 'jagged_ss', hue = "Jagged Initial", data = df)
ax22.set(ylim = [-0.1, 1.1])

ax23 = sns.catplot(x = 'dext', y = 'vegfr_ss', hue = "Jagged Initial", data = df)
ax23.set(ylim = [-0.1, 1.1])
ax24 = sns.catplot(x = 'jext', y = 'vegfr_ss', hue = "Jagged Initial", data = df)
ax24.set(ylim = [-0.1, 1.1])
ax25 = sns.catplot(x = 'next', y = 'vegfr_ss', hue = "Jagged Initial", data = df)
ax25.set(ylim = [-0.1, 1.1])
ax26 = sns.catplot(x = 'vext', y = 'vegfr_ss', hue = "Jagged Initial", data = df)
ax26.set(ylim = [-0.1, 1.1])

ax33 = sns.catplot(x = 'dext', y = 'vegf_ss', hue = "Jagged Initial", data = df)
ax33.set(ylim = [-0.1, 1.1])
ax34 = sns.catplot(x = 'jext', y = 'vegf_ss', hue = "Jagged Initial", data = df)
ax34.set(ylim = [-0.1, 1.1])
ax35 = sns.catplot(x = 'next', y = 'vegf_ss', hue = "Jagged Initial", data = df)
ax35.set(ylim = [-0.1, 1.1])
ax36 = sns.catplot(x = 'vext', y = 'vegf_ss', hue = "Jagged Initial", data = df)
ax36.set(ylim = [-0.1, 1.1])"""

#ax37 = sns.catplot(x='Tip', y='Stalk', data=dfj)
#ax37.set(ylim=[-0.1,1.1])


"""ax34 = sns.catplot(x='tip', y='stalk',  hue = 'Jagged Initial', col='dext', data=df)
ax34.set(ylim=[-0.1,1.1])

ax38 = sns.catplot(x='tip', y='hybrid',  hue = 'Jagged Initial', col='dext', data=df)
ax38.set(ylim=[-0.1,1.1])

ax39 = sns.catplot(x='stalk', y='hybrid',  hue = 'Jagged Initial', col='dext', data=df)
ax39.set(ylim=[-0.1,1.1])

ax40 = sns.catplot(x='tip/stalk', y='tip',  hue = 'Jagged Initial', col='dext', data=df)
ax40.set(ylim=[-0.1,1.1])

ax41 = sns.catplot(x='tip/stalk', y='stalk',  hue = 'Jagged Initial', col='dext', data=df)
ax41.set(ylim=[-0.1,1.1])

ax99 = sns.catplot(x='hybrid', y = 'tip', hue = 'Jagged Initial', col = 'dext', data = df)
plt.show()"""

tip_prop = total_tip * 1.0 / ss_num
stalk_prop = total_stalk * 1.0 / ss_num
tip_stalk_prop = total_tip_stalk * 1.0 / ss_num
hybrid_prop = total_hybrid * 1.0 / ss_num

tip_prop_j = total_tip_j * 1.0 / ss_num2
stalk_prop_j = total_stalk_j * 1.0 / ss_num2
tip_stalk_prop_j = total_tip_stalk_j * 1.0 / ss_num2
hybrid_prop_j = total_hybrid_j * 1.0 / ss_num2

tip_prop_k = total_tip_k * 1.0 / ss_num3
stalk_prop_k = total_stalk_k * 1.0 / ss_num3
tip_stalk_prop_k = total_tip_stalk_k * 1.0 / ss_num3
hybrid_prop_k = total_hybrid_k * 1.0 / ss_num3



# lets look at the distribution of the phenotypes according to their  distribution
label = ['Tip', 'Stalk', 'Tip/Stalk', 'Hybrid']
jag = [round(tip_prop, 2), round(stalk_prop, 2), round(tip_stalk_prop, 2), round(hybrid_prop, 2)]
jag_oe = [round(tip_prop_j, 2), round(stalk_prop_j, 2), round(tip_stalk_prop_j, 2), round(hybrid_prop_j, 2)]
jag_ko = [round(tip_prop_k, 2), round(stalk_prop_k, 2), round(tip_stalk_prop_k, 2), round(hybrid_prop_k, 2)]

x = np.arange(len(label))
width = 0.15
fig, ax5 = plt.subplots()

rects1 = ax5.bar(x, jag, width, label='No Jagged Mutation')
rects2 = ax5.bar(x + width, jag_oe, width, label='Jagged o/e')
rects3 = ax5.bar(x - width, jag_ko, width, label='Jagged k/o')


ax5.set_ylabel('Percentage Phenotype Adoption')
ax5.set_title('Asynchronous Phenotype by Jagged Mutation')
ax5.set_xticks(x)
ax5.set_xticklabels(label)
ax5.legend()


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax5.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 3, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


autolabel(rects1)
autolabel(rects2)
autolabel(rects3)

fig.tight_layout()

plt.show()
ax5.savefig('Asynchronous Proportion of Phenotype Adoption vs Jagged Mutation.pdf')


print("total tips:", total_tip, total_tip_j, total_tip_k, "stalk:", total_stalk, total_stalk_j, total_stalk_k,
      "total tip/stalk", total_tip_stalk, total_tip_stalk_j, total_tip_stalk_k, "hybrid", total_hybrid, total_hybrid_j, total_hybrid_k)
print("total steady states:", ss_num, ss_num2, ss_num3)