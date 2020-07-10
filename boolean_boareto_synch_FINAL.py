# testing new rules for boolean boareto
import boolean2
from boolean2 import Model, util
#import matplotlib; matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
import random
from pylab import *



#2: n *= i or not (n and dext) or not (n and jext) or not (n and j) or not (n and d)
#2: n *= i or not (n and (d or j))


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
    text_new = boolean2.modify_states(text = text, turnoff = off, turnon = on)
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

def asynch_threshold(protein, threshold = 0.30):
    """ asynch determine if on or off """
    if protein >= threshold:
        return 1
    else:
        return 0

def asynch_threshold_list(proteins, threshold = 0.25):
    bin = []
    if len(proteins) > 1:
        for i in range(len(proteins)):
            if proteins[i] >= threshold:
                bin.append(1)
            else:
                bin.append(0)
    return bin

def tip_f(icd, jagged, delta, vegfr, notch):
    if delta == 1 and vegfr == 1  and jagged == 0 and icd == 0 and notch == 0:
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
#while flag == True:

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

d_ics2 = []
j_ics2 = []
n_ics2 = []
i_ics2 = []
vr_ics2 = []
v_ics2 = []
dext_ics2 = []
jext_ics2 = []
next_ics2 = []
vext_ics2 = []

d_ics3 = []
j_ics3 = []
n_ics3 = []
i_ics3 = []
vr_ics3 = []
v_ics3 = []
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


n_bin =[]
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
R = 2**(len(ics_bin))
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

ics_bin = []
tip_ics_bin = []
stalk_ics_bin = []
tipstalk_ics_bin = []
hybrid_ics_bin = []

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
                                        ics = []
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
                                        model = Model(text = text_ics, mode = 'sync')
                                        on = ['jagged_oe', 'j']
                                        off = ['jagged_ko', 'j']
                                        text_mod = boolean2.modify_states(text = text_ics, turnon = on)
                                        model_jag_oe = Model(text = text_mod, mode = 'sync')
                                        text_mod2 = boolean2.modify_states(text = text_ics, turnoff = off)
                                        model_jag_ko = Model(text = text_mod2, mode = 'sync')
                                        model.initialize()
                                        model_jag_oe.initialize()
                                        model_jag_ko.initialize()
                                        model.iterate(steps = 20)
                                        model_jag_oe.iterate(steps = 20)
                                        model_jag_ko.iterate(steps = 20)

                                        n = model.data['n']
                                        d = model.data['d']
                                        j = model.data['j']
                                        i = model.data['i']
                                        vr = model.data['vr']
                                        v = model.data['v']
                                        t = model.data['tip']
                                        s = model.data['stalk']
                                        ts = model.data['tip_stalk']
                                        h = model.data['hybrid']

                                        n2 = model_jag_oe.data['n']
                                        d2 = model_jag_oe.data['d']
                                        j2 = model_jag_oe.data['j']
                                        i2 = model_jag_oe.data['i']
                                        vr2 = model_jag_oe.data['vr']
                                        v2 = model_jag_oe.data['v']
                                        t2 = model_jag_oe.data['tip']
                                        s2 = model_jag_oe.data['stalk']
                                        ts2 = model_jag_oe.data['tip_stalk']
                                        h2 = model_jag_oe.data['hybrid']

                                        n3 = model_jag_ko.data['n']
                                        d3 = model_jag_ko.data['d']
                                        j3 = model_jag_ko.data['j']
                                        i3 = model_jag_ko.data['i']
                                        vr3 = model_jag_ko.data['vr']
                                        v3 = model_jag_ko.data['v']
                                        t3 = model_jag_ko.data['tip']
                                        s3 = model_jag_ko.data['stalk']
                                        ts3 = model_jag_ko.data['tip_stalk']
                                        h3 = model_jag_ko.data['hybrid']

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
                                            ics = [notch, delta, jagged, icd, vegfr, vegf, d_ext, j_ext, n_ext, v_ext]
                                            notch_ss.append(n[-1])
                                            delta_ss.append(d[-1])
                                            jagged_ss.append(j[-1])
                                            icd_ss.append(i[-1])
                                            vegfr_ss.append(vr[-1])
                                            vegf_ss.append(v[-1])
                                            tip_ss.append(t[-1])
                                            stalk_ss.append(s[-1])
                                            tip_stalk_ss.append(ts[-1])
                                            hybrid_ss.append(h[-1])
                                            n_ics.append(notch)
                                            d_ics.append(delta)
                                            j_ics.append(jagged)
                                            i_ics.append(icd)
                                            vr_ics.append(vegfr)
                                            v_ics.append(vegf)
                                            dext_ics.append(d_ext)
                                            jext_ics.append(j_ext)
                                            next_ics.append(n_ext)
                                            vext_ics.append(v_ext)
                                            jagged_ics.append(j0)
                                            ss_num += 1
                                            tip_j_bin.append(tip_ss[-1])
                                            stalk_j_bin.append(stalk_ss[-1])
                                            hybrid_j_bin.append(tip_stalk_ss[-1])
                                            n_ss_j.append(n[-1])
                                            ics_n_j.append(notch)
                                            d_ss_j.append(d[-1])
                                            ics_d_j.append(delta)
                                            j_ss_j.append(j[-1])
                                            ics_j_j.append(jagged)
                                            i_ss_j.append(i[-1])
                                            ics_i_j.append(icd)
                                            vr_ss_j.append(vr[-1])
                                            ics_vr_j.append(vegfr)
                                            v_ss_j.append(v[-1])
                                            ics_v_j.append(vegf)
                                            ics_dext_j.append(d_ext)
                                            ics_jext_j.append(j_ext)
                                            ics_vext_j.append(v_ext)
                                            ics_next_j.append(n_ext)
                                            other_bin.append(h[-1])
                                            if tip_ss[-1] == 1:
                                                total_tip += 1
                                                tip_ics_bin.append(ics)
                                            if stalk_ss[-1] == 1:
                                                total_stalk += 1
                                                stalk_ics_bin.append(ics)
                                            if hybrid_ss[-1] == 1:
                                                total_hybrid += 1
                                                hybrid_ics_bin.append(ics)
                                            if tip_stalk_ss[-1]:
                                                total_tip_stalk += 1
                                                tipstalk_ics_bin.append(ics)
                                        else:
                                            print("ss not good")
                                        ss2 = steady_state(x2, y2, 0.05)
                                        if ss2:
                                            notch_ss2.append(n2[-1])
                                            delta_ss2.append(d2[-1])
                                            jagged_ss2.append(j2[-1])
                                            icd_ss2.append(i2[-1])
                                            vegfr_ss2.append(v2[-1])
                                            vegf_ss2.append(v2[-1])
                                            tip_ss2.append(t2[-1])
                                            stalk_ss2.append(s2[-1])
                                            tip_stalk_ss2.append(ts2[-1])
                                            hybrid_ss2.append(h2[-1])
                                            n_ics2.append(notch)
                                            d_ics2.append(delta)
                                            j_ics2.append(jagged)
                                            i_ics2.append(icd)
                                            vr_ics2.append(vegfr)
                                            v_ics2.append(vegf)
                                            dext_ics2.append(d_ext)
                                            jext_ics2.append(j_ext)
                                            next_ics2.append(n_ext)
                                            vext_ics2.append(v_ext)
                                            #jagged_ics.append(j0)
                                            ss_num2 += 1
                                            if tip_ss2[-1] == 1:
                                                total_tip_j += 1
                                            if stalk_ss2[-1] == 1:
                                                total_stalk_j += 1
                                            if hybrid_ss2[-1] == 1:
                                                total_hybrid_j += 1
                                            if tip_stalk_ss2[-1]:
                                                total_tip_stalk_j += 1
                                        else:
                                            print("ss not good")
                                        ss3 = steady_state(x3, y3, 0.05)
                                        if ss3:
                                            notch_ss3.append(n3[-1])
                                            delta_ss3.append(d3[-1])
                                            jagged_ss3.append(j3[-1])
                                            icd_ss3.append(i3[-1])
                                            vegfr_ss3.append(v3[-1])
                                            vegf_ss3.append(v3[-1])
                                            tip_ss3.append(t3[-1])
                                            stalk_ss3.append(s3[-1])
                                            tip_stalk_ss3.append(ts3[-1])
                                            hybrid_ss3.append(h3[-1])
                                            n_ics3.append(notch)
                                            d_ics3.append(delta)
                                            j_ics3.append(jagged)
                                            i_ics3.append(icd)
                                            vr_ics3.append(vegfr)
                                            v_ics3.append(vegf)
                                            dext_ics3.append(d_ext)
                                            jext_ics3.append(j_ext)
                                            next_ics3.append(n_ext)
                                            vext_ics3.append(v_ext)
                                            #jagged_ics.append(j0)
                                            ss_num3 += 1
                                            if tip_ss3[-1] == 1:
                                                total_tip_k += 1
                                            if stalk_ss3[-1] == 1:
                                                total_stalk_k += 1
                                            if hybrid_ss3[-1] == 1:
                                                total_hybrid_k += 1
                                            if tip_stalk_ss3[-1]:
                                                total_tip_stalk_k += 1
                                        else:
                                            print("ss not good")

print(ss_num, "steady states achieved normal", ss_num2, "steady jag oe", ss_num3, "steady jagged ko", step, "total runs")
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
        "N0": n_ics,
        'D0': d_ics,
        'J0': j_ics,
        'I0': i_ics,
        'VR0': vr_ics,
        'V0': v_ics,
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
         "N0": n_ics2,
         'D0': d_ics2,
         'J0': j_ics2,
         'I0': i_ics2,
         'VR0': vr_ics2,
         'V0': v_ics2,
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
         "N0": n_ics3,
         'D0': d_ics3,
         'J0': j_ics3,
         'I0': i_ics3,
         'VR0': vr_ics3,
         'V0': v_ics3,
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

#df.to_csv('ndj.csv', sep = ",")
#dfp.to_csv('ndj_phenotype_transition_mat.csv', sep = ',')
"""dfj.to_csv('ndj_phenotype_jagged.csv', sep = ',')
dfr.to_csv('ndj_transition_rules.csv', sep=',')"""
#df_stalk.to_csv('Tip_no_stalk.csv', sep = ',')
#df_stalk_jagged.to_csv('Tip_w_stalk.csv', sep = ',')


import seaborn as sns; sns.set()
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


data_ics_tip = {'Tip': tip_ics_bin}
data_ics_stalk = {'Stalk': stalk_ics_bin}
data_ics_tipstalk = {'Tip/Stalk': tipstalk_ics_bin}
data_ics_hybrid = {'Hybrid': hybrid_ics_bin}

df_ics_tip = pd.DataFrame(data_ics_tip)
df_ics_tip.to_csv('Tip Intitial Conditions Synchronous Bin', sep = ',')

df_ics_stalk = pd.DataFrame(data_ics_stalk)
df_ics_stalk.to_csv('Stalk Intitial Conditions Synchronous Bin', sep = ',')

df_ics_hybrid = pd.DataFrame(data_ics_hybrid)
df_ics_hybrid.to_csv('Hybrid Intitial Conditions Synchronous Bin', sep = ',')

df_ics_tipstalk = pd.DataFrame(data_ics_tipstalk)
df_ics_tipstalk.to_csv('Tip_Stalk Intitial Conditions Synchronous Bin', sep = ',')

df.to_csv('Steady States Synchronous Bin')

#----------------------
# new corrected version for obtaining all the data one time
jagged_roam_free = pd.DataFrame(df)
jagged_oe = pd.DataFrame(df2)
jagged_ko = pd.DataFrame(df3)

jagged_roam_free.to_csv('jagged roam free', sep = ',')
jagged_oe.to_csv('jagged oe', sep = ',')
jagged_ko.to_csv('jagged ko', sep = ',')

#---------------------------------
ax1 = sns.catplot(x = "Tip", y = "Stalk",  hue = "Jagged Expression", data = df)
ax1.set(ylim = [-0.1, 1.1])
ax1.set(title = 'No Jagged Mutation')
ax1.savefig('No Jagged Mutation Tip vs Stalk.pdf')
plt.show()
ax2 = sns.catplot(x = "Tip", y = "Stalk", data = df2)
ax2.set(ylim = [-0.1, 1.1])
ax2.set(title = 'Jagged Overexpressed')
ax2.savefig('Jagged Overexpressed Tip vs Stalk.pdf')
plt.show()
ax3 = sns.catplot(x = "Tip", y = "Stalk", data = df3)
ax3.set(ylim = [-0.1, 1.1])
ax3.set(title = 'Jagged Knockout')
ax3.savefig('Jagged Knockout Tip vs Stalk.pdf')
plt.show()
#------------------------------------


#------------------------------------
ax10 = sns.catplot(x = 'Dext', y = 'Tip', hue = 'Vext', data = df)
ax10.set(ylim = [-0.1, 1.1])
ax10.set(title = 'No Jagged Mutation')
ax10.savefig('No Jagged Mutation Tip vs Dext.pdf')
plt.show()
ax13 = sns.catplot(x = 'Dext', y = 'Tip', hue = 'Vext', data = df2)
ax13.set(ylim = [-0.1, 1.1])
ax13.set(title = 'Jagged Overrexpressed')
ax13.savefig('Jagged Overexpressed Tip vs Dext.pdf')
plt.show()
ax15 = sns.catplot(x = 'Dext', y = 'Tip', hue = 'Vext', data = df3)
ax15.set(ylim = [-0.1, 1.1])
ax15.set(title = 'Jagged Knockout')
ax15.savefig('Jagged Knockout Tiip vs Dext.pdf')
plt.show()
#-------------------------------------


#---------------------------------------
ax12 = sns.catplot(x = 'Dext', y = 'Tip/Stalk', hue = 'Vext', data = df)
ax12.set(ylim = [-0.1, 1.1])
ax12.set(title = 'No Jagged Mutation')
ax12.savefig('No Jagged Mutation Tip Stalk vs Dext.pdf')
plt.show()
ax18 = sns.catplot(x = 'Dext', y = 'Tip/Stalk', hue = 'Vext', data = df2)
ax18.set(ylim = [-0.1, 1.1])
ax18.set(title = 'Jagged Overexpressed')
ax18.savefig('Jagged Overexpressed Tip Stalk vs Dext.pdf')
plt.show()
ax17 = sns.catplot(x = 'Dext', y = 'Tip/Stalk', hue = 'Vext', data = df3)
ax17.set(ylim = [-0.1, 1.1])
ax17.set(title = 'Jagged Knockout')
ax17.savefig('Jagged Knockout Tip Stalk vs Dext.pdf')
plt.show()
#---------------------------------


#---------------------------------
ax11 = sns.catplot(x = 'Dext', y = 'Stalk', hue = 'Vext', data = df)
ax11.set(ylim = [-0.1, 1.1])
ax11.set(title = 'No Jagged Mutation')
ax11.savefig('No Jagged Mutation Stalk vs Dext.pdf')
plt.show()
ax14 = sns.catplot(x = 'Dext', y = 'Stalk', hue = 'Vext', data = df2)
ax14.set(ylim = [-0.1, 1.1])
ax14.set(title = 'Jagged Overexpressed')
ax14.savefig('Jagged Overexpressed Stalk vs Dext.pdf')
plt.show()
ax16 = sns.catplot(x = 'Dext', y = 'Stalk', hue = 'Vext', data = df3)
ax16.set(ylim = [-0.1, 1.1])
ax16.set(title = 'Jagged Knockout')
ax16.savefig('Jagged Knockout Stalk vs Dext.pdf')
plt.show()
#-------------------------------------


#---------------------------------
ax21 = sns.catplot(x = 'Dext', y = 'Hybrid', hue = 'Vext', data = df)
ax21.set(ylim = [-0.1, 1.1])
ax21.set(title = 'No Jagged Mutation')
ax21.savefig('No Jagged Mutation Stalk vs Dext.pdf')
plt.show()
ax24 = sns.catplot(x = 'Dext', y = 'Hybrid', hue = 'Vext', data = df2)
ax24.set(ylim = [-0.1, 1.1])
ax24.set(title = 'Jagged Overexpressed')
ax24.savefig('Jagged Overexpressed Hybrid vs Dext.pdf')
plt.show()
ax26 = sns.catplot(x = 'Dext', y = 'Hybrid', hue = 'Vext', data = df3)
ax26.set(ylim = [-0.1, 1.1])
ax26.set(title = 'Jagged Knockout')
ax26.savefig('Jagged Knockout Hybrid vs Dext.pdf')
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
ax5.set_title('Phenotype by Jagged Mutation')
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
ax5.savefig('Proportion of Phenotype Adoption vs Jagged Mutation.pdf')


print("total tips:", total_tip, total_tip_j, total_tip_k, "stalk:", total_stalk, total_stalk_j, total_stalk_k,
      "total tip/stalk", total_tip_stalk, total_tip_stalk_j, total_tip_stalk_k, "hybrid", total_hybrid, total_hybrid_j, total_hybrid_k)
print("total steady states:", ss_num, ss_num2, ss_num3)