import graphviz

from graphviz import Graph, Digraph
import numpy as np
import pandas as pd
import scipy
import random


def hamming_distance(string_1, string_2):
    ham_dist = 0
    # zip takes a string or input and returns a tuple the
    # size of the number of characters
    # if two strings or input of different lenghts then
    # it will only count the first k in common
    # where k is the smallest size of the tuple
    #assert len(string_1) == len(string_2)
    for ch1, ch2 in zip(string_1, string_2):
        if ch1 != ch2:
            ham_dist += 1
    return ham_dist


"""f = Digraph("Boareto's Framework", engine = 'dot')#, filename = 'BF.gv')

f.attr('node', shape = 'doublecircle')
f.node('Delta')
f.node('Jagged')
f.node('Notch')
f.node('Vegfr')
f.node('Vegf')
f.node('Nicd')

f.attr('node', shape = 'square')
f.node('Dext')
f.node('Vext')
f.node('Jext')
f.node('Next')

f.attr('node', shape = 'circle')
f.edge('Vext', 'Vegfr', label = 'Interaction')
f.edge('Vegfr', 'Vegf', label = 'Activate')
f.edge('Vegf', 'Vegfr', label = 'Activate')
f.edge('Next', 'Jagged', label = 'Interaction')
f.edge('Next', 'Delta', label = 'Interaction')
f.edge('Dext', 'Notch', label = 'Cleavage')
f.edge('Jext', 'Notch', label = 'Cleavage')
f.edge('Cleavage', 'Trans-Interaction', label = 'Activate')
f.edge('Trans-Interaction', 'Nicd', label = 'Activate')
f.edge('Nicd', 'Notch', label = 'Activate')
f.edge('Nicd', 'Vegfr', label = 'Inhibit')
f.edge('Nicd', 'Delta', label = 'Inhibit')
f.edge('Nicd', 'Jagged', label = 'Activate')

f.edge('Vegf', 'Delta', label = 'Activate')

f.view()"""


f = Digraph("Boareto's Framework", engine = 'dot')

f.attr('node', shape = 'doublecircle')
f.node('D')
f.node('J')
f.node('N')
f.node('Vr')
f.node('V')
f.node('I')

f.attr('node', shape = 'square')
f.node('Dext')
f.node('Vext')
f.node('Jext')
f.node('Next')

f.attr('node', shape = 'circle')
f.edge('Vext', 'Vr', label = 'Interaction')
f.edge('Vr', 'V', label = 'Activate')
f.edge('V', 'D', label = 'Activate')
f.edge('Next', 'J', label = 'Interaction')
f.edge('Next', 'D', label = 'Interaction')
f.edge('Dext', 'N', label = 'Cleavage')
f.edge('Jext', 'N', label = 'Cleavage')
f.edge('Cleavage', 'I', label = 'Activate')
f.edge('I', 'N', label = 'Activate')
f.edge('I', 'Vr', label = 'Inhibit')
f.edge('I', 'D', label = 'Inhibit')
f.edge('I', 'J', label = 'Activate')

#f.view()


tips = pd.read_csv('Tip Intitial Conditions Synchronous Bin', usecols = [1])
stalks = pd.read_csv('Stalk Intitial Conditions Synchronous Bin', usecols = [1])
hybrids = pd.read_csv('Hybrid Intitial Conditions Synchronous Bin', usecols = [1])
tipstalks = pd.read_csv('Tip_Stalk Intitial Conditions Synchronous Bin', usecols = [1])

tip_ics = tips.to_numpy()
stalk_ics = stalks.to_numpy()
hybrid_ics = hybrids.to_numpy()
tip_stalk_ics = tipstalks.to_numpy()

"""tip_ics = tips.get_values()
stalk_ics = stalks.get_values()
hybrid_ics = hybrids.get_values()
tip_stalk_ics = tipstalks.get_values()"""
print(tip_ics)

for k in range(len(tip_ics)):
    val = tip_ics[k]
    print("value:", val)
    print(val[0])


    def reorder(ics_bin):
        """A function that reorders all the intiial conditons for
        a phenotype where each ics is seperated by the minimum hamming distance"""
        N = len(ics_bin)
        for k in range(N):
            rowplace = ics_bin[k]
            row = rowplace[0]
            # place holder value:
            val = 0
            # hamming_array collects all the hamming distances
            # between the kth row and the rows beneath
            # we wish to find the smallest hamming distance to the kth row
            # this will be the neighbor to that row
            # calculated hamming distance between for all the rows beneath it
            hamming_array = []
            for j in range(k, N):
                nrowplace = ics_bin[j]
                nrow = nrowplace[0]
                nval = hamming_distance(row, nrow)
                hamming_array.append(nval)
            minimum, argmin = np.min(hamming_array), np.argmin(hamming_array)
            # swap the smallest with the next one
            if k < N - 1:
                ics_bin[[k + 1, argmin]] = ics_bin[[argmin, k + 1]]
        return ics_bin


# ics = [notch, delta, jagged, icd, vegfr, vegf, d_ext, j_ext, n_ext, v_ext] -- bar code definition
g = Digraph('Phenotype tip')# node_attr={'color': 'red', 'style': 'filled'})
g.attr('node', shape='square')
g.attr(size = '2,2')


K = len(tip_ics)
K2 = len(stalk_ics)
K3 = len(hybrid_ics)
K4 = len(tip_stalk_ics)
v = int(3 * K/4)
v2 = int(3 * K2/4)
v3 = int(3 * K3/4)
v4 = int(3 * K4/4)

""" graph all the steady states: """

""" do the tip first """
ics_bin = tip_ics[v:]
N = len(ics_bin)
choice_list = np.arange(1,8,1)
for k in range(N):
    rowplace = ics_bin[k]
    row = rowplace[0]
    row_str = "{}".format(row)
    row_str = row_str.replace(', ', '').replace('[','').replace(']','')
    hamming_array = []
    randval = random.choice(choice_list)
    for j in range(k, N):
        nrowplace = ics_bin[j]
        nrow = nrowplace[0]
        nrow_str = "{}".format(nrow)
        nrow_str = nrow_str.replace(', ', '').replace('[', '').replace(']', '')
        nval = hamming_distance(row, nrow)
        # here we will look at neighbor states and create an arrow
        if nval == 1 and randval == 3:
            g.edge(row_str, nrow_str)
            g.attr(label=r'\n\nTip Stable State Initial Conditions')
            g.attr(fontsize='50')

g.view()

r = Digraph('Phenotype stalk')# node_attr={'color': 'red', 'style': 'filled'})
r.attr('node', shape='square')
r.attr(size = '2,2')


""" next the stalk """
""" do the tip first """
ics_bin = stalk_ics[v2:]
N = len(ics_bin)
choice_list = np.arange(1,8,1)
for k in range(N):
    rowplace = ics_bin[k]
    row = rowplace[0]
    row_str = "{}".format(row)
    row_str = row_str.replace(', ', '').replace('[','').replace(']','')
    hamming_array = []
    randval = random.choice(choice_list)
    for j in range(k, N):
        nrowplace = ics_bin[j]
        nrow = nrowplace[0]
        nrow_str = "{}".format(nrow)
        nrow_str = nrow_str.replace(', ', '').replace('[', '').replace(']', '')
        nval = hamming_distance(row, nrow)
        # here we will look at neighbor states and create an arrow
        if nval == 1: #and randval == 3:
            r.edge(row_str, nrow_str)
            r.attr(label=r'\n\Stalk Stable State Initial Conditions')
            r.attr(fontsize='30')

#r.view()

s = Digraph('Phenotype tip_stalk')# node_attr={'color': 'red', 'style': 'filled'})
s.attr('node', shape='square')
s.attr(size = '2,2')


""" do the tip/stalk """
ics_bin = tip_stalk_ics[v4:]
N = len(ics_bin)
choice_list = np.arange(1,5,1)
for k in range(N):
    rowplace = ics_bin[k]
    row = rowplace[0]
    row_str = "{}".format(row)
    row_str = row_str.replace(', ', '').replace('[','').replace(']','')
    hamming_array = []
    randval = random.choice(choice_list)
    for j in range(k, N):
        nrowplace = ics_bin[j]
        nrow = nrowplace[0]
        nrow_str = "{}".format(nrow)
        nrow_str = nrow_str.replace(', ', '').replace('[', '').replace(']', '')
        nval = hamming_distance(row, nrow)
        # here we will look at neighbor states and create an arrow
        if nval == 1:# and randval == 3:
            s.edge(row_str, nrow_str)
            s.attr(label=r'\n\nTip/Stalk Stable State Initial Conditions')
            s.attr(fontsize='20')

#s.view()

g = Digraph('G', filename='cluster.gv')

# NOTE: the subgraph name needs to begin with 'cluster' (all lowercase)
#       so that Graphviz recognizes it as a special cluster subgraph



with g.subgraph(name='cluster_0') as c:
    c.attr(style='filled', color='lightgrey')
    c.node_attr.update(style='filled', color='white')
    c.attr(label='process #1')

with g.subgraph(name='cluster_1') as c:
    c.attr(color='blue')
    c.node_attr['style'] = 'filled'
    c.edges([('b0', 'b1'), ('b1', 'b2'), ('b2', 'b3')])
    c.attr(label='process #2')

g.edge('start', 'a0')
g.edge('start', 'b0')
g.edge('a1', 'b3')
g.edge('b2', 'a3')
g.edge('a3', 'a0')
g.edge('a3', 'end')
g.edge('b3', 'end')

g.node('start', shape='Mdiamond')
g.node('end', shape='Msquare')

#g.view()

t = Digraph('G', filename='tip stalk crosstalk.gv')
t.attr('node', shape='square')
t.attr(size = '2,2')


" lets look at tip stalk crosstalk:"
ics_bin1, ics_bin2 = tip_ics[30:50], stalk_ics[30:50]
N = len(ics_bin1)
choice_list = np.arange(1,5,1)


tip_arrow = []
stalk_arrow = []
for k in range(N):
    tip_rowplace = ics_bin1[k]
    stalk_rowplace = ics_bin2[k]
    tip_row = tip_rowplace[0]
    stalk_row = stalk_rowplace[0]

    tip_row_str = "{}".format(tip_row)
    tip_row_str = tip_row_str.replace(', ', '').replace('[', '').replace(']', '')
    stalk_row_str = "{}".format(stalk_row)
    stalk_row_str = stalk_row_str.replace(', ', '').replace('[', '').replace(']', '')
    randval = random.choice(choice_list)
    for j in range(k, N):
        tip_nrowplace = ics_bin1[j]
        stalk_nrowplace = ics_bin2[j]
        tip_nrow = tip_nrowplace[0]
        stalk_nrow = stalk_nrowplace[0]
        tip_nrow_str = "{}".format(tip_nrow)
        tip_nrow_str = tip_nrow_str.replace(', ', '').replace('[', '').replace(']', '')
        tip_nval = hamming_distance(tip_row, tip_nrow)
        stalk_nrow_str = "{}".format(stalk_nrow)
        stalk_nrow_str = stalk_nrow_str.replace(', ', '').replace('[', '').replace(']', '')
        stalk_nval = hamming_distance(stalk_row, stalk_nrow)
        tip_stalk_nval = hamming_distance(tip_row, stalk_nrow)
        stalk_tip_nval = hamming_distance(stalk_row, tip_nrow)
        # here we will look at neighbor states and create an arrow
        if tip_nval == 1:# and randval == 3:
            t.attr(color='red')
            string_tip = "({},{})".format(tip_row_str, tip_nrow_str)
            t.edge(tip_row_str, tip_nrow_str)
            print(string_tip)
            tip_arrow.append(string_tip)
            t.attr(label=r'\n\nTip and Stalk Stable State Initial Conditions')
            t.attr(fontsize='50')
        if stalk_nval == 1:
            t.attr(color='blue')
            string_stalk = "({},{})".format(stalk_row_str, stalk_nrow_str)
            stalk_arrow.append(string_stalk)
            t.edge(stalk_row_str, stalk_nrow_str)
        if tip_stalk_nval == 1:
            t.edge(tip_row_str, stalk_nrow_str, label = 'ts')
        if stalk_tip_nval == 1:
            t.edge(stalk_row_str, tip_nrow_str, label = 'ts')
"""""
with t.subgraph(name='tip cluster') as c:
    c.attr(label='Tip')
    c.attr('node', shape='square')
    c.edges(tip_arrow)
with t.subgraph(name='stalk cluster') as c:
    c.attr(label='Stalk')
    c.attr('node', shape='square')
    c.edges(stalk_arrow)
"""""




t.view()



#---------------------------------------------------------
# look at all the tip states and stalk states




t = Digraph('G', filename='tip and stalk .gv')
t.attr('node', shape='square')
t.attr(size = '3,3')


" lets look at tip stalk crosstalk:"
ics_bin1, ics_bin2 = tip_ics[:], stalk_ics[:]
N = len(ics_bin1)
Z = len(ics_bin2)

choice_list = np.arange(1, 3, 1)
for k in range(N):
    tip_rowplace = ics_bin1[k]
    tip_row = tip_rowplace[0]
    tip_row_str = "{}".format(tip_row)
    tip_row_str = tip_row_str.replace(', ', '').replace('[', '').replace(']', '')
    for j in range(Z):
        stalk_rowplace = ics_bin2[j]
        stalk_row = stalk_rowplace[0]
        stalk_row_str = "{}".format(stalk_row)
        stalk_row_str = stalk_row_str.replace(', ', '').replace('[', '').replace(']', '')
        tip_nrowplace = ics_bin1[j]
        stalk_nrowplace = ics_bin2[j]
        tip_stalk_nval = hamming_distance(tip_row, stalk_row)
        # here we will look at neighbor states and create an arrow
        randval = random.choice(choice_list)
        if tip_stalk_nval == 1 and randval == 2:
            t.edge(tip_row_str, stalk_row_str, label = '')
            t.attr(label=r'\n\nTip and Stalk Crosstalk: \nStable State Initial Conditions')
            t.attr(fontsize='50')

t.view()
"""""
with t.subgraph(name='tip cluster') as c:
    c.attr(label='Tip')
    c.attr('node', shape='square')
    c.edges(tip_arrow)
with t.subgraph(name='stalk cluster') as c:
    c.attr(label='Stalk')
    c.attr('node', shape='square')
    c.edges(stalk_arrow)
"""""




t.view()


