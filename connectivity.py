import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
# make a connectome network with annotate file and connectomic dataset
# digraph with weight
a = pd.read_csv('connectomic.csv')
adj = a.values
lst = []
connect = []
for i in adj:
    if i[2] == 'S' or i[2] == 'Sp':
        if [i[0],i[1]] not in connect and i[3] != 0:
            lst.append([i[0],i[1],1/i[3]])
            connect.append([i[0],i[1]])
    elif i[2] == 'R' or i[2] == 'Rp':
        if [i[1], i[0]] not in connect and i[3] != 0:
            lst.append([i[1], i[0], 1 / i[3]])
            connect.append([i[0], i[1]])
    else:
        if i[3] != 0:
            lst.append([i[1], i[0], 1 / i[3]])
            lst.append([i[0], i[1], 1 / i[3]])

lst = [n for n in lst if n[0] in name or n[1] in name]

def construct(lst):
    G = nx.DiGraph()
    for i in lst:
        G.add_edge(i[0],i[1],weight=i[2])
    return G

nx.draw_networkx(G,with_labels=False,node_size=20)
plt.show(block=True)

#calculate distance
pos = [n for n in pos if n[0] in a['Neuron 1'].to_list() and n[1] in a['Neuron 1'].to_list()]
neg = [n for n in neg if n[0] in a['Neuron 1'].to_list() and n[1] in a['Neuron 1'].to_list()]

#draw positive/negative network
def pro_lst(lst,ref,ran):
    p = list(zip(*lst))
    p = list(p[1]) + list(p[0])
    p = list(set(p))
    res = []
    for j in ref:
        if j[0] in p or j[1] in p:
            if j[0] in ran and j[1] in ran:
                res.append(j)
    return res

def calculate_path(g,lst):
    p = []
    for i in lst:
        p.append(nx.shortest_paths.shortest_path_length(g,i[0],i[1]))
    return p


# lstp = pro_lst(pos,lst,name)
# gp = construct(lstp)
# lstn = pro_lst(neg,lst,name)
# gn = construct(lstn)
# nx.draw_networkx(gp,with_labels=False,node_size=20,edge_color='coral',label='positive')
# nx.draw_networkx(gn,with_labels=False,node_size=20,edge_color='deepskyblue',label='negative')
# plt.show(block=True)
# pos_dist = calculate_path(g,pos)
# neg_dist = calculate_path(g,neg)
