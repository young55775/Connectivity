import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from collections import Counter
# make a connectome network with annotate file and connectomic dataset
# digraph with weight
a = pd.read_csv('connectome_all.csv')
adj = a.values
lst = []
connect = []
for i in adj:
    if i[2] == 'S' or i[2] == 'Sp' or i[2] == 'NMJ':
        if [i[0],i[1]] not in connect and i[3] != 0:
            lst.append([i[0],i[1],1/i[3]])
            connect.append([i[0],i[1]])
    elif i[2] == 'R' or i[2] == 'Rp':
        if [i[1], i[0]] not in connect and i[3] != 0:
            lst.append([i[1], i[0], 1 / i[3]])
            connect.append([i[0], i[1]])
    else:  #gap junction
        if i[3] != 0:
            lst.append([i[1], i[0], 1 / 0.7*i[3]])
            lst.append([i[0], i[1], 1 / 0.7*i[3]])



def construct(lst):
    G = nx.DiGraph()
    for i in lst:
        G.add_edge(i[0],i[1],weight=i[2])
    return G



#calculate distance


#draw positive/negative network
def pro_lst(lst,ref,g,remove=[]):
    p = list(zip(*lst))
    p = list(p[1]) + list(p[0])
    p = list(set(p))
    res = []
    ran = []
    lst = [n for n in lst if n[0] in g.nodes() and n[1] in g.nodes()]
    for i in lst:
        ran += nx.shortest_paths.shortest_path(g,i[0],i[1])
    ran = list(set(ran))
    ran = [n for n in ran if n not in remove]
    for j in ref:
        if j[0] in p or j[1] in p:
            if j[0] in ran and j[1] in ran:
                res.append(j)
    return res

def calculate_path(g,lst):
    p = []
    lst = [n for n in lst if n[0] in g.nodes() and n[1] in g.nodes()]
    for i in lst:
        try:
            p.append(nx.shortest_paths.dijkstra_path_length(g,i[0],i[1],'weight'))
        except:
            print()
    return p

def get_product_path(G,node1,node2):
    nodelst = nx.dijkstra_path(G,node1,node2)
    path = 1
    for i in range(len(nodelst)-1):
        path /= nx.dijkstra_path_length(G,nodelst[0],nodelst[1])
    return 1/path

def get_count(G,lst):
    nodelst = []
    for i in lst:
        try:
            nodelst.extend(nx.dijkstra_path(G,i[0],i[1]))
        except:
            print(i)
    return dict(Counter(nodelst))

# lst = [n for n in lst if n[0] in name or n[1] in name]
pos = [n for n in pos if n[0] in a['Neuron 1'].to_list() and n[1] in a['Neuron 1'].to_list()]
neg = [n for n in neg if n[0] in a['Neuron 1'].to_list() and n[1] in a['Neuron 1'].to_list()]
G=construct(lst)
nx.draw_networkx(G,with_labels=False,node_size=20,width=0.5)
plt.show(block=True)
# lstp = pro_lst(pos,lst,G)
# gp = construct(lstp)
# lstn = pro_lst(neg,lst,G)
# gn = construct(lstn)
# nx.draw_networkx(gp,with_labels=False,node_size=20,edge_color='coral',label='positive')
# nx.draw_networkx(gn,with_labels=False,node_size=20,edge_color='deepskyblue',label='negative')
# plt.show(block=True)
# pos_dist = calculate_path(g,pos)
# neg_dist = calculate_path(g,neg)

#find neuron importance
pos_count = get_count(G,pos)
neg_count = get_count(G,neg)
n = []
x = []
y = []
dic = [n[0] for n in neg]
neg_spe = [n[0] for n in neg if n[0] not in pos_count.keys()]
pos_spe = [n[0] for n in pos if n[0] not in neg_count.keys()]
for k, v in pos_count.items():
    if k in neg_count.keys():
        x.append(v)
        y.append(neg_count[k])
        n.append(k)

n = n + neg_spe
x = x + [0] * len(neg_spe)
y = y + [neg_count[k] for k in neg_spe]
n = n + pos_spe
x = x + [pos_count[k] for k in pos_spe]
y = y + [0] * len(pos_spe)
ax1 = plt.subplot(2,1,1)
plt.scatter(x, y)
ind_d = []
for i in range(len(x)):
    j = x[i]
    k = y[i]
    while (j, k) in ind_d:
        k += 0.22
        j += 0.12
    plt.text(j, k, n[i], {'size': 7})
    ind_d.append((j, k))
plt.xlabel('positive')
plt.ylabel('negative')

# correlationship between cor and path
corr = []
dist = []
name = list(ele.columns)
for i in name:
    for j in name:
        if i != j:
            try:
                dist.append(nx.shortest_paths.dijkstra_path_length(G, i, j))
                corr.append(stats.spearmanr(ele[i], ele[j])[0])
            except:
                print(i, j)

comp = list(zip(corr, dist))
c_pos = [n for n in comp if n[0] > 0.4]
c_neg = [n for n in comp if n[0] < -0.4]
cor_pos,d_pos = zip(*c_pos)
cor_neg,d_neg = zip(*c_neg)
ax2 = plt.subplot(2,1,2)
plt.scatter(corr,dist)
plt.scatter(cor_neg,d_neg,c='deepskyblue',label = 'neg cor = {}'.format(stats.spearmanr(cor_neg,d_neg)[0]))
plt.scatter(cor_pos,d_pos,c="coral",label='pos cor = {}'.format(stats.spearmanr(cor_pos,d_pos)[0]))
plt.xlabel('correlation')
plt.ylabel('Synaptic distance')
plt.legend(loc='upper right')
