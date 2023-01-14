import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
# make a connectome network with annotate file and connectomic dataset
a = pd.read_csv('connectomic.csv')
con = list(zip(a['Neuron 1'].to_list(),a['Neuron 2'].to_list()))
spot = a['Neuron 1'].to_list()
neuron = list(pd.read_csv('ch1.csv').columns)
con = [n for n in con if n[1] in neuron and n[0] in neuron]
con = list(set(con))
table = {}
for i in con:
    if i[0] not in table.keys():
        table[i[0]] = []
    table[i[0]].append(i[1])
art = []
for k, v in table.items():
    art.append(k + ' ' + str(len(v)) + '\n' + '\n'.join(v))
conn = '\n'.join(art)
with open('connectome_network', 'w+') as f:
    f.write(conn)
conn = nx.read_multiline_adjlist('connectome_network')
# nx.draw(conn)

def clean_lst(lst,spot):
    return [n for n in lst if n[0] in spot and n[1] in spot]
def dist_cal(net,lst):
    dist = []
    for i in lst:
        dist.append(nx.shortest_paths.shortest_path_length(net,i[0],i[1]))
    return dist
