import preprocess as pr
import connectome as connect
import networkx as nx
import pandas as pd
import matplotlib as mpl

mpl.use('TkAgg')
import matplotlib.pyplot as plt

ch1_p = r"ch1.csv"
ch2_p = r"ch2.csv"
y1 = main(ch1_p, ch2_p)  # generate normalized mat
sns.clustermap(y1, col_cluster=False)  # draw a signal fig
co_eff = np.corrcoef(y1.T)  # co_eff matrix
sns.clustermap(co_eff)  # heatmap of co_eff matrix
pr.histplot(co_eff, label='young')
plt.legend()
plt.show()
# elect 40 most activated neurons via std
# ele = pr.elect(fret,40)
co, pos, neg, null = pr.cor(fret, 0.3)
type_dict = {'sensory': 1.3, 'interneuron': 1.1, 'motorneuron': 0.8, 'muscle': 0.1}
transmitter_dict = {'unknown': 0.5, 'GABA': -1.4, 'Glu': 1,
                    'Ach': 1.1, 'ACh': 1.1, 'Unknown (orphan)': 0.5,
                    'unknown MA (cat-1)': 0.5, 'DA': 0.7, 'ACh (minor)': 0.3,
                    'Octopamine': 1.3, '': 1}
network = connect.csv2net('connectome_all.csv', 'celltype.csv', 'transmitter.csv', transmitter_dict, type_dict)
Graph = network.G_chem
path_dict = network.count_path(4, pos)

plt.bar(range(len(path_dict)), list(path_dict.values()), align='center')
plt.xticks(range(len(path_dict)), list(path_dict.keys()))
plt.show()

a = list(path_dict.items())
a = sorted(a, key=lambda x: x[1], reverse=True)
dict(a)
