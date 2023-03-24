import preprocess as pr
import connectome as connect
import networkx as nx
import pandas as pd

ch1_p = r"ch1.csv"
ch2_p = r"ch2.csv"
ch1 = pd.read_csv(ch1_p, index_col=0)
ch2 = pd.read_csv(ch2_p, index_col=0)
name = ch1.columns
# or can manually remove error values
# ch1 = ch1[120:]
# ch2 = ch2[120:]
ch1_d = ch1.values
ch2_d = ch2.values
fret_hm = ch2_d / ch1_d
fret_hm = fret_hm.T
fret_hm_dn = pr.de_noise(fret_hm)
fret_hm_dt, curve = pr.de_trend(fret_hm_dn)
fret_hm_nm = pr.norm(fret_hm_dt.T)
fret = pd.DataFrame(fret_hm_nm, columns=name)
#elect 40 most activated neurons via std
#ele = pr.elect(fret,40)
co, pos, neg, null = pr.cor(fret, 0.3)
type_dict = {'sensory': 1.3, 'interneuron': 1.1, 'motorneuron': 0.8, 'muscle': 0.1}
transmitter_dict = {'unknown': 0.5, 'GABA': -1.4, 'Glu': 1,
                    'Ach': 1.1, 'ACh': 1.1, 'Unknown (orphan)': 0.5,
                    'unknown MA (cat-1)': 0.5, 'DA': 0.7, 'ACh (minor)': 0.3,
                    'Octopamine': 1.3, '': 1}
network = connect.csv2net('connectome_all.csv', 'celltype.csv', 'transmitter.csv', transmitter_dict, type_dict)
Graph = network.G_chem
path_dict = network.count_path(4,pos)

plt.bar(range(len(path_dict)), list(path_dict.values()), align='center')
plt.xticks(range(len(path_dict)), list(path_dict.keys()))
plt.show()

a = list(path_dict.items())
a = sorted(a, key = lambda x:x[1], reverse = True)
dict(a)