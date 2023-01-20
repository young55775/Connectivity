import sys
from collections import Counter
import matplotlib as mpl
mpl.use('TkAgg')
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns




class Neuron:  # define a neuron with post,pre,and ej
    def __init__(self, name: str, transmitter_weight=1, cat='unknown', transmitter='unknown', co_transmitter=None):
        self.name = name
        self.transmitter_weight = transmitter_weight
        self.post = []  # post_synapse
        self.pre = []  # pre_synapse
        self.ej = []  # gap_junction
        self.cat = cat
        self.transmitter = transmitter
        self.co_transmitter = co_transmitter
        self.input = []
        self.output = []

    def add_post(self, n1, n2, weight):
        self.post.append(Synapse(weight, n1, n2, self.transmitter, self.transmitter_weight, self.co_transmitter))

    def add_pre(self, n1, n2, weight):
        self.pre.append(Synapse(weight, n2, n1, self.transmitter, self.transmitter_weight, self.co_transmitter))

    def add_ej(self, n1, n2, weight):
        self.ej.append(GapJunc(weight=weight, pre=n1, post=n2))

    def input_signal(self, input: np.array):  #
        self.input.append(input)

    def generate_output(self):  # transduct the signal of neurons
        output = []
        for i in self.input:
            i = list(i)
            i = [0]*30 + i[:-30] # consider time frame
            output.append(i)
        output = np.mean(np.array(output), axis=0)
        m = output.copy()
        ind = np.abs(m)>0.5
        self.output = output*ind


class Synapse:  # synapse description including receptors and transmitters
    def __init__(self, n: float, pre, post, transmitter='unknown', transmitter_weight=1, co_transmitter='unknown'):
        self.transmitter = transmitter
        self.transmitter_weight = transmitter_weight
        self.weight = n
        self.pre = pre
        self.post = post
        self.ele = False
        self.co_transmitter = co_transmitter
        self.post_signals = []


class GapJunc(Synapse):
    def __init__(self, weight, pre, post, transmitter='ele', transmitter_weight=1, co='ele'):
        super().__init__(weight, pre, post, transmitter, transmitter_weight, co)
        self.ele = True


class NetWork:  # a neuron network created by neuron list
    def __init__(self, neuron_lst):
        self.neurons = neuron_lst
        self.weight_dict = {}
        self.G = nx.MultiDiGraph()
        self.G_ej = nx.DiGraph()
        self.G_chem = nx.MultiDiGraph()

    def add_neuron(self, neuron_lst):
        new = {}
        new.update(neuron_lst)
        new.update(self.neurons)
        self.neurons = new

    def construct(self):  # construct whole edges
        for neu in self.neurons.values():
            if neu.post:  # synapse list
                for syn in neu.post:
                    # if syn.pre in self.neurons.keys() and syn.post in self.neurons.keys():
                    self.G.add_edge(syn.pre, syn.post, weight=syn.weight, type='chem', transmitter=syn.transmitter,
                                    co_transmitter=syn.co_transmitter,transmitter_weight = syn.transmitter_weight)
                    self.G_chem.add_edge(syn.pre, syn.post, weight=syn.weight, transmitter=syn.transmitter,
                                         co_transmitter=syn.co_transmitter,transmitter_weight = syn.transmitter_weight)
            if neu.pre:
                for syn in neu.pre:
                    # if syn.pre in self.neurons.keys() and syn.post in self.neurons.keys():
                    self.G.add_edge(syn.pre, syn.post, weight=syn.weight, type='chem', transmitter=syn.transmitter,
                                    co_transmitter=syn.co_transmitter,transmitter_weight = syn.transmitter_weight)
                    self.G_chem.add_edge(syn.pre, syn.post, weight=syn.weight, transmitter=syn.transmitter,
                                         co_transmitter=syn.co_transmitter,transmitter_weight = syn.transmitter_weight)
            if neu.ej:
                for syn in neu.ej:
                # if syn.pre in self.neurons.keys() and syn.post in self.neurons.keys():
                    self.G.add_edge(syn.pre, syn.post, weight=syn.weight, type='gap', transmitter=syn.transmitter,
                                    co_transmitter=syn.co_transmitter,transmitter_weight = syn.transmitter_weight)
                    self.G.add_edge(syn.post, syn.pre, weight=syn.weight, type='gap', transmitter=syn.transmitter,
                                    co_transmitter=syn.co_transmitter,transmitter_weight = syn.transmitter_weight)
                    self.G_ej.add_edge(syn.pre, syn.post, weight=syn.weight, transmitter=syn.transmitter,
                                       co_transmitter=syn.co_transmitter,transmitter_weight = syn.transmitter_weight)
                    self.G_ej.add_edge(syn.post, syn.pre, weight=syn.weight, transmitter=syn.transmitter,
                                       co_transmitter=syn.co_transmitter,transmitter_weight = syn.transmitter_weight)

    def chem_plot(self):
        nx.draw_networkx(self.G_chem, with_labels=False, node_size=20, edge_color='coral',
                         label='chemical junction', width=0.3)

    def gap_plot(self):
        nx.draw_networkx(self.G_ej, with_labels=False, node_size=20, edge_color='deepskyblue', label='gap junction',
                         width=0.3)

    def total_plot(self, with_labels=False, node_size=20, width=0.3):
        edges, type = zip(*nx.get_edge_attributes(self.G, 'type').items())
        c_list = []
        for i in type:
            if i == 'chem':
                c_list.append('coral')
            else:
                c_list.append('deepskyblue')
        nx.draw_networkx(self.G, with_labels=with_labels, edge_color=c_list, node_size=node_size, label='network',
                         width=width)

    def reset_gap(self):  # reset the weight of gap junction
        pass

    def copy(self):
        return NetWork(neuron_lst=self.neurons.copy())

    def ablation(self, name):  # ablation a neuron and return a new network
        neuron_lst = self.neurons.copy()
        if name in neuron_lst.keys():
            del neuron_lst[name]
        else:
            print('Neuron {} not in the network'.format(name), file=sys.stderr)
        net = NetWork(neuron_lst)
        net.construct()
        return net

    def ablation(self, name_lst: list):
        new = self.copy()
        for i in name_lst:
            if i in new.neurons.keys():
                new = new.ablation(i)
        return new

    def activate(self,neuron,signal=np.sin(np.linspace(0, 5, 800))):
        self.neurons[neuron].input_signal(signal)
        self.neurons[neuron].generate_output()


    def sub_network(self, name_lst):  # subtract subnet from it
        neuron = {}
        for i in name_lst:
            neuron[i] = self.neurons[i]
        net = NetWork(neuron)
        net.construct()
        return net


def csv2net(path, class_path, trans_path, transmitter_dict={}):  # -> NetWork
    info = pd.read_csv(path)
    neu = info['Neuron 1'].to_list() + info['Neuron 2'].to_list()
    cat = pd.read_csv(class_path)
    cat = dict(zip(cat['name'].to_list(), cat['type'].to_list()))
    trans = pd.read_csv(trans_path)
    tran = dict(zip(trans['Neuron'].to_list(), trans['trans'].to_list()))
    for k in tran.keys():
        if not tran[k]:
            tran[k] = 'unknown'
    co_tran = dict(zip(trans['Neuron'].to_list(), trans['cotrans'].to_list()))
    neu = list(set(neu))
    neuron_list = {}
    for i in neu:
        if i in tran.keys() and i in cat.keys() and i in co_tran.keys():
            neuron_list[i] = Neuron(i, transmitter_dict[tran[i]], cat[i], tran[i], co_tran[i])
        else:
            neuron_list[i] = Neuron(i, transmitter_dict['unknown'], cat[i], 'unknown', 'unknown')
    for m in info.iterrows():
        n1 = m[1]['Neuron 1']
        n2 = m[1]['Neuron 2']
        t = m[1]['Type']
        v = m[1]['Nbr']
        if v != 0:
            if t == 'S' or t == 'Sb' or t == 'NMJ':
                neuron_list[n1].add_post(n1, n2, 1 / v)
            elif t == 'R' or t == 'Rp':
                neuron_list[n1].add_pre(n1, n2, 1 / v)
            elif t == 'EJ':
                neuron_list[n1].add_ej(n1, n2, 1 / v)
    connect = NetWork(neuron_list)
    connect.construct()
    return connect


def forward(neuron_list,network):
    '''
    :param neuron_list:list, graph:nx.MultiDiGraph():
    process the signal in every target neuron and give the output list
    :return: output list
    '''
    dic = network.neurons
    output = []
    for i in neuron_list:
        try:
            rec = dict(Counter(list(zip(*list(network.G.edges(i))))[1])) # find the numbers of edges and the post cells
            if len(dic[i].output) != 0:
                signal = dic[i].output
                for k,v in rec.items():
                    for n in range(v):
                        dt = network.G.edges[i,k,n]
                        w = (1/dt['weight'])*dt['transmitter_weight']
                        dic[k].input_signal(w*signal) # neuron i receive the signal from a synapse
                    dic[k].generate_output()
                    output.append(k)
        except:
            print(i)
    network.neurons = dic
    return output

def norm(ar):
    return (ar - np.min(ar))/(np.max(ar)-np.min(ar))

def wgn(sig,rate):
    noise = np.random.randn(len(sig))
    noise *= rate
    noise -= np.mean(noise)
    return sig+noise


if __name__ == '__main__':
    signal_asjl = 3*wgn(np.sin(5*np.linspace(0, 40, 800))+0.65*np.cos(np.linspace(12, 89, 800))*np.random.randn(800),0.1)
    signal_asjr = 3*wgn(np.sin(5*np.linspace(0, 40, 800)) + 0.65 * np.cos(np.linspace(12, 89, 800)) * np.random.randn(800), 0.1)
    signal_asel = 2*wgn(np.sin(3*np.linspace(0, 40, 800)),0.1)
    signal_aser = 2*wgn(np.cos(3*np.linspace(0, 40, 800)),0.1)
    signal_aval = 1.7*wgn(np.sin(np.linspace(0, 40, 800))+0.03*np.cos(0.005*np.linspace(12, 89, 800))*np.random.randn(800),0.1)
    signal_avar = 1.7*wgn(np.cos(np.linspace(0, 40, 800)) + 0.03 * np.sin(0.005*np.linspace(12, 89, 800)) * np.random.randn(800), 0.1)
    transmitter_dict = {'unknown': 0.5, 'GABA': -1.3, 'Glu': 0.7,
                        'Ach': 1, 'ACh': 1, 'Unknown (orphan)': 0.5,
                        'unknown MA (cat-1)': 0.5, 'DA': 0.7, 'ACh (minor)': 0.1,
                        'Octopamine': 1, '': 1}
    network = csv2net('connectome_all.csv', 'celltype.csv', 'transmitter.csv', transmitter_dict)
    network.activate('AVAL',signal_aval)
    network.activate('AVAR', signal_avar)
    network.activate('ASEL', signal_asel)
    network.activate('ASER', signal_aser)
    network.activate('ASJL', signal_asjl)
    network.activate('ASJR', signal_asjr)
    ax = []
    ns = ['AVAL','AVAR','ASEL','ASER','ASJL','ASJR']
    for i in range(3):
        ns = forward(ns,network)
        ns = list(set(ns))
        ax.extend(ns)
    m = {}
    ax = list(set(ax))
    for j in ax:
        m[j] = norm(network.neurons[j].output)
        # m['INPUT'] = norm(wgn(np.sin(np.linspace(0, 40, 800)),0.05))
    # a = sns.clustermap(pd.DataFrame(m).T, col_cluster=False, cmap='jet')
    plt.subplot(2,1,1)
    plt.plot(signal_asel,label = 'ASEL')
    plt.plot(signal_aser,label = 'ASER')
    plt.plot(signal_aval,label = 'AVAL')
    plt.plot(signal_avar,label = 'AVAR')
    plt.plot(signal_asjl,label = 'ASJL')
    plt.plot(signal_asjr,label = 'ASJR')
    plt.title('INPUT SIGNAL')
    plt.legend()
    plt.subplot(2,1,2)
    plt.plot(pd.DataFrame(m))
    plt.title('OUTPUT SIGNAL')