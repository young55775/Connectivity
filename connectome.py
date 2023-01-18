import sys

import matplotlib as mpl
import networkx as nx
import pandas as pd

mpl.use('TkAgg')
import matplotlib.pyplot as plt


class Neuron:  # define a neuron with post,pre,and ej
    def __init__(self, name: str, cat = 'unknown',trans = 'unknown',co_trans = None):
        self.name = name
        self.post = []  # post_synapse
        self.pre = []  # pre_synapse
        self.ej = []  # gap_junction
        self.cat = cat
        self.transmitter = trans
        self.co_tranmitter = co_trans

    def add_post(self, n1, n2, weight,tran,co):
        self.post.append(Synapse(weight, n1, n2,tran,co))

    def add_pre(self, n1, n2, weight,tran,co):
        self.pre.append(Synapse(weight, n2, n1,tran,co))

    def add_ej(self, n1, n2, weight):
        self.ej.append(GapJunc(factor=1, weight=weight, pre=n1, post=n2))


class Synapse:  # synapse description including receptors and transmitters
    def __init__(self, n: float, pre, post,transmitter='unknown',co='unknown'):
        self.transmitter = transmitter
        self.weight = n
        self.pre = pre
        self.post = post
        self.ele = False
        self.co_transmitter = co


class GapJunc(Synapse):
    def __init__(self, factor, weight, pre, post,transmitter='ele',co='ele'):
        super().__init__(weight, pre, post,transmitter,co)
        self.gap_weight = factor * self.weight
        self.ele = True


class NetWork:  # a neuron network created by neuron list
    def __init__(self,neuron_lst):
        self.neurons = neuron_lst
        self.G = nx.MultiDiGraph()
        self.G_ej = nx.MultiDiGraph()
        self.G_chem = nx.MultiDiGraph()

    def add_neuron(self, neuron_lst):
        new = {}
        new.update(neuron_lst)
        new.update(self.neurons)
        self.neurons = new

    def construct(self): # construct whole edges
        for neu in self.neurons.values():
            if neu.post:  # synapse list
                for syn in neu.post:
                    # if syn.pre in self.neurons.keys() and syn.post in self.neurons.keys():
                        self.G.add_edge(syn.pre, syn.post, weight=syn.weight, type='chem',transmitter = syn.transmitter,co_transmitter = syn.co_transmitter)
                        self.G_chem.add_edge(syn.pre, syn.post, weight=syn.weight,transmitter = syn.transmitter,co_transmitter = syn.co_transmitter)
            if neu.pre:
                for syn in neu.pre:
                    # if syn.pre in self.neurons.keys() and syn.post in self.neurons.keys():
                        self.G.add_edge(syn.pre, syn.post, weight=syn.weight, type='chem',transmitter = syn.transmitter,co_transmitter = syn.co_transmitter)
                        self.G_chem.add_edge(syn.pre, syn.post, weight=syn.weight,transmitter = syn.transmitter,co_transmitter = syn.co_transmitter)
            if neu.ej:
                # if syn.pre in self.neurons.keys() and syn.post in self.neurons.keys():
                    self.G.add_edge(syn.pre, syn.post, weight=syn.weight, type='gap',transmitter = syn.transmitter,co_transmitter = syn.co_transmitter)
                    self.G.add_edge(syn.post, syn.pre, weight=syn.weight, type='gap',transmitter = syn.transmitter,co_transmitter = syn.co_transmitter)
                    self.G_ej.add_edge(syn.pre, syn.post, weight=syn.weight,transmitter = syn.transmitter,co_transmitter = syn.co_transmitter)
                    self.G_ej.add_edge(syn.post, syn.pre, weight=syn.weight,transmitter = syn.transmitter,co_transmitter = syn.co_transmitter)

    def chem_plot(self):
        nx.draw_networkx(self.G_chem, with_labels=False, node_size=20, edge_color='coral',
                         label='chemical junction', width=0.3)

    def gap_plot(self):
        nx.draw_networkx(self.G_ej, with_labels=False, node_size=20, edge_color='deepskyblue', label='gap junction',
                         width=0.3)

    def total_plot(self,with_labels=False,node_size=20,width=0.3):
        edges, type = zip(*nx.get_edge_attributes(self.G, 'type').items())
        c_list = []
        for i in type:
            if i == 'chem':
                c_list.append('coral')
            else:
                c_list.append('deepskyblue')
        nx.draw_networkx(self.G, with_labels=with_labels, edge_color=c_list, node_size=node_size, label='network',
                         width=width)

    def reset_gap(self): # reset the weight of gap junction
        pass

    def copy(self):
        pass

    def ablation(self, name): # ablation a neuron and return a new network
        neuron_lst = self.neurons.copy()
        if name in neuron_lst.keys():
            del neuron_lst[name]
        else:
            print('Neuron {} not in the network'.format(name), file=sys.stderr)
        net = NetWork(neuron_lst)
        net.construct()
        return net

    def ablation_chemical(self, name):
        pass

    def ablation_gap(self, name):
        pass

    def sub_network(self,name_lst): #subtract subnet from it
        neuron = {}
        for i in name_lst:
            neuron[i] = self.neurons[i]
        net = NetWork(neuron)
        net.construct()
        return net


def csv2net(path,class_path,trans_path):  # -> NetWork
    info = pd.read_csv(path)
    neu = info['Neuron 1'].to_list() + info['Neuron 2'].to_list()
    cat = pd.read_csv(class_path)
    cat = dict(zip(cat['name'].to_list(), cat['type'].to_list()))
    trans = pd.read_csv(trans_path)
    tran = dict(zip(trans['Neuron'].to_list(), trans['trans'].to_list()))
    co_tran = dict(zip(trans['Neuron'].to_list(), trans['cotrans'].to_list()))
    neu = list(set(neu))
    neuron_list = {}
    for i in neu:
        if i in tran.keys() and i in cat.keys() and i in co_tran.keys():
            neuron_list[i] = Neuron(i,cat[i],tran[i],co_tran[i])
        else:
            neuron_list[i] = Neuron(i,cat[i])
    for m in info.iterrows():
        n1 = m[1]['Neuron 1']
        n2 = m[1]['Neuron 2']
        t = m[1]['Type']
        v = m[1]['Nbr']
        if v != 0:
            if t == 'S' or t == 'Sb' or t == 'NMJ':
                neuron_list[n1].add_post(n1, n2, 1 / v,neuron_list[n1].transmitter,neuron_list[n1].co_tranmitter)
            elif t == 'R' or t == 'Rp':
                neuron_list[n1].add_pre(n1, n2, 1 / v,neuron_list[n2].transmitter,neuron_list[n2].co_tranmitter)
            else:
                neuron_list[n1].add_ej(n1, n2, 1 / v)
    connect = NetWork(neuron_list)
    connect.construct()
    return connect


if __name__ == '__main__':
    network = csv2net('connectome_all.csv','celltype.csv','transmitter.csv')
