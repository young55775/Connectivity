import matplotlib as mpl
import networkx as nx
import pandas as pd

mpl.use('TkAgg')
import matplotlib.pyplot as plt


class Neuron:  # define a neuron with post,pre,and ej
    def __init__(self, name: str):
        self.name = name
        self.post = []  # post_synapse
        self.pre = []  # pre_synapse
        self.ej = []  # gap_junction

    def add_post(self, n1, n2, weight):
        self.post.append(Synapse(weight, n1, n2))

    def add_pre(self, n1, n2, weight):
        self.pre.append(Synapse(weight, n2, n1))

    def add_ej(self, n1, n2, weight):
        self.ej.append(GapJunc(factor=1, weight=weight, pre=n1, post=n2))


class SensoryNeuron(Neuron):  # Different neuron class for further information gain
    def __init__(self):
        super().__init__()


class MotorNeuron(Neuron):
    def __init__(self):
        super().__init__()


class InterNeuron(Neuron):
    def __init__(self):
        super().__init__()


class Synapse:  # synapse description including receptors and transmitters
    def __init__(self, n: float, pre, post):
        self.receptor = []
        self.transmitter = []
        self.weight = n
        self.pre = pre
        self.post = post
        self.ele = False

    def signal(self):
        pass


class GapJunc(Synapse):
    def __init__(self, factor, weight, pre, post):
        super().__init__(weight, pre, post)
        self.gap_weight = factor * self.weight
        self.ele = True


class NetWork:  # a neuron network created by neuron list
    def __init__(self):
        self.neurons = []
        self.G = []  # total graph
        self.G_chem = []  # without ej
        self.G_ej = []  # ej only

    def add_neuron(self, neuron_lst):
        self.neurons += neuron_lst

    def construct(self):
        self.G = nx.DiGraph()
        self.G_ej = nx.DiGraph()
        self.G_chem = nx.DiGraph()
        for neu in self.neurons:
            if neu.post:  # synapse list
                for syn in neu.post:
                    self.G.add_edge(syn.pre, syn.post, weight=syn.weight, type='chem')
                    self.G_chem.add_edge(syn.pre, syn.post, weight=syn.weight)
            if neu.pre:
                for syn in neu.pre:
                    self.G.add_edge(syn.pre, syn.post, weight=syn.weight, type='chem')
                    self.G_chem.add_edge(syn.pre, syn.post, weight=syn.weight)
            if neu.ej:
                self.G.add_edge(syn.pre, syn.post, weight=syn.weight, type='gap')
                self.G.add_edge(syn.post, syn.pre, weight=syn.weight, type='gap')
                self.G_ej.add_edge(syn.pre, syn.post, weight=syn.weight)
                self.G_ej.add_edge(syn.post, syn.pre, weight=syn.weight)

    def chem_plot(self):
        nx.draw_networkx(self.G_chem, with_labels=False, node_size=20, edge_color='coral',
                         label='chemical junction', width=0.3)

    def gap_plot(self):
        nx.draw_networkx(self.G_ej, with_labels=False, node_size=20, edge_color='deepskyblue', label='gap junction',
                         width=0.3)

    def total_plot(self):
        edges, type = zip(*nx.get_edge_attributes(self.G, 'type').items())
        c_list = []
        for i in type:
            if i == 'chem':
                c_list.append('coral')
            else:
                c_list.append('deepskyblue')
        nx.draw_networkx(self.G, with_labels=False, edge_color=c_list, node_size=20, label='network',
                         width=0.3)

    def reset_gap(self): # reset the weight of gap junction
        pass

    def copy(self):
        pass

    def ablation(self, name):
        pass

    def ablation_chemical(self, name):
        pass

    def ablation_gap(self, name):
        pass


def csv2net(path):  # -> NetWork
    info = pd.read_csv(path)
    neu = info['Neuron 1'].to_list() + info['Neuron 2'].to_list()
    neu = list(set(neu))
    neuron_list = {}
    for i in neu:
        neuron_list[i] = Neuron(i)
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
            else:
                neuron_list[n1].add_ej(n1, n2, 1 / v)
    connect = NetWork()
    connect.add_neuron(neuron_list.values())
    connect.construct()
    return connect


if __name__ == '__main__':
    network = csv2net(r'D:\annotation\Connectivity\connectome_all.csv')
