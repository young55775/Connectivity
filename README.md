# A *C. elegans* neuron connectome simulator

Analyze C. elegans neuron activity  
based on networkx 3.0
___
# Connectome.py
## 1 Class Neuron(*name:str | cat:str | trans:str | co_trans:str*)
### 1.1 Attribute
- **name:str**| Name of the neuron  
- **cat:str** | Category of 'sensory','interneuron','motorneuron'  
- **transmitter,co_transmitter:str** | neuron transmitter  
- **pre,post:list** | list of pre/post neuron synapses (see "Class Synapse")  
- **ej:list** | list of gap junctions (see "Class GapJunc")  
____
### 1.2 Method
pass
___
## Class Synapse(*weight:float | pre:str | post:str | transmitter:str | co:str*)
### 1.1 Attribute
- **pre/post:str** | name of pre/post synaptic neurons
- **transmitter/co_transmitter":str** | neuron transmitter of the post synaptic cell
- **weight:float** | 1/(number of synapses between pre-post neurons)
- **ele:bool** | is electric synapse or not
____
### 1.2 Method
pass
___
## Class GapJunc(*factor:float | weight:float | pre:str | post:str | transmitter:str | co:str*)

define this class for the switch of the weight of gap junction 
### 1.1 Attribute
```
class GapJunc(Synapse):
    super().__init__
```
- **gap_weight:float** | factor(default=1) * self.weight

### 1.2 Method
pass
___
## Class NetWork(*neuron_list:dict*)
### 1.1 Attribute
- **neurons:dict** | dictionary {name:class Neuron}  
- **G:MultiDiGraph** | networkx MultiDiGraph of all synapse  
- **G_ej:MultiDiGraph** | networkx MultiDiGraph of gap junction
- **G_chem:MultiDiGraph** | networkx MultiDiGraph of chemical synapse
___
### 1.2 Method
#### 1.2.1 add_neuron(*self*,neuron_lst)
append a novel neuron list to the attribute-neurons
_____
#### 1.2.2 construct(*self*)
construct the **self.G/self.G_ej/self.G_chem**
___
#### 1.2.3 total_plot(*self*,with_labels:bool = False,node_size:float = 20,width:float = 0.3)
plot self.G with gap junction colored 'deepskyblue' and chemical synapses colored 'coral' 
___
#### 1.2.4 gap_plot/chem_plot(*self*)
plot **self.G_ej/self.G_chem**
___
#### 1.2.5 ablation(*self*,name:str)
return a copy of the network with a neuron deleted
___
#### 1.2.6 sub_network(*self*,name_lst:list)
return a copy of the network constructed with the selected neurons in the name_lst

