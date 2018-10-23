import itertools as it

class DATA(object):
        #A 'DATA' object stores all initial info needed to create ReadGraph G
        #data contains read info
        #states is the state of each node, in diploid this is always 1, counts the number of mutant alleles
        #k = ploidy
        #error = error rate assumed
        #read_list is list of reads
        #positions is the actual position on the genome (1-3bb) in base pair count
        def __init__(self,data,states,k,error,read_list,positions,names,chroms):
                self.data = data
                self.states = states
                self.error = error
                self.read_list = read_list
                self.k = k
                self.E = data.keys()

                self.nodekeys = sorted(list(self.node_keys()))
                self.ADJ = self.adj()

                self.positions = positions
                self.names = names
                self.chroms = chroms
                d={}
                for num in range(len(self.nodekeys)):
                        index = self.nodekeys[num]
                        d[index] = node(self.k,index,self.states[index],self.ADJ[index],self.positions[index],self.names[index],self.chroms[index],num)
                self.nodes = d

                N = self.nodes
                d={}
                for e in self.E:
                        d[e] = edge([N[e[0]],N[e[1]]],self.k)
                self.edges = d

        def node_keys(self):
                #computes all nodes actually seen by reads
                s = set()
                for edge in self.E:
                        s.add(edge[0])
                        s.add(edge[1])
                return s

        def adj(self):
                #returns dictionary mapping NODE to those nodes adjacent to NODE
                d={}
                for i in self.nodekeys:
                        d[i] = []

                for e in self.E:
                        d[e[0]].append(e[1])
                        d[e[1]].append(e[0])
                return d


def make_data_from_RNA_data(RD):
        read_list = RD.long_read_list
        data = edges_from_readlist(read_list)
        return DATA(data,RD.states,RD.k,RD.error,read_list,RD.positions,RD.names,RD.chroms)



def edges_from_readlist(read_list):
        #returns dictionary with keys those edges that appear in G and values all 1
        d={}
        for zz in read_list:
                r = read_list[zz]
                seen_nodes = r.keys
                pairs = list(it.combinations(seen_nodes,2))
                for pair in pairs:
                        i,j = pair
                        pair = order(i,j)
                        d[pair]=1
        return d

def order(i,j):
        #correctly orders i<j
        if i<j:
                return (i,j)
        elif i>j:
                return (j,i)
        else:
                print 'Error'



class easy_graph(object):
        #ReadGraph object, takes in data object
        def __init__(self,data):
                print 'Making ReadGraph'
                self.data = data
                self.error = data.error
                self.nodes = data.nodes
                self.edges = data.edges
                self.nodekeys = self.data.nodekeys
                self.size = len(self.data.nodekeys)
                print str(self.size) + ' SNPs in non-trivial connected components'
                self.k = self.data.k
                self.states = self.data.states

                self.components,self.comps_dict = build_all_comps(self)
                self.comp_mins = mins_of_comps(self.components)
                self.names = data.names
                self.chroms = data.chroms


                self.read_list = data.read_list
                self.POs = [[(0, 0)], [(0, 1), (1, 0)], [(1, 1)]]
                self.OPs = {(0, 1): 0, (1, 0): 1}
                self.POs
                self.read_dict,self.comp_reads = self.make_read_dict()

        def __iter__(self):
            return self.nodes.values().__iter__()

        def make_position_dict(self):
            position_to_SNP_index_dict = {}
            for s in self.nodekeys:
                S = self.nodes[s]
                p = S.position
                position_to_SNP_index_dict[p] = S.index
            return position_to_SNP_index_dict


        def make_read_dict(self):

            comp_reads = {m:[] for m in self.comp_mins}
            read_dicts = {s:[] for s in self.nodekeys}

            for r in self.read_list.values():
                m = self.comps_dict[r.keys[0]]
                comp_reads[m].append(r)
                for key in r.keys[1:]:
                    read_dicts[key].append(r)


            return read_dicts,comp_reads





class node(object):
        def __init__(self,k,index,shape,neighbors,pos,name,chrom,num):
                self.k = k
                self.shape = shape
                self.index = index
                self.neighbors = neighbors
                self.position = pos
                self.name = name
                self.chrom = chrom
                self.num = num



        def __eq__(self,other):
                return self.index == other.index

        def __lt__(self,other):
                if self.index < other.index:
                        return True
                return False




class edge(object):
        def __init__(self,node_pair,k):
                self.nodes = node_pair
                self.n0 = node_pair[0]
                self.n1 = node_pair[1]
                self.indices = (self.n0.index,self.n1.index)
                self.k = k
                self.shapes = (self.n0.shape,self.n1.shape)

        def __eq__(self,other):
                return self.indices == other.indices





def build_all_comps(G):
        #builds connected component dictionaries
        #comps maps each node to a list of nodes it is connected to
        #comps_dict maps each node to the smallest node it is connected to
            comps = {}
            comps_dict = {}
            unbuilt = set(G.nodekeys)
            while len(unbuilt) > 0:
                    start = unbuilt.pop()
                    s = build_comp(G,start)
                    comps[start] = s
                    comps_dict[start] = min(s)
                    unbuilt.difference_update(s) # we build_comp once for each set of components
                    for key in s:  #and apply it to each nodekey
                        comps[key] = s
                        comps_dict[key]=comps_dict[start]
            return comps,comps_dict

def build_comp(G,start):
        #builds connected component containing the node start
        #returns ordered list of nodes in component
    s = set([start])
    old_batch = set([start])
    while len(old_batch) > 0:
            new_batch = set([])
            for x in old_batch:
                    new_batch.update(set(G.nodes[x].neighbors))
            new_batch = new_batch.difference(s)
            s.update(new_batch)
            old_batch = new_batch

    return sorted(list(s))


def mins_of_comps(comps):
        mins = {}
        for comp in comps.values():
                mins[min(comp)] = 1
        return sorted(mins.keys())


def edges_from_readlist(read_list):
        #returns dictionary with keys those edges that appear in G and values all 1
        d={}
        for zz in read_list:
                r = read_list[zz]
                seen_nodes = r.keys
                pairs = list(it.combinations(seen_nodes,2))
                for pair in pairs:
                        i,j = pair
                        pair = order(i,j)
                        d[pair]=1
        return d
