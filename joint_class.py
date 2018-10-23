import basic_class

class joint_graph(object):
        def __init__(self,RD,G):
                self.RD = RD
                self.G = G
                self.read_list  = RD.all_reads
                self.k = 2
                self.error = G.error
                self.short_components = RD.components  #components induced by genes  / DASE
                self.short_comp_mins = sorted(RD.comp_mins)
                temp_D = {}
                for x in self.short_comp_mins:
                    for y in self.short_components[x]:
                        temp_D[y] = x
                
                self.short_comps_dict = temp_D
                
                self.long_comp_mins = sorted(G.comp_mins) #components induced by 2+-reads (HapTree)
                self.long_components = G.components
                self.long_comps_dict = G.comps_dict
                

                self.nodes = self.organize_nodes()
                self.nodekeys = sorted(self.nodes.keys())

                self.weird_components,self.weird_comps_dict = basic_class.build_all_comps(self)
                self.weird_comp_mins = basic_class.mins_of_comps(self.weird_components)

                
                self.components, self.comp_mins, self.comps_dict = self.translate_components()
                self.comp_mins = basic_class.mins_of_comps(self.components)
                self.read_dict,self.comp_reads = self.make_read_dict()

                
        def organize_nodes(self):
                nodes = {}
                Ls = len(self.short_comp_mins)
                self.Ls = Ls
                Ll = len(self.long_comp_mins)
                self.Ll = Ll

                forward_short = {}
                forward_long = {}
                back_short = {}
                back_long = {}
                for i in range(Ls):
                    Ms = self.short_comp_mins[i]
                    forward_short[i] = self.short_components[Ms]
                    back_short[Ms] = i
                for i in range(Ls,Ls+Ll):
                    Ml = self.long_comp_mins[i-Ls]
                    forward_long[i] = self.long_components[Ml]
                    back_long[Ml] = i

                for i in range(Ls):
                    temp_neighbors = set()
                    for lil_node in self.short_components[self.short_comp_mins[i]]:
                        if self.long_components.has_key(lil_node):
                            temp_neighbors.add(back_long[self.long_comps_dict[lil_node]])

                    nodes[i] = mini_node(i,sorted(list(temp_neighbors)))

                for i in range(Ls,Ll+Ls):
                    temp_neighbors = set()
                    for lil_node in self.long_components[self.long_comp_mins[i-Ls]]:
                        if self.short_components.has_key(lil_node):
                            temp_neighbors.add(back_short[self.short_comps_dict[lil_node]])

                    nodes[i] = mini_node(i,sorted(list(temp_neighbors)))

                self.forward_short = forward_short
                self.forward_long = forward_long
                self.back_short = back_short
                self.back_long = back_long

                return nodes

        def translate_components(self):
                comp_mins = []
                components = {}
                comps_dict = {}
                for c in self.weird_components:
                    temp_comp = set()
                    for C in self.weird_components[c]:
                            if C < self.Ls:
                                    for x in self.forward_short[C]:
                                        temp_comp.add(x)
                            else:
                                    for x in self.forward_long[C]:
                                        temp_comp.add(x)
                                        
                    comp = sorted(list(temp_comp))
                    m = min(comp)
                    comp_mins.append(m)
                    for s in comp:
                        components[s] = comp
                        comps_dict[s] = m
                        
                return components, sorted(comp_mins), comps_dict
                
                    

        def make_read_dict(self):
                #This does not work correctly if G is not the 2-reads from RNA-seq data.
                read_dict = {}
                for s in self.RD.components:
                    read_dict[s] = []
                    for r in self.RD.read_dict[s]: #1-reads
                        read_dict[s].append(r)
                comp_reads = {m:[] for m in self.comp_mins}
                for r in self.G.read_list.values():
                    m = self.comps_dict[r.keys[0]]
                    comp_reads[m].append(r)
                    for k in r.keys:
                        if read_dict.has_key(k):
                            read_dict[k].append(r)
                        else:
                            read_dict[k] = [r]

                            
                return read_dict,comp_reads
                

class mini_node(object):
        def __init__(self,index,neighbors):
                self.index = index
                self.neighbors = neighbors
                             
