import random

class baby_graph(object):
    def __init__(self,tree_sol):
        self.tree_sol = tree_sol
        self.comp_mins = []
        self.components = {}
        self.comp_dict = {}
        self.phase = {}
        self.MEC = 0

    def make_info(self):
        f = open(self.tree_sol,'r')
        a = f.readlines()
        f.close()
        comp_mins = []
        components = {}
        comp_dict = {}
        phase = {}
        mec = 0
        for k in a:
            if k[0] =='B':
                line = k.split(' ')
                m = int(line[2])-1
                comp_mins.append(m)
                components[m] = []
                phase[m] = {0:{},1:{}}
                mec +=int(round(float(line[10])))
            elif k[0] =='*':
                None
            else:
                line = k.split('\t')
                i = int(line[0])-1
                components[m].append(i)
                comp_dict[i] = m
                phase[m][0][i] = int(line[1])
                phase[m][1][i] = int(line[2])
        self.components = components
        self.comp_dict = comp_dict
        self.phase = phase
        self.comp_mins = comp_mins
        self.MEC = mec

    def switches(self,truegenes):
        return m_swcounter(self.phase,truegenes,self)
        
def make_golden_from_true2(filename):
    f = open(filename,'r')
    a = f.readlines()
    f.close()
    i = 0
    while a[i][0] == '#':
        i+=1
    c = [x.split()[9][:3] for x in a[i:]]
    vChroms = [x.split()[0] for x in a[i:]]  
    vPositions = [x.split()[1] for x in a[i:]]  
    d = {0:{},1:{}}
    for j in range(len(c)):
            if c[j][1] == '|':
                d[0][j] = int(c[j][0])
                d[1][j] = int(c[j][2])
            else:
                d[0][j] = '.'
                d[1][j] = '.'
    return d,vChroms,vPositions

def m_swcounter(sol,gold,G):
    switches = 0
    for start in G.comp_mins:
        num,l = switches_in_comp(start,sol,gold,G)
        switches +=num
                
    return switches

def switches_in_comp(start,sol,gold,G):
    m0,switches0 = switches_comp_strand(start,0,sol,gold,G)
    m1,switches1 = switches_comp_strand(start,1,sol,gold,G)
    if m0 <m1:
        return m0,switches0
    else:
        return m1,switches1


def switches_comp_strand(start,strand,sol,gold,G):
    m = 0
    switches = []
    for i in G.components[start]:
        if (sol[start][strand][i] == gold[start][0][i] or gold[start][0][i]=='.'):
            None
        else:
            switches.append(i)  
            m+=1
            strand = 1-strand
    return m,switches


def m_swcounter2(sol,gold,G):
    switches = 0
    for start in G.comp_mins:
        num,l = switches_in_comp2(start,sol,gold,G)
        switches +=num
                
    return switches

def switches_in_comp2(start,sol,gold,G):
    m0,switches0 = switches_comp_strand2(start,0,sol,gold,G)
    m1,switches1 = switches_comp_strand2(start,1,sol,gold,G)
    if m0 <m1:
        return m0,switches0
    else:
        return m1,switches1


def switches_comp_strand2(start,strand,sol,gold,G):
    m = 0
    switches = []
    for i in G.components[start]:
        if (sol[start][strand][i] == gold[0][i] or gold[0][i]=='.'):
            None
        else:
            switches.append(i)
            m+=1
            strand = 1-strand
    return m,switches    
