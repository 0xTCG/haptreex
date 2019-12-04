import read_class
import math

def find_rates(snps,reads,rate):
    ##if sufficient coverage, approximation works well
    ##snps must be sorted
    lb_for_approx = 200
    counts = get_counts(snps,reads)
    m = max(map(max,counts.values()))

    if m == 0: print("should not happen",snps,reads)
    elif m > lb_for_approx: r = approx_rate(list(counts.values()))
    else: r = NW(counts,rate)
    
    return {0:r,1:1-r} ##had swapped
    
def approx_rate(counts):
        l = len(counts)
        maxes = map(max,counts)
        sums = map(sum,counts)
        #if sum(sums) == 0:
        #        print counts
        return float(sum(maxes))/float(sum(sums))

def get_counts(snps,reads,shuffle = True):
    if shuffle:
        reads = read_class.sample_from_reads(reads)
    counts = {i:[0,0] for i in range(len(snps))}
    back = {}
    for i in range(len(snps)):
            back[snps[i]]=i
    for R in reads:
        for snp in R.keys:
            counts[back[snp]][R.read[snp]%2]+=R.count
    return counts



def der(f,x,dx):
    return (f(x+dx) - f(x-dx))/(2*float(dx))

def ld(f,x,dx):
    return f(x+dx)-f(x)

def rd(f,x,dx):
    return f(x)-f(x-dx)



def NW(vec,rate):
    left = .001
    right = 1-left
    dx = .001
    d = 1 
    def f(x):
        FD = {}
        return math.log(forward(FD,0,vec,len(vec)-1,rate,x)+forward(FD,1,vec,len(vec)-1,rate,x))

    if f(.501)==f(.499):
        if f(.5)>=f(.499):
            return .5
        else:
            left = .5 + .001
            
    while abs(d)>.001:
        mid = (right+left)/2.
        d = der(f,mid,dx)
        if d > 0:
            left = mid
        else:
            right = mid

        if mid > .98:
                return .98
        if mid < .02:
                return .02
    return mid

def trans(s1,s2,rate):
        if s1 == s2:
                return rate
        else:
                return 1-rate



def forward(FD,X,vec,i,rate,p):    
    if (X,i) in FD:
        return FD[(X,i)]
    e = .0001
    p = (p*(1-e) + (1-p)*(e))
    q = 1-p
    if i == 0:
        Y = vec[i]
        if X ==1:
            FD[(1,0)] = 0 
            return 0
        else: 
            val = choose(Y[0]+Y[1],Y[0])*(p**Y[0])*(q**Y[1])
            FD[(0,0)] = 10**10*val
            return 10**10*val
    else:
        Y=vec[i]
        f0 = forward(FD,0,vec,i-1,rate,p)  
        f1 = forward(FD,1,vec,i-1,rate,p)
        if (0,i-1) not in FD:
            FD[(0,i-1)] = (10**10)*f0
        if (1,i-1) not in FD:
            FD[(1,i-1)] = (10**10)*f1
                
        SUM_log = math.log(trans(0,X,rate)*f0 + trans(1,X,rate)*f1)
        sample_prob_log = (Y[X]*math.log(p)) + (Y[1-X]*math.log(q))
        binom_log = math.log(choose(sum(Y),Y[0]))
        f_val = math.exp(SUM_log + sample_prob_log + binom_log)
        
        if (X,i) not in FD:
            FD[(X,i)] = 10**10*f_val
        return 10**10*f_val



def choose(n,k):
    #n CHOOSE k
    j = min(k,n-k)
    numerator = 1
    for i in range(j):
        numerator = numerator * (n-i)
    return numerator//math.factorial(j)
    
