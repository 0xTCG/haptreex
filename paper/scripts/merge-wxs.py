import sys
from pprint import pprint

phases = {}
lines = []
with open(sys.argv[1]) as vcf: # RNAseq
    for line in vcf:
        if line[0] == '#':
            continue

        chr, pos, id, ref, alt, _, _, _, fmt_, sample = line.split('\t')
        if len(ref) != 1 or len(alt) != 1:  # Ignore indels
            continue
        alleles = [ref] + alt.split(',')
        if len(alleles) > 2: # only diploids
            continue
        if sample.startswith('0/1'):
            phases[chr, pos] = ref, alt, line

errors = []
with open(sys.argv[2]) as vcf: # WXS
    for line in vcf:
        if line[0] == '#':
            print(line.strip())
            continue

        chr, pos, id, ref, alt, _, _, _, fmt_, sample = line.split('\t')
        chr = chr[3:]
        if chr == 'M': chr = 'MT'
        if len(ref) != 1 or len(alt) != 1:  # Ignore indels
            continue
        alleles = [ref] + alt.split(',')
        if len(alleles) > 2: # only diploids
            continue
        if not sample.startswith('0/1'):
            continue
        if (chr, pos) in phases:
            r, a, p = phases[chr, pos]
            if r != ref or a != alt:
                errors.append((chr, pos))
        else:
            if chr == 'MT':
                phases[chr, pos] = ref, alt, 'MT'+line[4:]
            else:
                phases[chr, pos] = ref, alt, line[3:]


CC=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT', 'GL000207.1', 'GL000226.1', 'GL000229.1', 'GL000231.1', 'GL000210.1', 'GL000239.1', 'GL000235.1', 'GL000201.1', 'GL000247.1', 'GL000245.1', 'GL000197.1', 'GL000203.1', 'GL000246.1', 'GL000249.1', 'GL000196.1', 'GL000248.1', 'GL000244.1', 'GL000238.1', 'GL000202.1', 'GL000234.1', 'GL000232.1', 'GL000206.1', 'GL000240.1', 'GL000236.1', 'GL000241.1', 'GL000243.1', 'GL000242.1', 'GL000230.1', 'GL000237.1', 'GL000233.1', 'GL000204.1', 'GL000198.1', 'GL000208.1', 'GL000191.1', 'GL000227.1', 'GL000228.1', 'GL000214.1', 'GL000221.1', 'GL000209.1', 'GL000218.1', 'GL000220.1', 'GL000213.1', 'GL000211.1', 'GL000199.1', 'GL000217.1', 'GL000216.1', 'GL000215.1', 'GL000205.1', 'GL000219.1', 'GL000224.1', 'GL000223.1', 'GL000195.1', 'GL000212.1', 'GL000222.1', 'GL000200.1', 'GL000193.1', 'GL000194.1', 'GL000225.1', 'GL000192.1', 'NC_007605', 'hs37d5'] 
def ff(x):
    return (CC.index(x[0]), int(x[1]))

for i in sorted(phases.keys(), key=ff):
    print(phases[i][2].strip())

sys.stderr.write(f'total errors: {len(errors)}: {errors}\n')
