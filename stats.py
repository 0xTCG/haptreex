import sys, os, gzip
from pprint import pprint

def parse_phaser(phaser_path, vcf):
    phases = {}
    blocks, block_snps, edges, errors = 0, 0, 0, 0
    st, spans = 0, []
    prev, prev_vcf = None, None
    prev_pi = ''
    with gzip.open(phaser_path+".vcf.gz", 'rt') as hap:
        for line in hap:
            if line[0] == '#':
                continue

            chr, pos, _, _, _, _, _, _, fmt_, sample = line.split('\t')
            fmt = dict(zip(fmt_.split(':'), sample.split(':')))
            if 'PG' not in fmt or 'PI' not in fmt or '|' not in fmt['PG']:
                continue
            
            pg, pi = tuple(fmt['PG'].split('|')), fmt['PI']
            if pi != prev_pi:
                if prev_pi != '':
                    edges += block_snps - 1
                    spans.append((st, prev))

                blocks += 1
                block_snps = 0
                prev, prev_vcf = None, None
                prev_pi = pi
            
            gt0, gt1 = pg
            chr = chr[3:] if chr.startswith('chr') else chr
            k = chr, pos
            gt = (gt0, gt1)
            assert k not in phases
            phases[k] = gt
            block_snps += 1

            if not prev:
                st = int(pos)
            prev = int(pos)

            if k in vcf:
                if not prev_vcf:
                    prev_vcf = k
                else:
                    if vcf[k] != vcf[prev_vcf] and phases[k] == phases[prev_vcf]:
                        # print ( 'ERR', chr, pos )
                        errors += 1
                    elif vcf[k] == vcf[prev_vcf] and phases[k] != phases[prev_vcf]:
                        # print ( 'ERR', chr, pos )
                        errors += 1
                    prev_vcf = k
    edges += block_snps - 1
    spans.append((st, prev))
    blocks += 1

    spans = sorted(spans)
    for i in range(len(spans)):
        if not spans[i]: continue
        s, e = spans[i]
        for j in range(i + 1, len(spans)):
            if not spans[j]: continue
            sx, ex = spans[j]
            if sx > e: break
            if ex <= e: spans[j] = None
            elif sx <= e: spans[j] = (e+1, ex)
    if spans:
        span = spans[0][1] - spans[0][0] + 1
        prev_i = 0
        for i in range(1, len(spans)):
            assert (not spans[i]) or (spans[i][0] > spans[prev_i][1])
            if spans[i]:
                span += spans[i][1] - spans[i][0] + 1
                prev_i=i
    else:
        span = 0
    return phases, edges, blocks, errors, span

def parse_haptree(haptree_path, vcf):
    phases = {}
    blocks, block_snps, edges, errors = 0, 0, 0, 0
    st, spans = 0, []
    prev, prev_vcf = None, None
    with open(haptree_path) as hap:
        for line in hap:
            if line[:5] == 'BLOCK':
                blocks += 1
                block_snps = 0
                prev, prev_vcf = None, None
            elif line[:5] == '*****':
                edges += block_snps - 1
                spans.append((st, prev))
            else:
                _, gt0, gt1, chr, pos, *_ = line.split('\t')
                chr = chr[3:] if chr.startswith('chr') else chr
                k = chr, pos
                gt = (gt0, gt1)
                assert k not in phases
                phases[k] = gt
                block_snps += 1

                if not prev:
                    st = int(pos)
                prev = int(pos)

                if k in vcf:
                    if not prev_vcf:
                        prev_vcf = k
                    else:
                        if vcf[k] != vcf[prev_vcf] and phases[k] == phases[prev_vcf]:
                            # print ( 'ERR', chr, pos )
                            errors += 1
                        elif vcf[k] == vcf[prev_vcf] and phases[k] != phases[prev_vcf]:
                            # print ( 'ERR', chr, pos )
                            errors += 1
                        prev_vcf = k

    spans = sorted(spans)
    for i in range(len(spans)):
        if not spans[i]: continue
        s, e = spans[i]
        for j in range(i + 1, len(spans)):
            if not spans[j]: continue
            sx, ex = spans[j]
            if sx > e: break
            if ex <= e: spans[j] = None
            elif sx <= e: spans[j] = (e+1, ex)
    if spans:
        span = spans[0][1] - spans[0][0] + 1
        prev_i = 0
        for i in range(1, len(spans)):
            assert (not spans[i]) or (spans[i][0] > spans[prev_i][1])
            if spans[i]:
                span += spans[i][1] - spans[i][0] + 1
                prev_i=i
    else:
        span = 0
    return phases, edges, blocks, errors, span

def parse_vcf(vcf_path: str, phases = dict()):
    phases = {}
    with open(vcf_path) as vcf:
        for line in vcf:
            if line[0] == '#':
                continue

            chr, pos, id, ref, alt, _, _, _, fmt_, sample = line.split('\t')
            if len(ref) != 1 or len(alt) != 1:  # Ignore indels
                continue
            fmt = dict(zip(fmt_.split(':'), sample.split(':')))
            if '|' not in fmt['GT']: # not phased
                continue
            gt = tuple(fmt['GT'].split('|'))
            alleles = [ref] + alt.split(',')
            if len(alleles) > 2: # only diploids
                continue
            # Get only alleles that are specified in GT field
            a = [a for i, a in enumerate(alleles) if str(i) in gt and len(a) == 1]
            if len(a) <= 1:  # Ignore homozygous SNPs
                continue
            # gt = ''.join(gt if not invert else [{'0':'1','1':'0'}[g] for g in gt])
            k = chr, pos
            assert k not in phases
            phases[k] = gt
    return phases

if __name__ == "__main__":
    vcf_path = sys.argv[1]
    hap_path = sys.argv[2]
    vcf = parse_vcf(vcf_path)
    if hap_path[-6:] == 'phaser':
        hap, edges, blocks, errors, span = parse_phaser(hap_path, vcf)
    else:
        hap, edges, blocks, errors, span = parse_haptree(hap_path, vcf)
    print(','.join(map(str, [
        os.path.basename(hap_path),
        len(hap),
        errors, 
        100*errors/len(hap) if len(hap) else 0,
        edges,
        blocks,
        100*errors/edges,
        span])))


# # print(f"""stats for {hap_path} (vcf: {vcf_path}):
# #   snps{len(hap)}
# #   errors: {errors}
# #   error%: {100*errors/len(hap):.1f}%
# #   edges: {edges}
# #   blocks: {blocks}
# #   rate%: {100*errors/edges:.3f}%
# #   span: {span}
# # """)
