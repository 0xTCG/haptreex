import sys, os, gzip
from pprint import pprint

## to use
#python3 stats.py path_to_vcf_for_ground_truth path_to_hc/htx/phaser_output [path_to_hc/htx/phaser_output]


def inHLA(k):
    if k[0] == "6":
        if 29691116 <= int(k[1]) and int(k[1]) <= 33054976:
            return True
    return False

def parse_phaser(phaser_path, vcf):
    phases = {}
    blocks, edges, errors, pot_err = 0, 0, 0, 0
    st = 0
    chr_prev, prev, prev_vcf = 0, None, None
    prev_vcfs = [None] #None, then last SNP in vcf for each phase
    spans = [None] #None, then span for each phase
    chr_to_phase = {} #chr to phase mapping
    max_pi = 0
    HLA = [0,0,0] #SNPs, blocks, edges
    wasHLA = False
    with gzip.open(phaser_path+".vcf.gz", 'rt') as hap:
        for line in hap:
            if line[0] == '#':
                continue

            chr, pos, _, _, _, _, _, _, fmt_, sample = line.split('\t')
            fmt = dict(zip(fmt_.split(':'), sample.split(':')))
            if 'PG' not in fmt or 'PI' not in fmt or '|' not in fmt['PG']:
                continue
            
            pg, pi = tuple(fmt['PG'].split('|')), int(fmt['PI'])

            #make phases
            gt0, gt1 = pg
            chr = chr[3:] if chr.startswith('chr') else chr
            k = chr, pos
            gt = (gt0, gt1)
            assert k not in phases
            phases[k] = gt

            if inHLA(k):
                HLA[0] += 1
                if wasHLA:
                    HLA[2] += 1
                else:
                    wasHLA = True
                    HLA[1] += 1

            #count edges and blocks and make spans and count errors
            if pi < max_pi: #interleaved case; update span and errors
                edges += 1

                #update span
                if spans[pi][0] == None: #when the case is skipped
                    blocks += 1
                    edges -= 1
                    spans[pi] = (int(pos),int(pos))
                else:
                    spans[pi] = (spans[pi][0],int(pos))

                old_k = prev_vcfs[pi] #last snp location of the phase in the vcf

                if k in vcf: #check for errors
                    if old_k:
                        if vcf[k] != vcf[old_k] and phases[k] == phases[old_k]:
                            # print ( 'ERR', chr, pos )
                            errors += 1
                        elif vcf[k] == vcf[old_k] and phases[k] != phases[old_k]:
                            # print ( 'ERR', chr, pos )
                            errors += 1
                        pot_err += 1
                    prev_vcfs[pi] = k


            elif pi == max_pi: #continuing old phase
                edges += 1

                if k in vcf:
                    if prev_vcf:
                        if vcf[k] != vcf[prev_vcf] and phases[k] == phases[prev_vcf]:
                            # print ( 'ERR', chr, pos )
                            errors += 1
                        elif vcf[k] == vcf[prev_vcf] and phases[k] != phases[prev_vcf]:
                            # print ( 'ERR', chr, pos )
                            errors += 1
                        pot_err += 1
                    prev_vcf = k

                prev = int(pos)

            elif pi > max_pi: 
                if max_pi != 0: #conclude old phase
                    assert len(spans) == max_pi
                    spans.append((st, prev))
                    #print("added", max_pi)
                    assert len(prev_vcfs) == max_pi
                    prev_vcfs.append(prev_vcf)

                blocks += 1 #beginning new phase
                wasHLA = False
                st = int(pos)
                #print("beginning", pi, "last", max_pi)
                chr_prev, prev, prev_vcf = chr, int(pos), None
                if k in vcf:
                    prev_vcf = k
                if chr in chr_to_phase:
                    chr_to_phase[chr][1] = pi
                else: 
                    chr_to_phase[chr] = [pi,pi]

                #when phaser is poorly ordered; add in the skipped phases
                spans += [(None, None)] * (pi - max_pi -1) 
                prev_vcfs += [None] * (pi - max_pi -1)

                max_pi = pi

    assert len(spans) == max_pi
    spans.append((st, prev))
    #print("added", max_pi)
    assert len(prev_vcfs) == max_pi
    prev_vcfs.append(prev_vcf)

    assert blocks == max_pi


    
    span = 0
    for chr in chr_to_phase:
        # print("chr_to_phase[chr]", chr_to_phase[chr])
        # print(spans[chr_to_phase[chr][0]], spans[chr_to_phase[chr][1]+1])
        # print(spans[chr_to_phase[chr][0]:chr_to_phase[chr][1]+1])
        span += span_from_spans(spans[chr_to_phase[chr][0]:chr_to_phase[chr][1]+1])

    return phases, edges, blocks, errors, pot_err, span, HLA

def parse_haptree(haptree_path, vcf):
    phases = {}
    blocks, block_snps, edges, errors, pot_err = 0, 0, 0, 0, 0
    st, spans = 0, {}
    prev, prev_vcf = None, None
    HLA = [0,0,0] #SNPs, blocks, edges
    wasHLA = False
    with open(haptree_path) as hap:
        for line in hap:
            if line[:5] == 'BLOCK':
                blocks += 1
                block_snps = 0
                prev, prev_vcf = None, None
                wasHLA = False
            elif line[:5] == '*****':
                edges += block_snps - 1
                if chr in spans:
                    spans[chr].append((st, prev))
                else:
                    spans[chr] = [(st, prev)]
            else:
                _, gt0, gt1, chr, pos, *_ = line.split('\t')
                chr = chr[3:] if chr.startswith('chr') else chr
                k = chr, pos
                gt = (gt0, gt1)
                assert k not in phases
                phases[k] = gt
                block_snps += 1

                if inHLA(k):
                    HLA[0] += 1
                    if wasHLA:
                        HLA[2] += 1
                    else:
                        wasHLA = True
                        HLA[1] += 1



                if not prev:
                    st = int(pos)
                prev = int(pos)

                if k in vcf:
                    if prev_vcf:
                        if vcf[k] != vcf[prev_vcf] and phases[k] == phases[prev_vcf]:
                            # print ( 'ERR', chr, pos )
                            errors += 1
                        elif vcf[k] == vcf[prev_vcf] and phases[k] != phases[prev_vcf]:
                            # print ( 'ERR', chr, pos )
                            errors += 1
                        pot_err += 1
                    prev_vcf = k

    if line[:5] != '*****': #hc output, finish processing
        edges += block_snps - 1
        if chr in spans:
            spans[chr].append((st, prev))
        else:
            spans[chr] = [(st, prev)]
    span = 0
    for chr in spans:
        span += span_from_spans(spans[chr])

    return phases, edges, blocks, errors, pot_err, span, HLA

def span_from_spans(chr_spans):
    spans = sorted(chr_spans)
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

    ##### Equivalent way to calculate span
    # sorted_chr_spans = sorted(chr_spans)
    # new_span = 0
    # s, e = sorted_chr_spans[0]
    # for i in sorted_chr_spans[1:]:
    #     sx, ex = i
    #     if sx > e: 
    #         new_span += e - s + 1
    #         s, e = sx, ex
    #     else:
    #         if ex > e: e = ex
    # new_span += e - s + 1

    # assert span == new_span

    return span


def parse_vcf(vcf_path: str, phases = dict()):
    phases = {}
    with open(vcf_path) as vcf:
    #with gzip.open(vcf_path, 'rt') as vcf:
        for line in vcf:
            if line[0] == '#':
                continue

            chr, pos, id, ref, alt, _, _, _, fmt_, sample = line.split('\t')
            chr = chr[3:] if chr.startswith('chr') else chr
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
    vcf = parse_vcf(vcf_path)
    for hap_path in sys.argv[2:]:
        if hap_path[-6:] == 'phaser':
            hap, edges, blocks, errors, pot_err, span, HLA = parse_phaser(hap_path, vcf)
        else:
            hap, edges, blocks, errors, pot_err, span, HLA = parse_haptree(hap_path, vcf)
        # print(','.join(map(str, [
        #     #os.path.basename(hap_path),
        #     HLA[0], HLA[1], HLA[2]
        #     ])))
        print(','.join(map(str, [
            os.path.basename(hap_path),
            len(hap),
            errors, 
            pot_err,
            edges,
            blocks,
            100*errors/pot_err,
            span])))
    print(sys.argv[2:])


# # print(f"""stats for {hap_path} (vcf: {vcf_path}):
# #   snps{len(hap)}
# #   errors: {errors}
# #   error%: {100*errors/len(hap):.1f}%
# #   edges: {edges}
# #   blocks: {blocks}
# #   rate%: {100*errors/edges:.3f}%
# #   span: {span}
# # """)
