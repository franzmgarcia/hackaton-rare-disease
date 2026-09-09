"""Whole-genome alignment of retrieved pairs; local phase graphs and WhatsHap.

Run after the reference and all four lane extractions have completed. Alignment
and phase evidence remain local. No phase is inferred from lack of a graph path.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import mappy as mp
import pysam

ROOT=Path('data/phasing')
REF=Path('data/reference/hg38.fa.gz')
START,END=40180000,40245000
CANDIDATES=(40209701,40220612)

def segment(hit,seq,qual,name,header,mate_number):
    read=pysam.AlignedSegment(header)
    read.query_name=name
    read.flag=1|(64 if mate_number==1 else 128)
    if hit.strand<0:
        seq=mp.revcomp(seq);qual=qual[::-1];read.flag|=16
        left=len(seq)-hit.q_en;right=hit.q_st
    else:left=hit.q_st;right=len(seq)-hit.q_en
    read.query_sequence=seq;read.query_qualities=pysam.qualitystring_to_array(qual)
    read.reference_id=header.get_tid(hit.ctg);read.reference_start=hit.r_st
    read.mapping_quality=hit.mapq
    read.cigartuples=([(4,left)] if left else [])+[(op,n) for n,op in hit.cigar]+([(4,right)] if right else [])
    read.set_tag('NM',hit.NM);read.set_tag('RG','library1')
    assert read.infer_query_length()==len(seq)
    return read

def unclipped_end(r):
    if r.is_reverse:
        return r.reference_end+(r.cigartuples[-1][1] if r.cigartuples[-1][0]==4 else 0)
    return r.reference_start-(r.cigartuples[0][1] if r.cigartuples[0][0]==4 else 0)

def align_all(lanes):
    for lane in lanes:
        if not (ROOT/f'lane{lane}/complete.json').exists():raise RuntimeError(f'lane {lane} not yet extracted')
        if not (ROOT/f'lane{lane}/specific_complete.json').exists():raise RuntimeError(f'lane {lane} specificity extraction missing')
    if REF.stat().st_size!=983659424:raise RuntimeError('Reference download incomplete or changed: verify expected size')
    idx=ROOT/'hg38.sr.mmi'
    print('Loading/building whole-genome short-read index',flush=True)
    if idx.exists() and idx.stat().st_size>1000000:
        # Never pass the input index as fn_idx_out: minimap2 opens output first.
        aligner=mp.Aligner(str(idx),preset='sr',n_threads=3,best_n=5)
    else:
        temporary=ROOT/'hg38.sr.building.mmi'
        aligner=mp.Aligner(str(REF),preset='sr',n_threads=3,fn_idx_out=str(temporary),best_n=5)
        if aligner:temporary.replace(idx)
    if not aligner:raise RuntimeError('index build failed')
    names=aligner.seq_names
    if not {'chr1','chr15','chrX'}<=set(names):raise RuntimeError('whole-genome contigs missing')
    assert aligner.seq('chr15',START,END).upper()==json.loads((ROOT/'ucsc_region.json').read_text())['dna'].upper()
    sq=[{'SN':name,'LN':len(aligner.seq(name))} for name in names]
    header=pysam.AlignmentHeader.from_dict({'HD':{'VN':'1.6'},'SQ':sq,'RG':[{'ID':'library1','SM':'proband','LB':'library1','PL':'ILLUMINA'}]})
    reference=ROOT/'chr15.fa'
    if not reference.exists():
        seq=aligner.seq('chr15')
        with reference.open('w') as f:
            f.write('>chr15\n')
            for i in range(0,len(seq),80):f.write(seq[i:i+80]+'\n')
    pysam.faidx(str(reference))
    counts=Counter();spans=[];dedup=set()
    with pysam.AlignmentFile(str(ROOT/'target.unsorted.bam'),'wb',header=header) as out:
        for lane in lanes:
            reads1=mp.fastx_read(str(ROOT/f'lane{lane}/specific_R1.fastq'))
            reads2=mp.fastx_read(str(ROOT/f'lane{lane}/specific_R2.fastq'))
            for r1,r2 in itertools.zip_longest(reads1,reads2):
                if r1 is None or r2 is None or r1[0]!=r2[0]:raise RuntimeError('paired extraction mismatch')
                counts['retrieved_pairs']+=1
                # Independent-end mapping does not suppress long inserts based on a library-length prior.
                hits=[]
                for r in (r1,r2):
                    primary=[h for h in aligner.map(r[1]) if h.is_primary]
                    hits.append(primary[0] if len(primary)==1 else None)
                if not any(h and h.ctg=='chr15' and h.r_st<END and h.r_en>START for h in hits):continue
                counts['pairs_with_target_alignment']+=1
                name=hashlib.sha256((str(lane)+':'+r1[0]).encode()).hexdigest()[:24]
                mapped=[segment(h,r[1],r[2],name,header,i+1) if h else None for i,(h,r) in enumerate(zip(hits,(r1,r2)))]
                usable=[r for r in mapped if r is not None]
                key=tuple((r.reference_id,unclipped_end(r),r.is_reverse) if r is not None else None for r in mapped)
                duplicate=key in dedup
                dedup.add(key)
                if duplicate:counts['duplicate_endpoint_pairs']+=1
                a,b=mapped
                span=None;fr=False
                if a is not None and b is not None and a.reference_id==b.reference_id:
                    span=max(a.reference_end,b.reference_end)-min(a.reference_start,b.reference_start)
                    left,right=sorted((a,b),key=lambda r:r.reference_start)
                    fr=not left.is_reverse and right.is_reverse
                    if min(a.mapping_quality,b.mapping_quality)>=30 and not duplicate:
                        spans.append({'lane':lane,'span':span,'FR':fr,'read':name})
                for i,r in enumerate(mapped):
                    if r is None:continue
                    mate=mapped[1-i]
                    if mate is None:r.flag|=8
                    else:
                        r.next_reference_id=mate.reference_id;r.next_reference_start=mate.reference_start
                        if mate.is_reverse:r.flag|=32
                        if span is not None:
                            r.template_length=span if (r.reference_start<mate.reference_start or (r.reference_start==mate.reference_start and i==0)) else -span
                            if fr and span<=1000:r.flag|=2
                    if duplicate:r.flag|=1024
                    out.write(r)
            print(json.dumps({'lane_aligned':lane,'counts':counts}),flush=True)
    pysam.sort('-o',str(ROOT/'target.bam'),str(ROOT/'target.unsorted.bam'))
    pysam.index(str(ROOT/'target.bam'))
    vals=[s['span'] for s in spans if s['FR']]
    audit={'counts':counts,'lanes':lanes,'mappy_version':mp.__version__,'reference_contigs':len(sq),
           'FR_span_quantiles':dict(zip(['min','p50','p95','p99','max'],np.quantile(vals,[0,.5,.95,.99,1]).tolist())) if vals else {},
           'long_or_nonFR_fragments':[s for s in spans if s['span']>1000 or not s['FR']],
           'limitations':['Target-seed retrieval before whole-genome alignment','No UMI: endpoint duplicate marking approximates independent molecules','Template-length prior not imposed in end mapping','Long/chimeric pairs require review before accepting phase']}
    (ROOT/'alignment_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    return header

def local_vcf():
    source=pysam.VariantFile('data/input.vcf.gz')
    chrom='15' if '15' in source.header.contigs else 'chr15'
    header=pysam.VariantHeader();header.contigs.add('chr15',length=101991189)
    for name,number,kind in [('GT',1,'String'),('GQ',1,'Integer'),('DP',1,'Integer'),('AD','R','Integer')]:
        header.formats.add(name,number,kind,name)
    header.add_sample('proband')
    for key in source.header.filters:
        if key!='PASS':header.filters.add(key,None,None,source.header.filters[key].description or key)
    with pysam.VariantFile(str(ROOT/'original_local.vcf'),'w',header=header) as out:
        for r in source.fetch(chrom,START,END):
            nr=out.new_record(contig='chr15',start=r.start,alleles=r.alleles,qual=r.qual)
            for filt in r.filter:nr.filter.add(filt)
            sample=next(iter(r.samples.values()))
            for key in ('GT','GQ','DP','AD'):
                if key in sample:nr.samples['proband'][key]=sample[key]
            out.write(nr)

def novel_anchors():
    executable=ROOT/'bcftools-1.22/bcftools'
    with (ROOT/'local_recall.log').open('w') as log:
        subprocess.run([str(executable),'mpileup','-f',str(ROOT/'chr15.fa'),'-r',f'chr15:{START+1}-{END}',
                       '-q','30','-Q','30','-a','FORMAT/AD,FORMAT/DP','-Ob','-o',str(ROOT/'local_pileup.bcf'),str(ROOT/'target.bam')],
                       stdout=log,stderr=subprocess.STDOUT,check=True)
        subprocess.run([str(executable),'call','-m','-v','-a','GQ','-Ov','-o',str(ROOT/'local_recall.vcf'),str(ROOT/'local_pileup.bcf')],
                       stdout=log,stderr=subprocess.STDOUT,check=True)
    original=pysam.VariantFile(str(ROOT/'original_local.vcf'))
    header=original.header.copy();records=list(original);occupied={r.pos for r in records};new=[]
    for r in pysam.VariantFile(str(ROOT/'local_recall.vcf')):
        s=r.samples['proband'];gt=s.get('GT');ad=s.get('AD')
        if r.pos in occupied or not r.alts or len(r.alts)!=1:continue
        if any(len(a)>50 or not set(a)<=set('ACGT') for a in r.alleles):continue
        if gt not in ((0,1),(1,0)) or not ad or sum(ad)<10 or min(ad)<4 or (r.qual or 0)<30:continue
        if not .2<=ad[1]/sum(ad)<=.8:continue
        nr=header.new_record(contig='chr15',start=r.start,alleles=r.alleles,qual=r.qual)
        for key in ('GT','GQ','DP','AD'):
            if key in s:nr.samples['proband'][key]=s[key]
        new.append(nr)
    with pysam.VariantFile(str(ROOT/'augmented_local.vcf'),'w',header=header) as out:
        for r in sorted(records+new,key=lambda r:r.pos):out.write(r)
    (ROOT/'new_anchor_candidates.json').write_text(json.dumps([{'pos':r.pos,'ref':r.ref,'alt':r.alts[0],'AD':r.samples['proband']['AD']} for r in new],indent=2)+'\n')
    # New SNVs/indels are exploratory anchors: balanced reads do not establish a true variant or phase.

def graph(vcf,label,mapq=30,baseq=30):
    hets={r.pos:(r.ref,r.alts[0]) for r in pysam.VariantFile(str(vcf)) if len(r.ref)==1 and r.alts and len(r.alts)==1 and len(r.alts[0])==1 and r.samples['proband'].get('GT') in ((0,1),(1,0))}
    fragments=defaultdict(list);observations=defaultdict(Counter)
    with pysam.AlignmentFile(str(ROOT/'target.bam')) as bam:
        for r in bam.fetch('chr15',START,END):
            if r.is_duplicate or r.is_secondary or r.is_supplementary or r.mapping_quality<mapq:continue
            calls={}
            for q,p in r.get_aligned_pairs(matches_only=True):
                if p+1 not in hets or r.query_qualities[q]<baseq:continue
                alleles=hets[p+1];base=r.query_sequence[q]
                if base in alleles:calls[p+1]=alleles.index(base)
            fragments[r.query_name].append({'calls':calls,'start':r.reference_start,'end':r.reference_end,'reverse':r.is_reverse,'mate':1 if r.is_read1 else 2})
    edges=defaultdict(Counter);direct=[];overlap_conflicts=0
    for name,reads in fragments.items():
        joined=defaultdict(set)
        for r in reads:
            for p,a in r['calls'].items():joined[p].add(a)
        overlap_conflicts+=sum(len(v)>1 for v in joined.values())
        calls={p:next(iter(v)) for p,v in joined.items() if len(v)==1}
        for p,a in calls.items():observations[p][a]+=1
        for p,q in itertools.combinations(sorted(calls),2):
            genotype=f'{calls[p]}{calls[q]}'
            edges[(p,q)][genotype]+=1
            if (p,q)==CANDIDATES:direct.append({'fragment':name,'alleles':genotype,'reads':reads})
    rows=[{'positions':key,'allele_pairs':dict(value),'same':value['00']+value['11'],'opposite':value['01']+value['10']} for key,value in sorted(edges.items())]
    any_adjacency=defaultdict(set)
    for row in rows:
        a,b=row['positions'];any_adjacency[a].add(b);any_adjacency[b].add(a)
    any_visited={CANDIDATES[0]};any_todo=[CANDIDATES[0]]
    while any_todo:
        for b in any_adjacency[any_todo.pop()]:
            if b not in any_visited:any_visited.add(b);any_todo.append(b)
    # Conservative graph: minimum 3 independent endpoint-deduplicated fragments, no discordant high-quality evidence.
    adjacency=defaultdict(list)
    for row in rows:
        a,b=row['positions'];same,opp=row['same'],row['opposite']
        if max(same,opp)>=3 and min(same,opp)==0:
            parity=int(opp>same);adjacency[a].append((b,parity));adjacency[b].append((a,parity))
    visited={CANDIDATES[0]:0};todo=[CANDIDATES[0]];conflicts=[]
    while todo:
        a=todo.pop()
        for b,parity in adjacency[a]:
            state=visited[a]^parity
            if b in visited:
                if visited[b]!=state:conflicts.append((a,b))
            else:visited[b]=state;todo.append(b)
    reachable=CANDIDATES[1] in visited and not conflicts
    out={'label':label,'MAPQ_min':mapq,'base_quality_min':baseq,'heterozygous_SNVs':len(hets),'candidate_observations':{str(p):dict(observations[p]) for p in CANDIDATES},'edges':rows,
         'direct_candidate_fragments':direct,'overlapping_mate_conflicts':overlap_conflicts,'reachable_from_L737Ter':sorted(visited),
         'any_edge_candidate_connected':CANDIDATES[1] in any_visited,'any_edge_reachable_from_L737Ter':sorted(any_visited),
         'cycle_conflicts':conflicts,'candidate_connected':reachable,'graph_phase':('trans' if visited[CANDIDATES[1]] else 'cis') if reachable else 'unresolved',
         'limitations':['SNV graph only; thresholds recorded separately','WhatsHap separately includes indels','Edge threshold is exploratory, not calibrated phase confidence','Long inserts/chimeric fragments need inspection even if graph connected']}
    (ROOT/f'{label}_phase_graph.json').write_text(json.dumps(out,indent=2)+'\n')
    return out

def whatshap(vcf,label,mapq=30):
    command=[str(Path(sys.executable).parent/'whatshap'),'phase','--reference',str(ROOT/'chr15.fa'),'--mapping-quality',str(mapq),'--ignore-read-groups','--output-read-list',str(ROOT/f'{label}_whatshap_reads.txt'),'-o',str(ROOT/f'{label}_whatshap.vcf'),str(vcf),str(ROOT/'target.bam')]
    with (ROOT/f'{label}_whatshap.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    result=[]
    for r in pysam.VariantFile(str(ROOT/f'{label}_whatshap.vcf')):
        if r.pos in CANDIDATES:
            s=r.samples['proband'];result.append({'pos':r.pos,'GT':s['GT'],'phased':s.phased,'PS':s.get('PS')})
    (ROOT/f'{label}_whatshap_candidates.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--skip-alignment',action='store_true')
    parser.add_argument('--lanes',type=int,nargs='+',default=list(range(1,5)));args=parser.parse_args()
    if not args.skip_alignment:align_all(args.lanes)
    elif json.loads((ROOT/'alignment_audit.json').read_text())['lanes']!=args.lanes:
        raise RuntimeError('Requested lanes differ from the existing BAM audit')
    local_vcf();novel_anchors()
    for label,vcf in [('original',ROOT/'original_local.vcf'),('augmented',ROOT/'augmented_local.vcf'),('recall',ROOT/'local_recall.vcf')]:
        result=graph(vcf,label);whatshap(vcf,label)
        print(json.dumps({'analysis':label,'phase':result['graph_phase'],'direct_fragments':len(result['direct_candidate_fragments'])}),flush=True)
    (ROOT/'analysis_complete.json').write_text(json.dumps({'completed':True,'lanes':args.lanes,'all_lanes':sorted(args.lanes)==[1,2,3,4],'phase_call_requires_review_of_graph_and_whatshap':True})+'\n')
