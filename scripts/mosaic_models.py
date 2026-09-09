"""Exploratory unphased allele-count model comparison; not a validated mCA caller."""
import os
os.environ.setdefault('MPLCONFIGDIR',str(__import__('pathlib').Path('data/mosaic/mpl').resolve()))
import json
from pathlib import Path
import numpy as np
from scipy.special import betaln,gammaln
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path('data/mosaic')

def histogram(k,n):
    pairs,counts=np.unique(np.column_stack([k,n]),axis=0,return_counts=True)
    return pairs[:,0].astype(float),pairs[:,1].astype(float),counts

def logbb(k,n,p,rho):
    t=(1-rho)/rho
    return gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+betaln(k+p*t,n-k+(1-p)*t)-betaln(p*t,(1-p)*t)

def logprob(h,theta,mixture):
    k,n,_=h; p,rho=theta[:2]
    if not mixture:return logbb(k,n,p,rho)
    d=theta[2]
    return np.logaddexp(logbb(k,n,p-d,rho),logbb(k,n,p+d,rho))-np.log(2)

def fit(h,mixture=False):
    starts=[[.5,.01]] if not mixture else [[.5,.01,d] for d in (0,.03,.08,.15)]
    bounds=[(.4,.6),(.0001,.15)]+([(0,.25)] if mixture else [])
    fits=[minimize(lambda theta:-np.dot(h[2],logprob(h,theta,mixture)),start,bounds=bounds,method='L-BFGS-B') for start in starts]
    successful=[r for r in fits if r.success and np.isfinite(r.fun)]
    if not successful: raise RuntimeError('No model fit converged')
    return min(successful,key=lambda r:r.fun)

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,default=OUT)
    parser.add_argument('--no-plot',action='store_true')
    args=parser.parse_args(); OUT=args.directory
    data=np.load(OUT/'heterozygous_readouts.npz')
    c=data['chrom']; pos=data['pos']; k=data['alt_reads'].astype(float); n=k+data['ref_reads']
    fraction=k/n; records=[]; bins=[]
    for chrom in range(1,23):
        sel=c==chrom; train=sel&(((pos-1)//1000000)%2==0); test=sel&~train
        h=histogram(k[train],n[train]); hold=histogram(k[test],n[test])
        null=fit(h); mix=fit(h,True)
        scores=logprob(hold,mix.x,True)-logprob(hold,null.x,False)
        record=dict(chrom=chrom,sites=int(sel.sum()),train_sites=int(train.sum()),test_sites=int(test.sum()),
            pooled_alt_fraction=float(k[sel].sum()/n[sel].sum()),median_AD_depth=float(np.median(n[sel])),
            null_p=float(null.x[0]),null_rho=float(null.x[1]),mixture_p=float(mix.x[0]),
            mixture_rho=float(mix.x[1]),mixture_delta=float(mix.x[2]),
            train_loglik_gain=float(null.fun-mix.fun),holdout_loglik_gain=float(np.dot(hold[2],scores)))
        records.append(record)
        p=null.x[0]
        for b in np.unique((pos[sel]-1)//5000000):
            ix=sel&((pos-1)//5000000==b)
            if ix.sum()<100:continue
            excess=(fraction[ix]-p)**2-p*(1-p)/n[ix]
            bins.append(dict(chrom=chrom,start=int(b)*5000000+1,sites=int(ix.sum()),
                excess_allele_fraction_variance=float(excess.mean()),pooled_alt_fraction=float(k[ix].sum()/n[ix].sum())))
    (OUT/'chromosome_model_comparison.json').write_text(json.dumps(records,indent=2)+'\n')
    (OUT/'allele_balance_5mb_bins.json').write_text(json.dumps(bins,indent=2)+'\n')
    if args.no_plot:
        print(json.dumps(dict(sites=len(k),strongest_mixture_holdout=sorted(records,key=lambda r:r['holdout_loglik_gain'],reverse=True)[:3]),indent=2))
        raise SystemExit(0)
    depth=json.loads((OUT/'variant_site_depth_bins.json').read_text())
    medians=[r['median_variant_site_dp'] for r in depth if r['variant_sites']>=100]
    baseline=float(np.median(medians))
    fig,axes=plt.subplots(3,1,figsize=(14,10),layout='constrained')
    offsets={};cursor=0;ticks=[]
    for chrom in range(1,23):
        offsets[chrom]=cursor
        length=int(pos[c==chrom].max());ticks.append(cursor+length/2)
        cursor+=length+5000000
    x=np.array([offsets[int(a)] for a in c])+pos
    rng=np.random.default_rng(20260907)
    ids=rng.choice(len(x),size=min(150000,len(x)),replace=False)
    axes[0].scatter(x[ids]/1e6,fraction[ids],s=.15,alpha=.12,color='#174f73',rasterized=True)
    axes[0].axhline(.5,color='black',lw=.6)
    axes[0].set(ylabel='Alternate-read fraction',ylim=(0,1),title='Patient VCF: unphased allele counts (selected heterozygous SNVs)')
    xx=[(offsets[r['chrom']]+r['start']+2500000)/1e6 for r in bins]
    axes[1].scatter(xx,[r['excess_allele_fraction_variance'] for r in bins],s=5,color='#aa5432')
    axes[1].axhline(0,color='black',lw=.6)
    axes[1].set(ylabel='Excess fraction variance',title='5 Mb windows: observed variance minus binomial sampling variance; technical + biological effects')
    good=[r for r in depth if r['variant_sites']>=100]
    axes[2].scatter([(offsets[r['chrom']]+r['start']+500000)/1e6 for r in good],
        [np.log2(r['median_variant_site_dp']/baseline) for r in good],s=4,color='#4f7956')
    axes[2].axhline(0,color='black',lw=.6)
    axes[2].set(ylabel='log2 variant-site DP ratio',xlabel='Autosome',title='1 Mb windows: biased depth proxy, NOT corrected coverage or a copy-number call')
    for ax in axes:
        ax.set_xticks(np.asarray(ticks)/1e6,[str(i) for i in range(1,23)])
        ax.spines[['top','right']].set_visible(False)
    fig.savefig(OUT/'genomic_readouts.png',dpi=180)
    print(json.dumps(dict(sites=len(k),variant_site_depth_baseline=baseline,
        strongest_mixture_holdout=sorted(records,key=lambda r:r['holdout_loglik_gain'],reverse=True)[:5]),indent=2))
