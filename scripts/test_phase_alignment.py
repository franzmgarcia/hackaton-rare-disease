"""Meaningful synthetic checks of orientation and phase parity."""
import json
from pathlib import Path
import random
import tempfile
import unittest
import mappy as mp
import pysam
import phase_align as phase

class AlignmentTests(unittest.TestCase):
    def test_strands_clipping_and_variant_coordinate(self):
        rng=random.Random(1701)
        reference=''.join(rng.choices('ACGT',k=5000))
        aligner=mp.Aligner(seq=reference,preset='sr',k=15,w=5)
        header=pysam.AlignmentHeader.from_dict({'SQ':[{'SN':'N/A','LN':5000}]})
        sequence=reference[1700:1849]
        alt=next(b for b in 'ACGT' if b!=sequence[75])
        changed=sequence[:75]+alt+sequence[76:]
        for reverse in (False,True):
            query='NNNNN'+(mp.revcomp(changed) if reverse else changed)+'NNN'
            hits=[h for h in aligner.map(query) if h.is_primary]
            self.assertEqual(len(hits),1)
            read=phase.segment(hits[0],query,'I'*len(query),'test',header,1)
            bases={p:read.query_sequence[q] for q,p in read.get_aligned_pairs(matches_only=True)}
            self.assertEqual(bases[1775],alt)
            # Minimap2 may align terminal Ns as mismatches rather than soft-clip them.
            self.assertTrue(all(base==reference[p] for p,base in bases.items() if p!=1775 and base!='N'))
            self.assertEqual(read.is_reverse,reverse)
            from types import SimpleNamespace
            fake=SimpleNamespace(strand=-1 if reverse else 1,q_st=5,q_en=154,r_st=1700,
                                 ctg='N/A',mapq=60,cigar=[(149,0)],NM=1)
            clipped=phase.segment(fake,query,'I'*len(query),'clip',header,1)
            mapped={p:clipped.query_sequence[q] for q,p in clipped.get_aligned_pairs(matches_only=True)}
            self.assertEqual(mapped[1775],alt)
            self.assertEqual(clipped.cigartuples,[(4,3 if reverse else 5),(0,149),(4,5 if reverse else 3)])

    def test_cis_trans_and_disconnected_graph(self):
        oldroot=phase.ROOT;oldcandidates=phase.CANDIDATES
        try:
            for mode in ('cis','trans','disconnected'):
                with tempfile.TemporaryDirectory(prefix='phase-graph-test-') as name:
                    phase.ROOT=Path(name);phase.CANDIDATES=(101,201)
                    h=pysam.AlignmentHeader.from_dict({'SQ':[{'SN':'chr15','LN':101991189}]})
                    vh=pysam.VariantHeader();vh.contigs.add('chr15',length=101991189)
                    vh.formats.add('GT',1,'String','Genotype');vh.add_sample('proband')
                    vcf=phase.ROOT/'input.vcf'
                    with pysam.VariantFile(vcf,'w',header=vh) as out:
                        for pos in (101,201):
                            r=out.new_record(contig='chr15',start=pos-1,alleles=('A','G'))
                            r.samples['proband']['GT']=(0,1);out.write(r)
                    with pysam.AlignmentFile(phase.ROOT/'u.bam','wb',header=h) as out:
                        for i in range(6):
                            r=pysam.AlignedSegment(h);r.query_name=f'f{i}';r.reference_id=0
                            r.reference_start=90+i;r.mapping_quality=60
                            length=130 if mode!='disconnected' else 40
                            seq=list('A'*length);seq[100-r.reference_start]='G' if i%2 else 'A'
                            if mode!='disconnected':seq[200-r.reference_start]=('G' if i%2 else 'A') if mode=='cis' else ('A' if i%2 else 'G')
                            r.query_sequence=''.join(seq);r.query_qualities=pysam.qualitystring_to_array('I'*length);r.cigartuples=[(0,length)]
                            out.write(r)
                    pysam.sort('-o',str(phase.ROOT/'target.bam'),str(phase.ROOT/'u.bam'));pysam.index(str(phase.ROOT/'target.bam'))
                    # Synthetic coordinates need a matching fetch interval.
                    oldstart,oldend=phase.START,phase.END;phase.START,phase.END=0,1000
                    try:result=phase.graph(vcf,'test')
                    finally:phase.START,phase.END=oldstart,oldend
                    self.assertEqual(result['graph_phase'],mode if mode!='disconnected' else 'unresolved')
        finally:phase.ROOT=oldroot;phase.CANDIDATES=oldcandidates

if __name__=='__main__':unittest.main()
