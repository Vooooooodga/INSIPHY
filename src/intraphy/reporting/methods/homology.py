"""Show the correspondence algorithm and uncertainty rather than only its result."""
from .drawing import Canvas, INK, MUTED, TEAL, PURPLE, AMBER, GREY, PALE, LINE

def draw(path):
    p=Canvas('How are homologous regions identified?',
        'Exon boundaries may differ. Compare sequences first, then evaluate their genomic positions.',height=1315)
    p.section(48,142,'1','Generate and project sequence alignments')
    p.text(82,204,'Protein-coding regions',19,bold=True)
    p.lines(82,237,['MAFFT protein alignment', 'Project residues through codons to genomic coordinates'],16)
    for i,(label,seq) in enumerate([('Species A','K L V A D  E F G'),('Species B','K L I A D  E F G')]):
        y=310+i*58;p.text(87,y,label,17);p.text(230,y,seq,21,mono=True)
        p.line(386,y-23,386,y+6,stroke=PURPLE,width=3)
    p.text(86,414,'Aligned residues establish a common coordinate system.',16,fill=MUTED)
    p.line(787,190,787,425,stroke=LINE,width=1)
    p.text(833,204,'Genomic DNA, including unannotated intervals',19,bold=True)
    p.lines(833,237,['minimap2 finds candidate nucleotide matches.',
                      'Unique flanking matches delimit local realignment.'],16)
    for i,lab in enumerate(['Annotated region','Target genome']):
        yy=309+i*65;p.text(833,yy+17,lab,16)
        p.line(995,yy+12,1320,yy+12,stroke=INK,width=1)
        for x,w,c in [(1007,60,GREY),(1100,110,TEAL),(1248,62,GREY)]:p.rect(x,yy,w,23,fill=c)
    p.ribbon(1100,1210,332,1100,1210,374)
    p.text(833,424,'Bounded short alignments use affine-gap dynamic programming.',16,fill=MUTED)
    p.line(48,455,1630,455,stroke=LINE,width=1)
    p.section(48,499,'2','Select compatible matches; retain unresolved alternatives')
    p.text(82,545,'Complementary coverage: a split or fused exon',19,bold=True)
    for x,w,c in [(114,95,TEAL),(209,95,PURPLE)]:p.rect(x,585,w,25,fill=c,stroke='none')
    p.rect(114,585,190,25,fill='none')
    p.ribbon(114,209,610,114,209,685,TEAL)
    p.ribbon(209,304,610,243,338,685,PURPLE)
    p.rect(114,685,95,25,fill=TEAL);p.rect(243,685,95,25,fill=PURPLE)
    p.line(209,698,243,698)
    p.lines(375,600,['Distinct, ordered parts of the same region.',
                      'The same correspondence works in reverse:',
                      'split versus fused is not a direction of evolution.'],16)
    p.text(850,545,'Repeated coverage: do not force a one-to-one match',19,bold=True)
    p.rect(897,585,112,25,fill=TEAL)
    p.ribbon(897,1009,610,847,959,685,TEAL)
    p.ribbon(897,1009,610,1037,1149,685,TEAL)
    p.rect(847,685,112,25,fill=TEAL);p.rect(1037,685,112,25,fill=TEAL)
    p.text(1210,622,'Alternative locations',17,fill=AMBER,bold=True)
    p.text(1210,650,'remain ambiguous.',17,fill=AMBER)
    p.rect(82,745,1516,205,fill=PALE,stroke='none',r=5)
    p.text(104,779,'Ordered-chain dynamic programming',19,bold=True)
    # A repeat supplies alternative paths, not two independent ordered matches.
    for y in (831,904):
        p.line(210,867,420,y,stroke=TEAL,arrow=True)
        p.line(570,y,785,867,stroke=TEAL,arrow=True)
    for x,y,w,label in [(114,852,116,'Left flank'),(417,816,160,'Match 1'),
                         (417,889,160,'Match 2'),(785,852,126,'Right flank')]:
        p.rect(x,y,w,30,fill='white',stroke=TEAL,r=3)
        p.text(x+w/2,y+21,label,16,anchor='middle')
    p.lines(980,821,['Connect matches with compatible order and strand.',
                     'Find the chain with the largest sum of match scores.',
                     'Keep near-optimal alternatives: no forced copy choice.'],16,step=32)
    p.line(48,978,1630,978,stroke=LINE,width=1)
    p.section(48,1020,'3','When an exon is missing from the annotation')
    rows=[(90,['DNA match in the expected interval'],
           ['Sequence present; exonic status unknown*.', '* No informative transcript annotation.'],1),
          (633,['Deletion supported', 'by flanks and alignment'],
           ['Sequence absent; exonic status inapplicable.'],0),
          (1170,['Assembly gap or competing matches'],
           ['Unknown; no confirmed loss.'],None)]
    for x,title,result,state in rows:
        p.lines(x,1066,title,15,bold=True,step=20)
        p.local(x+75,1122,'sequence',state,w=220)
        p.lines(x,1190,result,14,step=25,fill=AMBER if state is None else MUTED)
    return p.write(path)
