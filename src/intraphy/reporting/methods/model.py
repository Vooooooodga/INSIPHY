"""Biological character coding, a two-state process, and tree likelihood."""
from .drawing import Canvas, INK, MUTED, TEAL, PURPLE, AMBER, GREY, PALE, LINE
from .example import CHARACTERS, probability_example

def draw(path):
    p=Canvas('From gene structures to an evolutionary model',
        'Each aligned region defines a character. Its states are observed at species and inferred at ancestral nodes.',height=1140)
    p.section(48,141,'1','Encode comparable biological states')
    p.text(74,203,'Species tree and gene structures',19,bold=True)
    xy=p.tree(80,290,width=110,step=62)
    for s in 'ABCD':p.gene(241,xy[s][1]-11,width=337,segment=s!='D',exonic=None if s=='C' else True,intron=s!='B')
    p.line(610,381,670,381,arrow=True)
    for x,labels,kind in [(816,['Homologous DNA','present?'],'sequence'),
                         (1109,['That region','annotated as exon?'],'exonic'),
                         (1430,['Intron at this','aligned position?'],'intron')]:
        p.local(x-70,187,kind,1,w=140)
        p.lines(x,244,labels,16,step=21,anchor='middle',bold=True)
    p.line(706,277,1604,277,stroke=LINE,width=1)
    for i,s in enumerate('ABCD'):
        y=290+i*62
        if i%2==0:p.rect(705,y-18,900,39,fill=PALE,stroke='none')
        p.text(718,y+7,s,18,bold=True)
        for x,key in [(816,'sequence_presence'),(1109,'exonic_status'),(1430,'intron_presence')]:
            v=CHARACTERS[key][s];label='?' if v=='unknown' else 'N/A' if v=='inapplicable' else str(v)
            p.text(x,y+7,label,23,anchor='middle',fill=AMBER if isinstance(v,str) else INK)
    p.lines(74,542,['1 = present / exonic; 0 = supported absence / intronic.',
        '? = unresolved state; N/A = no corresponding DNA on which to define exonic status.'],16,step=25)
    p.line(48,594,1630,594,stroke=LINE,width=1)
    p.section(48,641,'2','Model a character along a branch')
    p.section(846,641,'3','Infer ancestral states on the species tree')
    p.text(80,690,'Example: intron presence at the highlighted position',17,bold=True)
    p.local(99,742,'intron',0,w=210);p.local(512,742,'intron',1,w=210)
    p.line(328,733,493,733,arrow=True)
    p.line(493,776,328,776,arrow=True)
    p.text(410,717,'gain rate g',17,anchor='middle',fill=PURPLE)
    p.text(410,806,'loss rate l',17,anchor='middle',fill=PURPLE)
    p.text(204,798,'0: uninterrupted exon',16,anchor='middle')
    p.text(617,798,'1: intron present',16,anchor='middle')
    p.text(85,859,'Q = [ -g   g ;  l  -l ]',21,mono=True)
    p.text(85,901,'P(t) = exp(Qt)',24,bold=True)
    p.lines(85,934,['Branch length t and rates determine transition probabilities.',
                    'Within a gene family, each character class shares rates.',
                    'Whole transcripts are not the states of this model.'],16,step=24)
    p.lines(85,1003,['Sum over alternative ancestral states; do not choose a single history',
                    'before calculating the likelihood. Unknown tips allow both states.'],16,step=25,fill=MUTED)
    data=probability_example()
    xy=p.tree(904,731,width=224,step=84,tip_states=CHARACTERS['intron_presence'],posterior=data['node_presence'])
    for s in 'ABCD':p.gene(1189,xy[s][1]-11,width=336,segment=s!='D',exonic=None if s=='C' else True,intron=s!='B')
    for node in ['AB','CD']:
        x,y=xy[node];p.text(x-8,y-25,f"{data['node_presence'][node]:.2f}",16,anchor='end',fill=PURPLE)
    x,y=xy['root'];p.text(x-25,y+5,f"{data['node_presence']['root']:.2f}",16,anchor='end',fill=PURPLE)
    p.pie(919,1041,.75,r=10)
    p.text(940,1047,'Filled fraction: probability of ancestral intron presence.',15)
    p.lines(845,1078,['Illustrative fixed rates, not estimates: g = 0.4, l = 0.2;',
                    'internal branches = 0.4, terminal branches = 1; root P(1) = 2/3. Tree not to scale.'],13,step=19,fill=MUTED)
    return p.write(path)
