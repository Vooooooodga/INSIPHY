"""Parsimony, conditional-likelihood pruning and across-gene rate fitting."""
from .drawing import box, line, text, math_text, tree, track, document, PALE, RULE, TEAL, PURPLE, MUTED


def draw():
    b = [text(42, 50, "Ancestral-state reconstruction on a supplied tree", 28, bold=True),
         text(42, 87, "Local configurations or declared repertoires are the modeled states", 15, color=MUTED)]
    # A: local tip observations and the supplied rooted species tree.
    b += [box(42,120,515,525,"white",RULE),
          text(66,163,"A   Tree + tip observations",19,bold=True)]
    b += tree(86,235,.86)
    for i,s in enumerate("ABCD"):
        yy=235+i*60
        b.append(line(86+145*.86+42,yy,300,yy,RULE,1.5))
        b += track(300,yy,215,100,[(4,24),(37,59),(76,96)],height=18)
    b += [text(66, 510, "Tip compatibility", 15, color=MUTED),
          text(66, 548, "Root, topology and branch lengths", 14),
          text(66, 582, "are supplied model conditions", 14)]
    # B: conditional likelihood messages travel from tips toward the root.
    b += [box(590,120,1158,525,PALE,RULE),
          text(614,163,"B   Continuous-time conditional likelihood",19,bold=True)]
    tx,ty,scale=655,251,1.16
    b += tree(tx,ty,scale)
    xa,xb=tx+70*scale,tx+145*scale
    # Branch messages point from each child toward its parent.
    for yy,dy in ((ty+1,-7),(ty+60,8),(ty+120,-7),(ty+180,8)):
        b.append(line(xb-5,yy+dy,xa+8,yy+dy,PURPLE,2.5,arrow=True))
    for yy,dy in ((ty+30,-9),(ty+150,9)):
        b.append(line(xa-4,yy+dy,tx+8,yy+dy,PURPLE,2.5,arrow=True))
    b += [text(932,222,"Tip observations define compatibility weights",14),
          math_text(932,267,[("m",False),("c",True),("(i) = Σ",False),("j",True),(" P",False),("ij",True),("(t",False),("c",True),(") L",False),("c",True),("(j)",False)],15,color=PURPLE,bold=True),
          math_text(932,315,[("L",False),("v",True),("(i) = ∏",False),("c∈children(v)",True),(" m",False),("c",True),("(i)",False)],15,color=PURPLE,bold=True),
          math_text(932,363,[("ℒ = Σ",False),("i",True),(" π",False),("i",True),(" L",False),("root",True),("(i)",False)],15,color=PURPLE,bold=True),
          text(932,421,"P(t) = exp(Qt)",17,color=PURPLE,bold=True),
          text(932,472,"Child messages combine by sum-product",14),
          text(932,508,"Node outputs: marginal ancestral states",14),
          text(932,544,"Branches: conditional event expectations",14)]
    # C–E: complementary analyses and rate treatment.
    b += [box(42,690,530,260,"white",RULE),text(66,733,"C   Maximum parsimony",18,bold=True)]
    for cx,cy,label in ((135,812,"C₁"),(300,812,"C₂"),(465,812,"C₃")):
        b += [f'<circle cx="{cx}" cy="{cy}" r="30" fill="{PALE}" stroke="{TEAL}" stroke-width="2"/>',
              text(cx,cy+7,label,15,anchor="middle",bold=True)]
    b += [line(165,812,270,812,arrow=True),line(330,812,435,812,arrow=True),
          text(66,884,"Minimize edit cost; retain ties",14,color=TEAL)]
    b += [box(630,690,530,260,PALE,RULE),text(654,733,"D   Shared-rate estimation",18,bold=True),
          text(654,786,"≥2 independent genes",15,bold=True),
          text(654,824,"Pooled likelihood; one scale ρ",14),
          text(654,860,"Relative rates · roots · trees fixed",14),
          text(654,903,"Prespecified independent genes",13,color=MUTED)]
    b += [box(1218,690,530,260,"white",RULE),text(1242,733,"E   Assumption sensitivity",18,bold=True),
          text(1242,786,"Recalculate across declared",15),
          text(1242,824,"root-state distributions",14),
          text(1242,860,"branch lengths or catalogues",14),
          text(1242,903,"Conditional structural summaries",13,color=MUTED)]
    b += [math_text(48,990,[("q",False),("ij",True),(" = ρ Σ",False),("e:i→j",True),(" λ",False),("k(e)",True),(" w",False),("e",True),("  (i ≠ j)",False)],15,color=PURPLE,bold=True),
          math_text(48,1027,[("q",False),("ii",True),(" = −Σ",False),("j≠i",True),(" q",False),("ij",True)],15,color=PURPLE,bold=True),
          text(48,1063,"ρ is one shared scale; λ are fixed relative event rates; w allocates outcomes.",14),
          text(48,1098,"Tip constraints depend on annotation coverage; molecular mechanisms are outside this model.",14,color=MUTED)]
    return document(1790, 1125, "Ancestral-state reconstruction on a supplied tree", b,
                    "Panel A shows local tip observations on a supplied rooted species tree. Panel B shows child conditional-likelihood messages flowing toward the root under sum-product pruning, with a root-weighted likelihood sum. Panel C shows minimum-cost parsimony. Panel D shows pooled likelihood fitting of one shared scale across a prespecified collection of at least two independent genes while relative rates, roots and trees are fixed. Panel E shows sensitivity to declared assumptions. Transition rates and branch summaries condition on the declared generator, root conditions and tree. Observation likelihoods use compatibility constraints conditional on supplied annotation coverage. Transcript-detection and ascertainment models are not fitted.")
