"""A two-row biological overview from homologous loci to phylogenetic summaries."""
from .drawing import box, line, text, tree, track, document, RULE, TEAL, MUTED, GOLD, PURPLE, PALE


def panel(x, y, title, number):
    return [box(x, y, 540, 455, "white", RULE),
            text(x+24, y+43, number, 18, color=TEAL, bold=True),
            text(x+65, y+43, title, 19, bold=True)]


def draw():
    b = [text(42, 54, "From genomic loci to local structural evolution", 28, bold=True),
         text(42, 88, "Evidence, states and evolutionary summaries at the same homologous interval", 15, color=MUTED)]
    cells = ((42,120),(625,120),(1208,120),(42,625),(625,625),(1208,625))
    titles = ("Orthologous loci", "Sequence correspondence", "Local configurations",
              "Admissible state graph", "Tree-based inference", "Ancestral reconstruction")
    for i, ((x,y), title) in enumerate(zip(cells,titles),1):
        b += panel(x,y,title,"ABCDEF"[i-1])
    # 1: Input loci and the supplied species tree.
    b += tree(80, 233, .75, labels=False)
    for i,s in enumerate("ABCD"):
        yy=233+i*60
        b.append(line(189,yy,241,yy,RULE,1.5))
        b += [text(247,yy+7,s,14,bold=True)]
        b += track(276,yy,240,100,[(4,22),(35,48),(66,80),(88,98)],height=17)
    b += [text(68, 511, "Genome + annotation", 14, color=MUTED),
          text(68, 544, "Rooted species tree supplied", 14, color=MUTED)]
    # 2: Correspondence evidence on genomic and coding tracks.
    x,y=cells[1]
    b += [text(x+30,y+105,"Genomic locus",14,color=MUTED)]
    b += track(x+30,y+133,420,100,[(3,24),(29,50),(67,83),(88,98)],height=24)
    b += [line(x+150,y+157,x+150,y+220,GOLD,2),line(x+294,y+157,x+294,y+220,GOLD,2),
          line(x+399,y+157,x+399,y+220,GOLD,2),text(x+30,y+239,"Coding alignment → genome",14,color=MUTED)]
    b += track(x+30,y+267,420,100,[(3,24),(29,50),(67,83),(88,98)],color=PURPLE,height=16)
    b += [text(x+30,y+409,"Copy · order · strand · flanks",14,color=MUTED)]
    # 3: Local ordered structures.
    x,y=cells[2]
    for i,(s,exons) in enumerate((("A",[(3,26),(39,61),(74,97)]),
                                  ("B",[(3,26),(39,97)]),
                                  ("C",[(3,26),(56,61),(74,97)]),
                                  ("D",[(3,26),(74,97)]))):
        yy=y+153+i*61
        b += [text(x+30,yy+6,s,14,bold=True)] + track(x+70,yy,405,100,exons,height=19)
    b += [text(x+30,y+400,"Observed · partial · unknown",14,color=MUTED),
          text(x+30,y+434,"Coexisting configurations",14,color=MUTED)]
    # 4: Admissible states and structural edit edges.
    x,y=cells[3]
    for cx,cy,label in ((x+145,y+190,"C₁"),(x+365,y+190,"C₂"),(x+255,y+330,"C₃")):
        b += [f'<circle cx="{cx}" cy="{cy}" r="47" fill="{PALE}" stroke="{TEAL}" stroke-width="3"/>',
              text(cx,cy+8,label,18,anchor="middle",bold=True)]
    b += [line(x+192,y+190,x+318,y+190,arrow=True),line(x+354,y+230,x+288,y+288,arrow=True),
          text(x+35,y+405,"Observed + admissible states",14,color=MUTED)]
    # 5: Parsimony and CTMC on the supplied tree.
    x,y=cells[4]
    b += tree(x+63,y+135,.78)
    b += [text(x+270,y+164,"Parsimony",17,bold=True),text(x+270,y+198,"min-cost paths",14,color=TEAL),
          text(x+270,y+263,"CTMC",17,bold=True),text(x+270,y+297,"sum over states",14,color=PURPLE),
          text(x+35,y+404,"Species tree supplied",14,color=MUTED),
          text(x+35,y+438,"Local structure reconstructed",14,color=MUTED)]
    # 6: Outputs are symbolic summaries rather than plotted results.
    x,y=cells[5]
    b += tree(x+60,y+150,.78)
    b += [text(x+157,y+242,"pᵢ",17,anchor="middle",color=PURPLE,bold=True),
          text(x+35,y+397,"Marginal state probabilities",14,bold=True),
          text(x+35,y+431,"Expected events along branches",13,color=MUTED)]
    # Arrows show reading order within rows; panel letters link the rows.
    for ax,ay,bx,by in ((582,350,612,350),(1165,350,1195,350),
                        (582,850,612,850),
                        (1165,850,1195,850)):
        b.append(line(ax,ay,bx,by,TEAL,2,arrow=True))
    return document(1790, 1110, "From genomic loci to local structural evolution", b,
                    "Panels A–F show orthologous loci with a supplied rooted species tree, sequence correspondence, ordered local configurations, admissible edit transitions, parsimony and continuous-time inference, and symbolic marginal ancestral-state and event summaries.")
