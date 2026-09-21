"""Biology-facing V0.19.1 teaching plates; numerical examples come from saved inference."""
from __future__ import annotations
from html import escape
import numpy as np
import math
from .exon_drawing import text, line, rect, exon_track, document, paragraphs


def panel(body, letter, title, x, y, width, height):
    body += [rect(x, y, width, height, "#F7F9FA"), text(x+18, y+31, f"{letter}   {title}", 20)]


def arrow(body, x1, y1, x2, y2, label=""):
    body += [line(x1, y1, x2, y2, "#2B5361")]
    length = max(1., math.hypot(x2-x1,y2-y1))
    dx,dy = (x2-x1)/length,(y2-y1)/length
    body += [line(x2-9*dx+5*dy,y2-9*dy-5*dx,x2,y2),
             line(x2-9*dx-5*dy,y2-9*dy+5*dx,x2,y2)]
    if label: body.append(text((x1+x2)/2, min(y1,y2)-10, label, 12, "middle"))


def tiny_tree(body, x, y, step=60, width=170):
    ys = {s: y+i*step for i,s in enumerate("ABCD")}
    ys.update(AB=(ys["A"]+ys["B"])/2, CD=(ys["C"]+ys["D"])/2)
    ys["root"] = (ys["AB"]+ys["CD"])/2
    xs = {"root": x, "AB": x+width*.48, "CD": x+width*.48, **{s:x+width for s in "ABCD"}}
    for p,c in (("root","AB"),("root","CD"),("AB","A"),("AB","B"),("CD","C"),("CD","D")):
        body += [line(xs[p],ys[p],xs[p],ys[c]),line(xs[p],ys[c],xs[c],ys[c])]
    for s in "ABCD": body.append(text(xs[s]+8,ys[s]+5,s,14))
    return xs,ys


def overview(data):
    body=[text(30,38,"IntraPhy 0.19.1 | Exon structural evolution",27),
          text(30,66,"Genomic sequence + exon annotation + rooted species tree. No transcriptome or exon usage measurements.",15)]
    panel(body,"A","Compare the same local exon structures across species",30,90,1540,330)
    _,ys=tiny_tree(body,70,171,56,190)
    units={u["unit_id"]:u for u in data["units"] if u["family_id"]=="Example_gene"}
    for j,(name,label) in enumerate((("split","One exon / two exons"),("fusion","Two exons / one exon"),
                                   ("deletion","Sequence deletion"),("boundary","Boundary difference"))):
        xx=350+j*285
        body.append(text(xx+110,141,label,14,"middle"))
        obs={o["species"]:o for o in units[name]["catalogue"]["observations"]}
        for s in "ABCD":
            c=obs[s]["configurations"][0]
            body+=exon_track(xx,ys[s],220,180,c["exons"],fill=("#AD783C" if (name=="fusion" and s=="B") or (name=="boundary" and s=="C") or (name=="split" and s=="D") else "#427788"))
            if not c["exons"]: body.append(text(xx+110,ys[s]-12,"DNA absent",12,"middle"))
    body.append(text(65,393,"These are four local comparisons within a schematic gene. They are not four mutually exclusive states of one exon.",13))
    panel(body,"B","Establish correspondence; keep uncertain annotation separate",30,442,1540,240)
    body+=exon_track(75,520,310,180,[(0,180)])
    body+=exon_track(75,600,310,180,[(0,60),(90,180)])
    for a,b in ((0,60),(90,180)):
        l=75+310*a/180;r=75+310*b/180
        body.append(f'<path d="M{l},529 L{r},529 L{r},591 L{l},591 Z" fill="#BCD4DB" opacity="0.65"/>')
    body.append(text(75,646,"Complementary matches, not extra exon births",13))
    arrow(body,425,560,505,560)
    body+=[text(550,514,"Sequence correspondence and native exon boundaries",18),
        text(550,548,"Retain whole-structure alternatives for missing, shifted, split or fused annotations.",15),
        text(550,581,"A candidate supported by DNA is not a confirmed transcript or a fixed evolutionary change.",15),
        text(550,617,"Conflicting copies, orientation or missing genomic context remain unresolved.",15)]
    panel(body,"C","Use the same configurations and elementary edits for history inference",30,704,1540,250)
    body+=exon_track(80,787,250,180,[(0,180)])
    body+=exon_track(80,890,250,180,[(0,60),(90,180)])
    arrow(body,175,812,175,865)
    arrow(body,220,865,220,812)
    body += [text(145,842,"split",13,"end"),text(240,842,"fusion",13)]
    arrow(body,385,837,475,837)
    body += [text(510,780,"Maximum parsimony",19),text(510,813,"Minimum legal edits and alternative branch placements",14),
        text(510,866,"Continuous-time model",19),text(510,899,"Conditional ancestral probabilities and expected edits",14),
        text(1130,786,"One edit can affect several exons.",14),text(1130,822,"Possible histories are not added.",14),
        text(1130,858,"Rates depend on declared opportunities.",14),text(1130,894,"Unresolved does not mean conserved.",14)]
    body.append(text(30,986,"Scope: collinear, locally comparable exon structures. Copy genealogy, shuffling and inversions are not reconstructed.",14))
    return document(1600,1010,body,"IntraPhy exon structural evolution — overview")


def homology(data):
    body=[text(30,38,"How exon correspondence and annotation uncertainty enter the analysis",25),
          text(30,67,"Current implementation: one genomic alignment with coding cross-checks, not joint integration over all alignments.",14)]
    panel(body,"A","Keep the full available locus, then check homologous coordinates",30,91,1540,275)
    _,ys=tiny_tree(body,65,158,48,170)
    for s in "ABCD":
        body+=exon_track(325,ys[s],440,180,[(0,60),(90,180)] if s=="D" else [(0,180)])
        body += [line(290,ys[s],320,ys[s],"#9AA8B0"),line(770,ys[s],810,ys[s],"#9AA8B0")]
    body += [text(860,154,"FASTA and GFF are kept in a reversible coordinate system.",16),
        text(860,194,"MAFFT: genomic columns; coding projection: conflict check.",16),
        text(860,234,"minimap2: competing copies and orientation candidates.",16),
        text(860,274,"Unannotated terminal DNA is not cropped from the search.",16),
        text(860,323,"CESAR2 remains a separate optional prediction adapter.",13)]
    panel(body,"B","Distinguish complementary coverage from competing copies",30,391,750,247)
    body+=exon_track(75,469,245,180,[(0,180)])
    body+=exon_track(75,568,245,180,[(0,60),(90,180)])
    for a,b in ((0,60),(90,180)):
        l=75+245*a/180;r=75+245*b/180
        body.append(f'<path d="M{l},478 L{r},478 L{r},560 L{l},560 Z" fill="#ADCBD4" opacity="0.6"/>')
    body+= [text(75,610,"One-to-many complementary correspondence",12)]
    body+=exon_track(425,469,220,180,[(0,180)])
    body+=exon_track(380,568,320,360,[(0,160),(200,360)])
    body += [line(450,480,420,551,"#AE6746"),line(620,480,650,551,"#AE6746"),
             text(410,521,"Competing matches of the same interval",12),text(390,610,"Retain ambiguity; no copy history inferred",12)]
    panel(body,"C","Compare complete alternative annotations",805,391,765,247)
    body+=exon_track(865,475,430,180,[(0,60),(90,180)])
    body+= [text(1315,480,"GFF structure",13)]
    body+=exon_track(865,544,430,180,[(0,180)],"none","5,3")
    body+= [text(1315,549,"candidate",13),text(843,589,"Require matched exon sequence, location flanks and preserved cut context.",14),
             text(843,619,"The original GFF is not rewritten; candidates carry a source and evidence.",14)]
    panel(body,"D","Report how the conclusion depends on annotation",30,663,1540,244)
    control=next(x for x in data["summary"] if x["family_id"]=="Annotation_control")
    a=control["annotation_conditional_minima"][0];b=control["evidence_compatible_minima"][0]
    body += [text(77,736,"Same synthetic DNA; one annotation is split",19),
        text(77,781,f"Exact GFF condition: minimum edits = {a}",18),
        text(77,823,f"Allow the sequence-supported whole-exon candidate: minimum edits = {b}",18),
        text(930,746,"This is annotation sensitivity, not proof that",17),text(930,779,"the more parsimonious annotation is correct.",17),
        text(930,827,"Assembly gaps and missing search context remain unknown.",14)]
    body.append(text(30,941,"Thresholds qualify candidates; they are not calibrated probabilities of homology or biological correctness.",14))
    return document(1600,970,body,"Exon homology and atomic annotation alternatives")


def probability_model(data):
    unit=next(u for u in data["units"] if u["unit_id"]=="split")
    c=unit["catalogue"];posterior=unit["ctmc"]
    body=[text(30,39,"From exon configurations to a phylogenetic probability model",26),
        text(30,69,"One local object, several possible structures. Example rates are fixed teaching values, not biological estimates.",14)]
    panel(body,"A","Define structures and the edits that connect them",30,92,640,334)
    body += exon_track(80,186,250,180,[(0,180)])
    body += exon_track(80,321,250,180,[(0,60),(90,180)])
    body += [text(385,190,"One exon",18),text(385,325,"Two exons",18)]
    arrow(body,170,211,170,286)
    arrow(body,210,286,210,211)
    body+= [text(145,242,"split",15,"end"),text(145,263,f'rate = {data["fixed_teaching_rates"]["split"]}',13,"end"),
            text(240,242,"fusion",15),text(240,263,f'rate = {data["fixed_teaching_rates"]["fusion"]}',13),
            text(80,377,"Reverse fusion has its own rate. Other legal edits remain in the model.",13),
            text(80,403,f"Full declared catalogue: {len(unit['states'])} states; not only observed endpoints.",13)]
    panel(body,"B","Build a transition process from these same edits",698,92,872,334)
    body += [text(735,169,"Each arrow connects two concrete exon configurations.",19),
        text(735,217,"q(S, S′) = sum of rates for allowed edits from S to S′",20),
        text(735,265,"P(t) = exp(Q t)",25),
        text(735,303,"Branch length controls the time available for zero, one or multiple edits.",15),
        text(735,343,"An insertion or deletion includes its structural consequences only once.",15),
        text(735,382,"This is a finite opportunity model, not a DNA-repair mechanism model.",14)]
    panel(body,"C","Combine observed tip structures with the species tree",30,450,1540,379)
    xs,ys=tiny_tree(body,100,525,67,440)
    obs={o["species"]:o for o in c["observations"]}
    for s in "ABCD": body+=exon_track(650,ys[s],310,180,obs[s]["configurations"][0]["exons"])
    for node in ("root","AB","CD"):
        values=np.asarray(posterior["nodes"][node]);index=int(np.argmax(values));p=float(values[index]);state=unit["states"][index]
        body += [f'<g data-node="{node}" data-probability="{p:.17g}">']
        body+=exon_track(xs[node]-30,ys[node]-40,120,180,state["exons"],height=11)
        body+=[text(xs[node]-30,ys[node]-59,f'Most probable: {state["state_id"]}',11),
               text(xs[node]-30,ys[node]-20,f"P = {p:.3f}",13),'</g>']
    body += [text(1030,525,"Ancestral configurations are unobserved.",17),
        text(1030,562,"Pruning sums over all legal ancestral states.",17),
        text(1030,608,"Node labels: conditional state probabilities.",15),
        text(1030,647,"An uncertain tip allows compatible states;",15),text(1030,675,"it is not forced to a missing exon.",15),
        text(1030,733,"The most probable marginal states need not",15),text(1030,759,"form one most probable joint history.",15)]
    panel(body,"D","Three different quantities for the branch CD → D",30,851,1540,177)
    b=next(b for b in posterior["branches"] if b["child"]=="D")
    body += [text(80,925,f'P(different endpoints) = {b["probability_different_endpoints"]:.4f}',18),
        text(615,925,f'P(at least one edit) = {b["probability_at_least_one_edit"]:.4f}',18),
        text(1150,925,f'Expected edits = {b["expected_edits"]:.4f}',18),
        text(80,990,"All values are calculated by this version of IntraPhy; full rates, states, observations and outputs accompany the figures.",14)]
    return document(1600,1060,body,"Exon configurations and their conditional CTMC history")


def architecture(data):
    body=[text(30,40,"IntraPhy 0.19.1 | Data and algorithm responsibilities",26)]
    boxes=[("Genomic FASTA + GFF + species tree","inputs/; preparation/",80,120),
           ("Whole-locus alignment and native exon coordinates","structure/alignment.py; corroboration.py",80,255),
           ("Atomic annotation alternatives + physical-unit validation","structure/alternatives.py; validation.py",80,390),
           ("Finite exon configurations + one elementary edit registry","structure/space.py; edits.py",80,525)]
    for title,path,x,y in boxes:
        body += [rect(x,y,1370,98,"#F7F9FA"),text(x+25,y+36,title,21),text(x+25,y+72,path,15)]
    for y in (218,353,488): arrow(body,760,y+5,760,y+32)
    body += [rect(80,696,650,130,"#EDF4F5"),rect(800,696,650,130,"#EDF4F5"),
        text(115,734,"Generalized Sankoff",23),text(115,770,"Sparse shortest paths; all-optimal histories",16),
        text(835,734,"Finite-state CTMC",23),text(835,770,"Canonical tree; bounded kernels; conditional probabilities",15)]
    arrow(body,760,635,405,675);arrow(body,760,635,1110,675)
    body += [rect(80,891,1370,118,"#F7F9FA"),text(110,928,"Saved states, edit IDs, uncertainty and probabilities → read-only result figures",21),
        text(110,969,"Reporting does not repair annotations, select a preferred biological history or refit rates.",16)]
    arrow(body,405,837,405,880);arrow(body,1110,837,1110,880)
    return document(1530,1040,body,"IntraPhy 0.19.1 software architecture")


PLATES={"method_overview": (overview,"Scope, inputs and structural inference"),
        "homology_inference": (homology,"Whole exons, paired sequence and atomic annotation alternatives"),
        "phylogenetic_model": (probability_model,"Configurations, elementary edits and calculated ancestral probabilities"),
        "architecture": (architecture,"Implemented data flow and code ownership")}
