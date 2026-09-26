"""Elementary structural edits on retained or source-supported DNA."""
from .drawing import box, line, text, track, document, RULE, PALE, GOLD, TEAL, PURPLE, MUTED


def track_event(x, y, title, before, after, annotation):
    """A wide before/after glyph with one concise interpretation line."""
    b = [box(x,y,530,242,"white",RULE), text(x+20,y+35,title,18,bold=True),
         text(x+18,y+105,"Before",13,color=MUTED), text(x+18,y+168,"After",13,color=MUTED)]
    gx=x+105; width=390
    b += track(gx,y+99,width,100,before,height=21)
    b += [line(gx+width/2,y+118,gx+width/2,y+143,GOLD,2,arrow=True)]
    b += track(gx,y+163,width,100,after,height=21)
    b += [text(x+20,y+218,annotation,13,color=MUTED)]
    return b


def draw():
    b = [text(42, 46, "Elementary evolutionary events", 27, bold=True),
         text(42, 80, "Local configuration edges describe supported structural changes", 15, color=MUTED)]
    b += track_event(42,105,"A   Split",[(8,92)],[(8,46),(54,92)],"One separator; DNA retained")
    b += track_event(625,105,"B   Fusion",[(8,43),(57,92)],[(8,92)],"Intervening DNA retained")
    b += track_event(1208,105,"C   Endpoint shift",[(15,77)],[(25,77)],"One supported endpoint shifts")
    b += track_event(42,370,"D   Exonization",[(7,28),(72,93)],[(7,28),(42,58),(72,93)],"Candidate exon; DNA retained")
    b += track_event(625,370,"E   Exon inactivation",[(7,28),(42,58),(72,93)],[(7,28),(72,93)],"DNA retained")
    # The outlined interval is source-supported DNA entering the locus (m: 0 -> 1).
    b += [box(1208,370,530,242,"white",RULE), text(1228,405,"F   DNA insertion",18,bold=True),
          text(1226,473,"Before",13,color=MUTED), text(1226,536,"After",13,color=MUTED)]
    gx=1305; width=390
    b += track(gx,464,width,100,[(8,42),(58,92)],height=21)
    b += [f'<rect x="{gx+width*.42}" y="{461}" width="{width*.16}" height="6" fill="white"/>',
          line(gx+width/2,483,gx+width/2,508,GOLD,2,arrow=True)]
    b += track(gx,528,width,100,[(8,42),(48,53),(58,92)],height=21)
    b += [f'<rect x="{gx+width*.42}" y="517" width="{width*.16}" height="28" fill="none" stroke="{GOLD}" stroke-width="2"/>',
          text(1228,598,"Material m: absent → present",13,color=MUTED)]
    # One highlighted DNA tract changes two coexisting structures jointly.
    b += [box(42,650,1696,300,PALE,RULE),
          text(64,690,"G   One continuous deletion affects coexisting configurations",19,bold=True),
          text(72,735,"Before",14,color=MUTED), text(838,735,"After",14,color=MUTED)]
    left_x,left_w=135,535
    source_x=left_x+left_w*.36; source_w=left_w*.19
    b += [f'<rect x="{source_x}" y="748" width="{source_w}" height="137" fill="{GOLD}" fill-opacity=".22"/>']
    b += track(left_x,766,left_w,100,[(4,30),(36,50),(80,96)],height=21)
    b += [text(left_x,810,"C₁",13,color=TEAL,bold=True)]
    b += track(left_x,857,left_w,100,[(4,30),(36,67),(80,96)],color=PURPLE,height=21)
    b += [text(left_x,901,"C₂",13,color=PURPLE,bold=True),
          line(702,825,790,825,GOLD,3,arrow=True),text(693,800,"same tract",13,color=GOLD)]
    right_x,right_w=885,500
    b += track(right_x,766,right_w,100,[(4,30),(80,96)],height=21)
    b += [text(right_x,810,"C₁′",13,color=TEAL,bold=True)]
    b += track(right_x,857,right_w,100,[(4,30),(55,67),(80,96)],color=PURPLE,height=21)
    # White masks interrupt the DNA backbone exactly across the deleted source tract.
    gap_x=right_x+right_w*.36; gap_w=right_w*.19
    b += [f'<rect x="{gap_x}" y="763" width="{gap_w}" height="6" fill="{PALE}"/>',
          f'<rect x="{gap_x}" y="854" width="{gap_w}" height="6" fill="{PALE}"/>',
          text(right_x,901,"C₂′",13,color=PURPLE,bold=True),
          text(1450,790,"One DNA event",16,bold=True),
          text(1450,829,"Distinct residual",13,color=MUTED),
          text(1450,863,"configurations",13,color=MUTED)]
    return document(1790, 985, "Elementary evolutionary events", b,
                    "Panels A–F show split, fusion on retained intervening DNA, endpoint shift, exonization, annotation-conditional inactivation, and source-supported DNA insertion. Panel G shows one continuous deletion across two coexisting configurations; their residual structures remain distinct. A deletion can create exon adjacency, with the resulting merge counted as part of that deletion event.")
