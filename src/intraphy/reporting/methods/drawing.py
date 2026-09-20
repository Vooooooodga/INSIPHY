"""SVG vocabulary: genomic intervals, exon boundaries, correspondences and trees."""
from __future__ import annotations
from pathlib import Path
from xml.etree.ElementTree import Element, SubElement, ElementTree

INK = '#202D38'
MUTED = '#576875'
LINE = '#B4C0C7'
TEAL = '#168577'
PURPLE = '#7554A3'
AMBER = '#B97716'
GREY = '#E3E8EB'
PALE = '#F4F7F8'

class Canvas:
    def __init__(self, title, subtitle, width=1680, height=1020):
        self.width, self.height = width, height
        self.root = Element('svg', xmlns='http://www.w3.org/2000/svg',
            width=str(width), height=str(height), viewBox=f'0 0 {width} {height}',
            **{'data-purpose':'synthetic-teaching-schematic'})
        SubElement(self.root, 'title').text = title
        SubElement(self.root, 'desc').text = subtitle + ' Original synthetic example; not a biological validation.'
        d=SubElement(self.root,'defs')
        mark=SubElement(d,'marker',id='arrow',viewBox='0 0 10 10',refX='9',refY='5',
                        markerWidth='6',markerHeight='6',orient='auto-start-reverse')
        SubElement(mark,'path',d='M0 0L10 5L0 10Z',fill=INK)
        pat=SubElement(d,'pattern',id='missing',width='7',height='7',patternUnits='userSpaceOnUse')
        SubElement(pat,'rect',width='7',height='7',fill='#FFF8EA')
        SubElement(pat,'path',d='M0 7L7 0',stroke=AMBER,**{'stroke-width':'.7'})
        self.rect(0,0,width,height,fill='white',stroke='none')
        self.text(48,49,title,30,bold=True)
        self.text(48,83,subtitle,18,fill=MUTED)
        self.text(48,height-22,'SYNTHETIC EXAMPLE  |  Genomic distances are schematic  |  IntraPhy 0.18',13,fill=MUTED)
    def text(self,x,y,s,size=18,bold=False,fill=INK,anchor='start',mono=False):
        e=SubElement(self.root,'text',x=str(x),y=str(y),fill=fill,
            **{'font-family':'DejaVu Sans Mono, monospace' if mono else 'DejaVu Sans, Arial, sans-serif',
               'font-size':str(size),'font-weight':'600' if bold else '400','text-anchor':anchor})
        e.text=str(s);return e
    def lines(self,x,y,items,size=17,step=26,**kw):
        for i,s in enumerate(items):self.text(x,y+i*step,s,size,**kw)
    def rect(self,x,y,w,h,fill='white',stroke=INK,dash=False,r=0,**kw):
        return SubElement(self.root,'rect',x=str(x),y=str(y),width=str(w),height=str(h),
            fill=fill,stroke=stroke,rx=str(r),**{'stroke-width':'1.4',
            'stroke-dasharray':'6 4' if dash else 'none',**{k:str(v) for k,v in kw.items()}})
    def line(self,x,y,u,v,stroke=INK,width=2,dash=False,arrow=False):
        a={'stroke-width':str(width),'stroke-dasharray':'5 4' if dash else 'none'}
        if arrow:a['marker-end']='url(#arrow)'
        return SubElement(self.root,'line',x1=str(x),y1=str(y),x2=str(u),y2=str(v),stroke=stroke,**a)
    def path(self,d,fill='none',stroke=INK,**kw):
        return SubElement(self.root,'path',d=d,fill=fill,stroke=stroke,
                          **{k:str(v) for k,v in kw.items()})
    def section(self,x,y,number,title,subtitle=None):
        self.text(x,y,number,23,bold=True,fill=TEAL)
        self.text(x+34,y,title,23,bold=True)
        if subtitle:self.text(x+34,y+29,subtitle,16,fill=MUTED)
    def ribbon(self,a,b,y,c,d,z,color=TEAL):
        mid=(y+z)/2
        self.path(f'M{a} {y}C{a} {mid} {c} {mid} {c} {z}L{d} {z}C{d} {mid} {b} {mid} {b} {y}Z',
            fill=color,stroke='none',**{'fill-opacity':'.20','data-evidence':'matched-subinterval-only'})
    def dot(self,x,y,fill='white',r=6):
        return SubElement(self.root,'circle',cx=str(x),cy=str(y),r=str(r),fill=fill,stroke=INK,
                          **{'stroke-width':'1.5'})
    def pie(self,x,y,p,r=17):
        from math import sin,cos,pi
        self.dot(x,y,'white',r)
        if p>=1:self.dot(x,y,PURPLE,r);return
        if p<=0:return
        end=(x+r*sin(2*pi*p),y-r*cos(2*pi*p))
        self.path(f'M{x} {y}L{x} {y-r}A{r} {r} 0 {int(p>.5)} 1 {end[0]} {end[1]}Z',
                  fill=PURPLE,stroke='none')
        SubElement(self.root,'circle',cx=str(x),cy=str(y),r=str(r),fill='none',stroke=INK)
    def gene(self,x,y,width=350,segment=True,exonic=True,intron=True,highlight=True):
        """Three exon regions: the middle DNA target and a separate intron position."""
        s=width/350;h=22
        pos=lambda a:x+a*s
        self.line(x,y+h/2,pos(350),y+h/2,width=1.7)
        self.rect(pos(8),y,48*s,h,fill=GREY)
        if segment:
            if exonic is True:self.rect(pos(104),y,64*s,h,fill=TEAL if highlight else GREY)
            elif exonic is False:self.rect(pos(104),y+8,64*s,6,fill=TEAL,stroke='none')
            else:self.rect(pos(104),y,64*s,h,fill='url(#missing)',stroke=AMBER,dash=True)
        else:
            self.line(pos(112),y+7,pos(125),y+16,stroke=TEAL)
            self.line(pos(127),y+7,pos(140),y+16,stroke=TEAL)
            self.line(pos(142),y+7,pos(155),y+16,stroke=TEAL)
        color=PURPLE if highlight else GREY
        if intron is True:
            self.rect(pos(226),y,43*s,h,fill=color)
            self.rect(pos(294),y,43*s,h,fill=color)
        elif intron is False:self.rect(pos(226),y,111*s,h,fill=color)
        else:self.rect(pos(226),y,111*s,h,fill='url(#missing)',stroke=AMBER,dash=True)
        return {'segment':(pos(104),pos(168)),'left':(pos(226),pos(269)),
                'right':(pos(294 if intron else 269),pos(337))}
    def local(self,x,y,kind,state,w=140,color=None):
        color=color or (PURPLE if kind=='intron' else TEAL)
        self.line(x,y+11,x+w,y+11,width=1.5)
        if kind=='intron':
            if state==1:
                self.rect(x+6,y,44,22,fill=color);self.rect(x+w-50,y,44,22,fill=color)
            elif state==0:self.rect(x+6,y,w-12,22,fill=color)
            else:self.rect(x+6,y,w-12,22,fill='url(#missing)',dash=True,stroke=AMBER)
        elif kind=='exonic':
            if state==1:self.rect(x+25,y,w-50,22,fill=color)
            elif state==0:self.rect(x+25,y+8,w-50,6,fill=color,stroke='none')
            else:self.rect(x+25,y,w-50,22,fill='url(#missing)',dash=True,stroke=AMBER)
        elif state==1:self.rect(x+25,y+5,w-50,12,fill=color,stroke='none')
        elif state==0:
            self.line(x+50,y+4,x+62,y+18,stroke=color);self.line(x+65,y+4,x+77,y+18,stroke=color)
        else:self.rect(x+25,y,w-50,22,fill='url(#missing)',dash=True,stroke=AMBER)
    def tree(self,x,y,width=145,step=76,labels=True,tip_states=None,posterior=None):
        ys={s:y+i*step for i,s in enumerate('ABCD')}
        xy={'root':(x,y+1.5*step),'AB':(x+width*.43,y+.5*step),
            'CD':(x+width*.43,y+2.5*step),**{s:(x+width,z) for s,z in ys.items()}}
        for par,kids in [('root',('AB','CD')),('AB',('A','B')),('CD',('C','D'))]:
            px,_=xy[par]
            self.line(px,xy[kids[0]][1],px,xy[kids[1]][1])
            for c in kids:self.line(px,xy[c][1],xy[c][0],xy[c][1])
        if posterior:
            for n,p in posterior.items():self.pie(*xy[n],p)
        if tip_states:
            for n,v in tip_states.items():self.dot(*xy[n],PURPLE if v==1 else 'white',r=7)
        if labels:
            for s,z in ys.items():self.text(x+width+14,z+6,s,17,bold=True)
        return xy
    def write(self,path):
        p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
        ElementTree(self.root).write(p,encoding='utf-8',xml_declaration=True);return p
