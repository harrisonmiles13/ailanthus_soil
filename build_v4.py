# Build v4 document.xml: figure surgery + red-flagged substantive edits + silent typo fixes.
import copy, shutil
from lxml import etree

NS = {
 'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
 'w14':'http://schemas.microsoft.com/office/word/2010/wordml',
 'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
 'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
 'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
 'pic':'http://schemas.openxmlformats.org/drawingml/2006/picture',
 'xml':'http://www.w3.org/XML/1998/namespace',
}
def w(t):  return '{%s}%s'%(NS['w'],t)
def q(p,t):return '{%s}%s'%(NS[p],t)

DOC='v4_work/word/document.xml'
tree=etree.parse(DOC); root=tree.getroot(); body=root.find(w('body'))

def ptext(p): return ''.join(t.text or '' for t in p.iter(w('t')))
def find(sub):
    for p in body.findall(w('p')):
        if sub in ptext(p): return p
    return None

_pid=[0x5E000000]
def fresh_ids(p):
    _pid[0]+=1
    p.set(q('w14','paraId'), '%08X'%_pid[0])
    p.set(q('w14','textId'), '77777777')

def make_run(text, bold=False, italic=False, red=False, sub=False, sup=False):
    r=etree.Element(w('r')); rpr=etree.SubElement(r,w('rPr'))
    rf=etree.SubElement(rpr,w('rFonts'))
    for a in ('ascii','eastAsia','hAnsi','cs'): rf.set(w(a),'Times New Roman')
    if bold: etree.SubElement(rpr,w('b')); etree.SubElement(rpr,w('bCs'))
    if italic: etree.SubElement(rpr,w('i')); etree.SubElement(rpr,w('iCs'))
    if red: etree.SubElement(rpr,w('color')).set(w('val'),'FF0000')
    if sub: etree.SubElement(rpr,w('vertAlign')).set(w('val'),'subscript')
    if sup: etree.SubElement(rpr,w('vertAlign')).set(w('val'),'superscript')
    t=etree.SubElement(r,w('t')); t.set(q('xml','space'),'preserve'); t.text=text
    return r

def runs_red(spec):
    # spec: list of (text, italic?) -> all red
    out=[]
    for item in spec:
        if isinstance(item,tuple): out.append(make_run(item[0], red=True, italic=item[1]))
        else: out.append(make_run(item, red=True))
    return out

def para_like(model, runs):
    newp=copy.deepcopy(model)
    for ch in list(newp):
        if ch.tag in (w('r'), w('hyperlink')): newp.remove(ch)
    fresh_ids(newp)
    for rn in runs: newp.append(rn)
    return newp

def insert_after(anchor, new): anchor.addnext(new)

def replace_plain(p, old, new):
    for r in p.findall(w('r')):
        t=r.find(w('t'))
        if t is not None and t.text and old in t.text:
            t.text=t.text.replace(old,new); return True
    return False

def replace_red(p, old, new):
    for r in p.findall(w('r')):
        t=r.find(w('t'))
        if t is not None and t.text and old in t.text:
            pre,post=t.text.split(old,1)
            t.text=pre; t.set(q('xml','space'),'preserve')
            idx=p.index(r)
            newrun=make_run(new, red=True)
            postrun=copy.deepcopy(r); pt=postrun.find(w('t')); pt.text=post; pt.set(q('xml','space'),'preserve')
            p.insert(idx+1,newrun); p.insert(idx+2,postrun)
            return True
    return False

log=[]
def L(ok,msg): log.append(('OK ' if ok else '!! ')+msg)

# ===== 1. FIGURE SURGERY: single Figure 5 (12-panel) =========================
bes_img=None; remove=[]
for p in body.findall(w('p')):
    embeds=[b.get(q('r','embed')) for b in p.iter(q('a','blip'))]
    txt=ptext(p)
    if any(e in ('rId24','rId25','rId26') for e in embeds): remove.append(p); continue
    if txt.startswith('Figure 2. As in') or txt.startswith('Figure 3. As in') or txt.startswith('Figure 4. As in'):
        remove.append(p); continue
    if 'rId23' in embeds: bes_img=p
# resize bes image -> combined aspect (1900x2200 px); width 6.5in=5943600 EMU
for ext in bes_img.iter(q('wp','extent')): ext.set('cx','5943600'); ext.set('cy','6882063')
for ext in bes_img.iter(q('a','ext')):     ext.set('cx','5943600'); ext.set('cy','6882063')
for d in bes_img.iter(q('wp','docPr')):    d.set('name','Figure 5')
for c in bes_img.iter(q('pic','cNvPr')):   c.set('name','combined_panel.png')
for p in remove: body.remove(p)
L(bes_img is not None, 'bes image found/resized; removed %d response paras'%len(remove))
shutil.copy('output/figures/combined_panel.png','v4_work/word/media/image8.png')
L(True,'media/image8.png replaced with combined_panel.png')

# rewrite bes caption -> Figure 5 (red)
cap=find('Within-plot relationships between the time-matched disease influence')
cap_runs=[make_run('Figure 5. ', bold=True, red=True),
  make_run('Within-plot relationships between the time-matched disease-influence index (ANI', red=True),
  make_run('disease', red=True, sub=True),
  make_run('; horizontal axis) and three bioassay responses—total seedling length, root-to-shoot ratio, and germination (as a proportion)—for the four bioassay species, with species in rows and responses in columns. Rows are (top to bottom) ', red=True),
  make_run('Rudbeckia hirta', red=True, italic=True), make_run(' (black-eyed Susan), ', red=True),
  make_run('Sorghum halepense', red=True, italic=True), make_run(' (johnsongrass), ', red=True),
  make_run('Ageratina altissima', red=True, italic=True), make_run(' (white snakeroot), and ', red=True),
  make_run('Achillea millefolium', red=True, italic=True),
  make_run(' (yarrow). Each point is one soil tray, with symbol shape denoting sampling month (May, July, or September) and color distinguishing the three plots, which enter the model only as a fixed nuisance block. Fitted lines are the within-plot relationships from a separate model of each response on the disease-influence index with plot held constant (plot partialled out): straight lines for total seedling length and root-to-shoot ratio (Gaussian linear models) and a logistic curve for germination (a binomial seed-count model with quasibinomial, overdispersion-corrected inference). Each panel is annotated with its within-plot slope per standard deviation of the index and the corresponding P value; the length and root-to-shoot-ratio slopes are in response units and the germination slope is on the log-odds scale. Because inoculation treatment is confounded with plot by design, the colors are not compared and no difference among plots or treatments is implied.', red=True)]
for ch in list(cap):
    if ch.tag==w('r'): cap.remove(ch)
for rn in cap_runs: cap.append(rn)
L(cap is not None,'bes caption rewritten -> Figure 5')

# ===== 2. IN-TEXT FIGURE REFS (red) ==========================================
bio=find('the lowest P among the twelve associations')
L(replace_red(bio,'(Figure 4)','(Figure 5)'),'in-text ref yarrow Fig4->Fig5')
L(replace_red(bio,'; Figure 1)','; Figure 5)'),'in-text ref bes Fig1->Fig5')
L(replace_red(bio,'Figures 2 and 3)','Figure 5)'),'in-text ref jg/wsr Figs2,3->Fig5')

# ===== 3. METHODS / RESULTS STAT FIXES (red) =================================
sens=find('more likely to be artifacts of model parameterization')
L(replace_red(sens,'primary kernel (α = 0,1,2; β = 0.15','primary kernel (α = 2, β = 0.15'),'fix primary-kernel alpha')
L(replace_red(sens,
  'the trends observed in black-eyed Susan and yarrow more likely to be artifacts of model parameterization than a replicable observation of their biology.',
  'the black-eyed Susan trend was insensitive to kernel choice whereas the yarrow trend was not, and given their marginal significance and the multiple-comparison context, neither is interpreted as an established biological effect.'),
  'fix artifacts sentence (bes robust, yarrow sensitive)')
dm=find('weeks after inoculation')
L(replace_red(dm,'and again 4, 8, and 12 weeks after inoculation.',
  'and again at four approximately monthly intervals after inoculation, the last coinciding with the September soil collection.'),
  'fix disease-monitoring reading count')

# ===== 4. MISSING-CITATION FLAGS (red entries in Literature Cited) ===========
ehr=find('Ehrenfeld JG. 2003.')
fradin=para_like(ehr,[make_run('Fradin EF, Thomma BPHJ. 2006. Physiology and molecular aspects of Verticillium wilt diseases caused by V. dahliae and V. albo-atrum. Mol Plant Pathol 7:71–86. ', red=True), make_run('[AI-supplied entry — verify before submission.]', red=True, italic=True)])
insert_after(ehr,fradin)
franz=para_like(ehr,[make_run('[Reference entry needed: “Franzluebbers et al. 2025” is cited in Soil characteristics but is absent from Literature Cited — please add the full citation.]', red=True)])
insert_after(fradin,franz)
medv=find('Medina-Villar S,')
miles=para_like(ehr,[make_run('[Reference entry needed: “Miles et al. 2024” is cited in the Introduction but is absent from Literature Cited — please add the full citation.]', red=True)])
insert_after(medv,miles)
schallb=find('Schall MJ, Davis DD. 2009b.')
shiv=para_like(ehr,[make_run('[Reference entry needed: “Shively et al. 2025” is cited in the Introduction but is absent from Literature Cited — please add the full citation.]', red=True)])
insert_after(schallb,shiv)
L(all(x is not None for x in (ehr,medv,schallb)),'inserted 4 red citation flags (Fradin, Franzluebbers, Miles, Shively)')

# ===== 5. EXPANDED DISCUSSION (red, after the soil-legacy paragraph) =========
leg=find('otherwise flat across the disease influence gradient')
exp1=para_like(leg,[
  make_run('Two considerations frame this single-season result. First, the pathways by which wilt could reshape these soils act gradually. As inoculated trees decline, the living inputs that distinguish ', red=True),
  make_run('Ailanthus', red=True, italic=True),
  make_run(' soils—root exudates and the allelochemical ailanthone delivered through root bark, litter, and fine-root turnover—diminish only slowly, while the countervailing pulse of decomposing root and shoot tissue accrues just as gradually as trees die and their tissues break down. Within a single growing season the inoculated trees were still standing and only beginning to express disease, so neither pathway would yet be expected to have substantially altered bulk soil. The flat bioassay response across the disease-influence gradient is therefore consistent with the timescale over which such a soil legacy would form or erode, rather than evidence that wilt leaves the soil unchanged.', red=True)])
exp2=para_like(leg,[
  make_run('Second, the absence of a soil-mediated signal bears on the safety case for ', red=True),
  make_run('V. nonalfalfae', red=True, italic=True),
  make_run(' as an augmentative biological control. Host-range testing and risk assessments establish that the fungus poses limited direct pathogenic risk to non-target plants (Kasson et al. 2015, O’Neal and Davis 2015); the present results complement that work by addressing an indirect, soil-mediated pathway, indicating that in the first season after inoculation the treatment did not detectably change the soil’s capacity to support germination and early growth of co-occurring species. For an agent intended for operational release, a neutral near-term soil legacy is a desirable property.', red=True)])
insert_after(leg,exp2); insert_after(leg,exp1)
L(leg is not None,'inserted 2 expanded-discussion paragraphs (red)')

# ===== 6. SILENT TYPO FIXES (black) =========================================
L(replace_plain(find('Department of Forest Resources and Environmental Cosnervation'),'Cosnervation','Conservation'),'typo Cosnervation->Conservation')
L(replace_plain(find('suggestive andprovisional'),'andprovisional','and provisional'),'typo andprovisional')
L(replace_plain(find('Opus 4.8) was used'),'Opus 4.8) was used','Opus 4.8) were used'),'typo was->were used')

tree.write(DOC, xml_declaration=True, encoding='UTF-8', standalone=True)
print('\n'.join(log))
print('\nWROTE', DOC)
