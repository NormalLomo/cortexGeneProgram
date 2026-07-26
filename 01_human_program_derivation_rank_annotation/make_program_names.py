#!/usr/bin/env python3
import csv, re, sys, io
import os
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()
IN = str(PROJECT_ROOT / 'results/crossregion_v1/program_annotation_gobp.tsv')
OUT = str(PROJECT_ROOT / 'results/crossregion_v1/program_names.tsv')
FDR_CUT = 0.25
with open(IN, newline='') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
PAREN_GO = re.compile('\\s*\\(GO:\\d+\\)\\s*$')

def strip_go(term):
    return PAREN_GO.sub('', (term or '').strip()).strip()

def strip_pid_prefix(s):
    return re.sub('^P\\d+\\s+', '', (s or '').strip()).strip()
ABBREV = {'\\bReg\\b': 'Regulation', '\\bReg of\\b': 'Regulation of'}

def clean_full(name):
    s = strip_go(name).strip()
    for (pat, rep) in ABBREV.items():
        s = re.sub(pat, rep, s)
    s = re.sub('\\s+', ' ', s).strip()
    return s
SHORT_BY_PROGRAM = {'1': 'Reg. mRNA splicing', '3': 'Cilium movement', '4': 'Fe-S cluster assembly', '5': 'LDL particle removal', '6': 'Pos. reg. secretion', '7': 'Glutamate metab.', '8': 'IT neuropil/cytoskeleton', '9': 'Ca2+ import to cytosol', '10': 'Reg. ECM org.', '11': 'Reg. cardiac contract.', '12': 'Sterol biosynth.', '13': 'Epithelial cell migr.', '14': 'Cell. resp. retinoic ac', '15': 'Eye development', '16': 'Resp. retinoic acid', '17': 'GABAergic synapse', '18': 'Smell chemodetection', '19': 'Retrograde ax. transp.', '20': 'Embryo skeletal dev.', '21': 'Chemosensory detection', '22': 'Neural crest migration', '23': 'Kidney development', '24': 'IT synaptic-vesicle neuropil', '25': 'Neuropeptide signaling', '26': 'Ear morphogenesis', '27': 'Cilium movement (P27)', '28': 'Anion transport', '29': 'Pos. reg. Tyr phosph.', '30': 'Neg. reg. endoth. migr.', '31': 'Negative chemotaxis', '32': 'Reg. muscle contract.', '33': 'Cholinergic synapse', '34': 'Endoderm formation', '35': 'B cell activation', '36': 'Sensory organ dev.', '37': 'Leukocyte adhesion', '38': 'Skeletal system dev.', '39': 'Cholinergic syn. (P39)', '40': 'BCR signaling', '41': 'Taste perception', '42': 'Pos. reg. muscle diff.', '43': 'Skeletal morphogenesis', '44': 'Amino acid transport', '45': 'Myelination', '46': 'Pos. reg. vasoconstr.', '47': 'Immune receptor signal.', '48': 'Sprouting angiogenesis', '49': 'Reg. T cell prolif.', '50': 'Ext. encaps. struct.', '51': 'Pos. reg. miRNA transc.', '52': 'Aromatic AA catab.', '53': 'Neg. reg. growth', '54': 'Pos. reg. lymphocyte', '55': 'Muscle fate commit.', '56': 'Blood vessel morphog.', '57': 'IT axon/neuropil cytoskeleton', '58': 'Antimicrob. peptide', '59': 'Response to zinc ion', '60': 'Antimicrob. humoral'}
SHORT_MAP = [('oxidative phosphoryl', 'Oxidative phosphoryl.'), ('proton transmembrane transport', 'IT synaptic-vesicle neuropil'), ('neuropeptide signaling', 'Neuropeptide signaling'), ('synaptic transmission, gabaergic', 'GABAergic synapse'), ('synaptic transmission, cholinergic', 'Cholinergic synapse'), ('\\bmyelination\\b', 'Myelination')]
STOP_LEAD = {'regulation', 'positive', 'negative', 'cellular', 'detection', 'response'}

def short_label(full):
    low = full.lower()
    for (pat, lab) in SHORT_MAP:
        if re.search(pat, low):
            return lab
    s = full
    s = re.sub('Positive Regulation of', 'Pos. reg.', s)
    s = re.sub('Negative Regulation of', 'Neg. reg.', s)
    s = re.sub('\\bRegulation of\\b', 'Reg.', s)
    s = re.sub('\\bTransmembrane Transport\\b', 'transport', s)
    s = re.sub('\\bTransmembrane\\b', 'transmemb.', s)
    s = re.sub('\\bDevelopment\\b', 'dev.', s)
    s = re.sub('\\bMorphogenesis\\b', 'morphogenesis', s)
    s = re.sub('\\bDifferentiation\\b', 'diff.', s)
    s = re.sub('\\bBiosynthetic Process\\b', 'biosynth.', s)
    s = re.sub('\\bMetabolic Process\\b', 'metab.', s)
    s = re.sub('\\bCatabolic Process\\b', 'catab.', s)
    s = re.sub('\\bOrganization\\b', 'org.', s)
    s = re.sub('\\bSignaling Pathway\\b', 'signaling', s)
    s = re.sub('\\s+', ' ', s).strip()
    if len(s) <= 24:
        return s
    words = s.split()
    out = []
    for w in words:
        cand = ' '.join(out + [w])
        if len(cand) > 24:
            break
        out.append(w)
    if not out:
        out = [words[0][:24]]
    return ' '.join(out).rstrip(',')
out_rows = []
for r in rows:
    prog = r['program'].strip()
    fdr_raw = r['fdr'].strip()
    try:
        fdr = float(fdr_raw)
    except ValueError:
        fdr = float('nan')
    top_term_full = clean_full(r['top_BP_term'])
    proposed = clean_full(strip_pid_prefix(r['proposed_name']))
    sig = fdr < FDR_CUT
    confidence = 'sig' if sig else 'unresolved'
    if proposed and top_term_full.lower().startswith(proposed.lower().rstrip()):
        proposed = top_term_full
    chosen = proposed if sig else top_term_full
    if not chosen:
        chosen = top_term_full or proposed
    name_full = chosen
    name_short = SHORT_BY_PROGRAM.get(prog) or short_label(name_full)
    out_rows.append({'program': prog, 'name_full': name_full, 'name_short': name_short, 'confidence': confidence, 'top_BP_term': top_term_full, 'fdr': fdr_raw, '_proposed': proposed, '_top3': r.get('top3_BP_terms', '')})

def second_term(top3):
    parts = [strip_go(p) for p in (top3 or '').split('|')]
    parts = [p for p in parts if p]
    return parts[1] if len(parts) > 1 else ''
seen = {}
for row in out_rows:
    s = row['name_short']
    if s not in seen:
        seen[s] = [row]
    else:
        seen[s].append(row)
for (s, group) in list(seen.items()):
    if len(group) == 1:
        continue
    for (i, row) in enumerate(group):
        if i == 0:
            continue
        extra = ''
        st = second_term(row['_top3'])
        full_words = [w for w in re.split('[\\s,]+', row['name_full']) if w]
        short_low = row['name_short'].lower()
        disamb = ''
        for w in full_words:
            if w.lower() not in short_low and len(w) > 2:
                disamb = w
                break
        if not disamb and st:
            stw = [w for w in re.split('[\\s,]+', st) if len(w) > 2]
            if stw:
                disamb = stw[0]
        if not disamb:
            disamb = 'P' + row['program']
        cand = (row['name_short'] + ' ' + disamb).strip()
        if len(cand) > 24:
            cand = (row['name_short'][:max(0, 24 - len(disamb) - 1)].rstrip() + ' ' + disamb).strip()
        row['name_short'] = cand
shorts = {}
for row in out_rows:
    if row['name_short'] in shorts:
        row['name_short'] = row['name_short'][:18].rstrip() + ' P' + row['program']
    shorts[row['name_short']] = True
cols = ['program', 'name_full', 'name_short', 'confidence', 'top_BP_term', 'fdr']
buf = io.StringIO()
w = csv.DictWriter(buf, fieldnames=cols, delimiter='\t', lineterminator='\n')
w.writeheader()
for row in out_rows:
    w.writerow({k: row[k] for k in cols})
data = buf.getvalue()
with open(OUT, 'w') as f:
    f.write(data)
n = len(out_rows)
empties = [r['program'] for r in out_rows if not r['name_short'].strip()]
allshort = [r['name_short'] for r in out_rows]
