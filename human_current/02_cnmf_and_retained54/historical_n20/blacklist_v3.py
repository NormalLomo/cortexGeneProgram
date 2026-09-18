"""V3 gene blacklist: V2 + nuclear-mito + ambient-RNA panel."""
import re

# V2 regex patterns (unchanged)
PATTERNS = [
    r'^RPS\d+', r'^RPL\d+', r'^RPLP\d+',
    r'^MRPS\d+', r'^MRPL\d+',
    r'^MT-',
    r'^EEF1[ABDG]', r'^EEF2$',
    r'^EIF[0-9]', r'^MTRNR[0-9]',
]
SUFFIX = [r'-AS\d*$', r'-DT$', r'-IT\d*$']
PREFIX = [r'^LINC\d', r'^MIR\d', r'^SNOR', r'^LOC\d', r'^RP\d+-']

# V2 housekeeping exact
EXACT_V2 = {
    'FAU','UBA52','RACK1','RPSA',
    'MALAT1','NEAT1','XIST','KCNQ1OT1','MEG3','MIAT','HOTAIR',
    'TSIX','FIRRE','HOTAIRM1','HCP5',
}

# V3 NEW: nuclear-encoded mitochondrial OXPHOS (MitoCarta 3.0 core OXPHOS members)
NUC_MITO = set()
# Complex I (NADH dehydrogenase) - NDUF family
for i in range(1, 14): NUC_MITO.add(f'NDUFA{i}')
NUC_MITO.add('NDUFAB1')
for i in range(1, 12): NUC_MITO.add(f'NDUFB{i}')
for i in [1,2]: NUC_MITO.add(f'NDUFC{i}')
for i in range(1, 9): NUC_MITO.add(f'NDUFS{i}')
for i in range(1, 4): NUC_MITO.add(f'NDUFV{i}')
for i in range(1, 9): NUC_MITO.add(f'NDUFAF{i}')
# Complex IV (COX)
COX_BASE = ['COX4I1','COX4I2','COX5A','COX5B',
            'COX6A1','COX6A2','COX6B1','COX6B2','COX6C',
            'COX7A1','COX7A2','COX7A2L','COX7B','COX7B2','COX7C','COX8A','COX8C',
            'COX10','COX11','COX14','COX15','COX16','COX17','COX18','COX19','COX20',
            'COXFA4','COA1','COA3','COA5','COA6','COA7']
NUC_MITO.update(COX_BASE)
# Complex III (UQCR)
NUC_MITO.update(['UQCR10','UQCR11','UQCRB','UQCRC1','UQCRC2','UQCRFS1','UQCRH','UQCRHL','UQCRQ','CYC1','CYCS'])
# Complex V (ATP synthase)
NUC_MITO.update(['ATP5F1A','ATP5F1B','ATP5F1C','ATP5F1D','ATP5F1E',
                 'ATP5MC1','ATP5MC2','ATP5MC3','ATP5MD','ATP5ME','ATP5MF',
                 'ATP5MG','ATP5MJ','ATP5MK','ATP5PB','ATP5PD','ATP5PF','ATP5PO','ATP5IF1'])
# Complex II (SDH)
NUC_MITO.update(['SDHA','SDHB','SDHC','SDHD','SDHAF1','SDHAF2','SDHAF3','SDHAF4'])
# TCA + FAO core mito enzymes (high-abundance)
NUC_MITO.update(['MDH2','ETFA','ETFB','ETFDH','ACAA2','HADHA','HADHB',
                 'CKMT1A','CKMT1B','CKMT2','IDH3A','IDH3B','IDH3G',
                 'OGDH','SUCLA2','SUCLG1','SUCLG2','FH','ACO2'])

# V3 NEW: neuronal-synaptic ambient panel (drop ONLY from non-neuronal subclasses)
AMBIENT_NEURONAL = {
    # Synaptic adhesion
    'SYT1','SYT2','NRXN1','NRXN2','NRXN3','NLGN1','NLGN2','NLGN3','NLGN4X',
    'SHANK1','SHANK2','SHANK3','PTPRD','LRRTM1','LRRTM2','LRRTM3',
    # Glutamate receptors
    'GRIN1','GRIN2A','GRIN2B','GRIN2C','GRIN2D','GRIN3A','GRIN3B',
    'GRIA1','GRIA2','GRIA3','GRIA4',
    'GRIK1','GRIK2','GRIK3','GRIK4','GRIK5',
    # GABA receptors (cortical neuron-enriched)
    'GABRA1','GABRA2','GABRA3','GABRA4','GABRA5','GABRB1','GABRB2','GABRB3',
    'GABRG1','GABRG2','GABRG3','GABBR1','GABBR2',
    # Synaptic vesicle / active zone
    'SNAP25','SYP','SYN1','SYN2','SYN3','VAMP1','VAMP2','BSN','PCLO',
    'STX1A','STX1B','STXBP1','STXBP5',
    # PSD
    'DLG4','HOMER1','HOMER2','CAMK2A','CAMK2B','CAMK2G',
    # Pan-neuron / cytoskeleton
    'MAP2','RBFOX3','NEFL','NEFM','NEFH','STMN2','TUBB3','UCHL1',
    # Neurotransmitter transporters
    'SLC17A6','SLC17A7','SLC17A8','GAD1','GAD2','SLC32A1',
    # Pan-neuron channels often leaking
    'CACNA1A','CACNA1B','CACNA1E','SCN1A','SCN2A','KCNQ2','KCNQ3',
}

_re = re.compile('|'.join(PATTERNS + SUFFIX + PREFIX))

def is_blacklisted(g, cell_group=None):
    if g in EXACT_V2 or g in NUC_MITO:
        return True
    if cell_group == 'non_neuronal' and g in AMBIENT_NEURONAL:
        return True
    return bool(_re.search(g))

def summarize():
    return {
        'patterns': len(PATTERNS) + len(SUFFIX) + len(PREFIX),
        'exact_v2': len(EXACT_V2),
        'nuc_mito': len(NUC_MITO),
        'ambient_neuronal_extra_nonneur': len(AMBIENT_NEURONAL),
    }

if __name__ == '__main__':
    import json
    print(json.dumps(summarize(), indent=2))
