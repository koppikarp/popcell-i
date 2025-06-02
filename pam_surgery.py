#!/usr/bin/env python3
"""
pam_surgery.py
==============

patch_pam(window, pam_type)  →  str | None
-------------------------------------------
  window   : 3- or 6-nt **sense-strand** DNA string that constitutes the PAM
  pam_type : "NGG"  |  "CCN"
             (choose the one that the window actually shows in *sense*)

Returns a synonymous replacement (same orientation) that removes
NGG/CCN, or None if impossible.
"""

import itertools, re

# ------------------------------------------------ codon tables
aa_to_codons = {
 'F':['TTT','TTC'],'L':['TTA','TTG','CTT','CTC','CTA','CTG'],
 'I':['ATT','ATC','ATA'],'M':['ATG'],'V':['GTT','GTC','GTA','GTG'],
 'S':['TCT','TCC','TCA','TCG','AGT','AGC'],
 'P':['CCT','CCC','CCA','CCG'],'T':['ACT','ACC','ACA','ACG'],
 'A':['GCT','GCC','GCA','GCG'],'Y':['TAT','TAC'],
 'H':['CAT','CAC'],'Q':['CAA','CAG'],'N':['AAT','AAC'],
 'K':['AAA','AAG'],'D':['GAT','GAC'],'E':['GAA','GAG'],
 'C':['TGT','TGC'],'W':['TGG'],
 'R':['CGT','CGC','CGA','CGG','AGA','AGG'],
 'G':['GGT','GGC','GGA','GGG'],'STOP':['TAA','TAG','TGA']
}
codon_to_aa = {c:aa for aa,cds in aa_to_codons.items() for c in cds}
# rough usage rank (smaller is better)
rank = {c:i for aa in aa_to_codons for i,c in enumerate(aa_to_codons[aa][::-1])}

pam_ngg = re.compile(r'[ACGT]GG')
pam_ccn = re.compile(r'CC[ACGT]')

def _erase_single(codon:str,is_ngg:bool)->str|None:
    if codon=="TGG":                          # Trp cannot change
        return None
    for alt in aa_to_codons[codon_to_aa[codon]]:
        if alt==codon: continue
        if is_ngg and not pam_ngg.search(alt):      return alt
        if (not is_ngg) and not pam_ccn.search(alt):return alt
    return None

def _erase_pair(c1,c2,is_ngg):
    best, best_rank = None, 1e9
    for a,b in itertools.product(aa_to_codons[codon_to_aa[c1]],
                                 aa_to_codons[codon_to_aa[c2]]):
        test = a+b
        bad  = pam_ngg.search(test) if is_ngg else pam_ccn.search(test)
        if bad: continue
        r = rank[a]+rank[b]
        if r < best_rank:
            best,best_rank = a+b, r
    return best

def patch_pam(window:str, pam_type:str) -> str|None:
    """
    Try to synonymously eliminate the PAM in `window`.
    Return edited window (sense) or None.
    """
    window=window.upper(); pam_type=pam_type.upper()
    is_ngg = pam_type=="NGG"
    if len(window)==3:
        return _erase_single(window,is_ngg)
    elif len(window)==6:
        return _erase_pair(window[:3],window[3:],is_ngg)
    else:
        raise ValueError("window must be 3 or 6 nt")
    
if __name__ == "__main__":
    # 8 cases: 3-nt & 6-nt, NGG & CCN
    test_cases = [
        ("AGG",     "NGG"),
        ("TGG",     "NGG"),
        ("CCG",     "CCN"), 
        ("CCA",     "CCN"),  
        ("AGGTTC",  "NGG"),
        ("TGGATG",  "NGG"), 
        ("CCTCCG",  "CCN"),  
        ("TCCAGT",  "CCN"), 
    ]

    def ok(window, pam_type, out):
        patt = pam_ngg if pam_type == "NGG" else pam_ccn
        return out is None or not patt.search(out)

    for win, pam in test_cases:
        res = patch_pam(win, pam)
        status = "PASS" if ok(win, pam, res) else "FAIL"
        print(f"{win:>6} {pam} -> {res}  {status}")

