seqfile = ./codon_guided_msa/codon_alignments/LOC101055557_prot_linear_multi_sl25_prt_codonAlignment.fasta
treefile = ./codon_guided_msa/raxml-ng/LOC101055557_prot_linear_multi_sl25_prt_codonAlignment.raxml.bestTree
outfile = mlc_2_LOC101055557_prot_linear_multi_sl25_prt_codonAlignment.txt

noisy = 0                         * 0,1,2,3,9: how much rubbish on the screen
verbose = 0                       * 0: concise; 1: detailed, 2: too much
runmode = 0                       * 0: user tree;  1: semi-automatic;  2: automatic
* 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

seqtype = 1                       * 1:codons; 2:AAs; 3:codons-->AAs                 <---
CodonFreq = 2                     * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table      <---

model = 0                         * models for codons:                              <---
* 0:one, 1:b, 2:2 or more dN/dS ratios for branches
* models for AAs or codon-translated AAs:
* 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
* 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

NSsites = 2
* 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
* 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
* 13:3normal>0

icode = 0                         * 0:universal code; 1:mammalian mt; 2-10:see below    <---

fix_kappa = 0                     * 1: kappa fixed, 0: kappa to be estimated                   <---
kappa = 2                         * initial or fixed kappa                                     <---
fix_omega = 0                     * 1: omega or omega_1 fixed, 0: estimate                     <---
omega = .4                        * initial or fixed omega, for codons or codon-based AAs      <---

getSE = 0                         * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0                  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

Small_Diff = .5e-6
cleandata = 0                     * remove sites with ambiguity data (1:yes, 0:no)?
fix_blength = 1                   * 0: ignore, -1: random, 1: initial, 2: fixed             <---
method = 1                        * Optimization method 0: simultaneous; 1: one branch a time

* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt.,
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt.,
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
