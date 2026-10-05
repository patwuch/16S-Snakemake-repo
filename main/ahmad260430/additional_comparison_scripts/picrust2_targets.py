"""Curated PICRUSt2 features for ahmad260430's two targeted panels: inosine and isoferulic acid.

Fixed from the feature annotations alone (picrust2/described), before any per-feature result
was read (2026-10-04), as reza250909's SCFA panels were. Each entry is
(database, feature, role, label).

Inosine: the routes that make inosine (from adenosine or IMP), break it down to hypoxanthine,
salvage it back to IMP, take it up, and carry hypoxanthine on, plus the MetaCyc pathways built on
those steps.

Isoferulic acid (3-hydroxy-4-methoxycinnamic acid): no KO, EC or MetaCyc feature is specific to it.
The list holds the bacterial routes annotated for hydroxycinnamic acids in general (ferulate,
caffeate, cinnamate): decarboxylation, CoA activation and side-chain cleavage to vanillin,
vanillin/vanillate demethylation and decarboxylation, and ring cleavage of the protocatechuate and
catechol products. Whether these enzymes act on the isomer is not known from the annotations.
"""

INOSINE = "inosine"
ISOFERULIC = "isoferulic"

I_FORM = "Inosine formation (from adenosine or IMP)"
I_BREAK = "Inosine breakdown to hypoxanthine"
I_SALV = "Salvage to IMP"
I_UPT = "Nucleoside uptake"
I_ON = "Hypoxanthine onward (xanthine, urate)"
I_PWY = "Pathways through inosine"

F_DEC = "Decarboxylation of hydroxycinnamic acids"
F_COA = "CoA activation and side-chain cleavage (ferulate to vanillin)"
F_RED = "Hydroxycinnamoyl-CoA reduction and methylation"
F_VAN = "Vanillin and vanillate onward (demethylation, decarboxylation)"
F_RING = "Ring cleavage of protocatechuate and catechol"
F_HCA = "Cinnamate and phenylpropanoate dioxygenase"

FEATURES = {
    INOSINE: [
        ("MetaCyc", "PWY-6609", I_PWY, "adenine and adenosine salvage III"),
        ("MetaCyc", "PWY-6353", I_PWY, "purine nucleotides degradation II (aerobic)"),
        ("MetaCyc", "PWY-6608", I_PWY, "guanosine nucleotides degradation III"),
        ("MetaCyc", "PWY-5695", I_PWY, "inosine 5'-phosphate degradation"),
        ("MetaCyc", "P164-PWY", I_PWY, "purine nucleobases degradation I (anaerobic)"),
        ("MetaCyc", "PWY-5497", I_PWY, "purine nucleobases degradation II (anaerobic)"),
        ("KO", "ko:K01488", I_FORM, "add adenosine deaminase"),
        ("KO", "ko:K01081", I_FORM, "5'-nucleotidase (K01081)"),
        ("KO", "ko:K03787", I_FORM, "surE 5'-nucleotidase"),
        ("KO", "ko:K08723", I_FORM, "5'-nucleotidase (K08723)"),
        ("KO", "ko:K11751", I_FORM, "ushA 5'-nucleotidase / UDP-sugar diphosphatase"),
        ("KO", "ko:K20881", I_FORM, "5'-nucleotidase (K20881)"),
        ("KO", "ko:K03783", I_BREAK, "punA purine-nucleoside phosphorylase"),
        ("KO", "ko:K03784", I_BREAK, "deoD purine-nucleoside phosphorylase"),
        ("KO", "ko:K09913", I_BREAK, "ppnP purine/pyrimidine-nucleoside phosphorylase"),
        ("KO", "ko:K01239", I_BREAK, "iunH purine nucleosidase"),
        ("KO", "ko:K00892", I_SALV, "gsk inosine kinase"),
        ("KO", "ko:K00760", I_SALV, "hpt hypoxanthine phosphoribosyltransferase"),
        ("KO", "ko:K03317", I_UPT, "nupC concentrative nucleoside transporter"),
        ("KO", "ko:K16323", I_UPT, "purine nucleoside transport protein"),
        ("KO", "ko:K11535", I_UPT, "nucleoside transport protein"),
        ("KO", "ko:K00087", I_ON, "xdhA xanthine dehydrogenase, molybdenum subunit"),
        ("KO", "ko:K13479", I_ON, "xanthine dehydrogenase, FAD subunit"),
        ("KO", "ko:K13481", I_ON, "xanthine dehydrogenase, small subunit"),
        ("KO", "ko:K13482", I_ON, "xanthine dehydrogenase, large subunit"),
        ("EC", "3.5.4.4", I_FORM, "adenosine deaminase"),
        ("EC", "3.1.3.5", I_FORM, "5'-nucleotidase"),
        ("EC", "2.4.2.1", I_BREAK, "purine-nucleoside phosphorylase"),
        ("EC", "3.2.2.1", I_BREAK, "purine nucleosidase"),
        ("EC", "2.7.1.73", I_SALV, "inosine kinase"),
        ("EC", "2.4.2.8", I_SALV, "hypoxanthine phosphoribosyltransferase"),
        ("EC", "1.17.1.4", I_ON, "xanthine dehydrogenase"),
    ],
    ISOFERULIC: [
        ("MetaCyc", "PROTOCATECHUATE-ORTHO-CLEAVAGE-PWY", F_RING, "protocatechuate degradation II (ortho-cleavage)"),
        ("MetaCyc", "PWY-5431", F_RING, "aromatic compounds degradation via beta-ketoadipate"),
        ("MetaCyc", "CATECHOL-ORTHO-CLEAVAGE-PWY", F_RING, "catechol degradation to beta-ketoadipate"),
        ("KO", "ko:K13727", F_DEC, "padC phenolic acid decarboxylase"),
        ("KO", "ko:K12508", F_COA, "fcs feruloyl-CoA synthase"),
        ("KO", "ko:K18383", F_COA, "ech feruloyl-CoA hydratase/lyase"),
        ("KO", "ko:K00588", F_RED, "caffeoyl-CoA O-methyltransferase"),
        ("KO", "ko:K13607", F_RED, "fldA cinnamoyl-CoA:phenyllactate CoA-transferase"),
        ("KO", "ko:K21802", F_VAN, "vdh vanillin dehydrogenase"),
        ("KO", "ko:K03862", F_VAN, "vanA vanillate monooxygenase"),
        ("KO", "ko:K03863", F_VAN, "vanB vanillate monooxygenase, ferredoxin subunit"),
        ("KO", "ko:K21759", F_VAN, "vanillate/4-hydroxybenzoate decarboxylase, subunit D"),
        ("KO", "ko:K00448", F_RING, "pcaG protocatechuate 3,4-dioxygenase, alpha"),
        ("KO", "ko:K00449", F_RING, "pcaH protocatechuate 3,4-dioxygenase, beta"),
        ("KO", "ko:K04100", F_RING, "ligA protocatechuate 4,5-dioxygenase, alpha"),
        ("KO", "ko:K04101", F_RING, "ligB protocatechuate 4,5-dioxygenase, beta"),
        ("KO", "ko:K05710", F_HCA, "hcaC phenylpropionate/cinnamate dioxygenase, ferredoxin"),
        ("KO", "ko:K00529", F_HCA, "hcaD phenylpropionate/cinnamate dioxygenase, reductase"),
        ("EC", "6.2.1.34", F_COA, "trans-feruloyl-CoA synthase"),
        ("EC", "4.2.1.101", F_COA, "trans-feruloyl-CoA hydratase"),
        ("EC", "4.1.2.41", F_COA, "vanillin synthase"),
        ("EC", "1.3.1.108", F_RED, "caffeoyl-CoA reductase"),
        ("EC", "2.1.1.104", F_RED, "caffeoyl-CoA O-methyltransferase"),
        ("EC", "1.2.1.67", F_VAN, "vanillin dehydrogenase"),
        ("EC", "1.14.13.82", F_VAN, "vanillate monooxygenase"),
        ("EC", "4.1.1.61", F_VAN, "4-hydroxybenzoate decarboxylase"),
        ("EC", "1.13.11.3", F_RING, "protocatechuate 3,4-dioxygenase"),
        ("EC", "1.13.11.8", F_RING, "protocatechuate 4,5-dioxygenase"),
    ],
}
TITLES = {INOSINE: "Inosine", ISOFERULIC: "Isoferulic acid (hydroxycinnamate routes)"}
