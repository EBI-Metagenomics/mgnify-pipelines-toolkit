#Protocluster Types

antiSMASH uses some abbreviations internally and in outputs to refer to the different
types of secondary metabolite clusters detected.

There are three strictness levels for cluster detection: strict, relaxed, and loose.
Each level includes the more strict levels.
The types of cluster within each ruleset are given below:

##Current Types


### Strict
|Label|Description|Added|Last updated|
|-----|-----------|:---:|:----------:|
|<span id="2dos">2dos</span>|2-deoxy-streptamine aminoglycoside|7.0|7.0|
|<span id="acyl_amino_acids">acyl_amino_acids</span>|N-acyl amino acid|4.0|4.1|
|<span id="amglyccycl">amglyccycl</span>|Aminoglycoside/aminocyclitol|&lt;= 3.0|&lt;= 3.0|
|<span id="aminocoumarin">aminocoumarin</span>|Aminocoumarin|&lt;= 3.0|&lt;= 3.0|
|<span id="aminopolycarboxylic-acid">aminopolycarboxylic-acid</span>|Aminopolycarboxylic acid metallophores (doi:10.1039/C8MT00009C)|7.0|7.0|
|<span id="archaeal-ripp">archaeal-RiPP</span>|Archaeal RiPPs (doi:10.1021/jacs.2c00521 supplemental)|7.1|7.1|
|<span id="arylpolyene">arylpolyene</span>|Aryl polyene|&lt;= 3.0|&lt;= 3.0|
|<span id="atropopeptide">atropopeptide</span>|Atropopeptide RiPPs, e.g. scabrirubin and tryptorubin|8.0 beta|8.0 beta|
|<span id="azoxy-crosslink">azoxy-crosslink</span>|axoxy compounds formed by carboxilic cross-link|8.0 beta|8.0 beta|
|<span id="azoxy-dimer">azoxy-dimer</span>|axoxy compounds formed by dimerisation|8.0 beta|8.0 beta|
|<span id="benzoxazole">benzoxazole</span>|Benzoxazoles|8.0 beta|8.0 beta|
|<span id="betalactone">betalactone</span>|Beta-lactone containing protease inhibitor|5.0|5.0|
|<span id="blactam">blactam</span>|&beta;-lactam|&lt;= 3.0|&lt;= 3.0|
|<span id="bottromycin">bottromycin</span>|Bottromycin|&lt;= 3.0|&lt;= 3.0|
|<span id="butyrolactone">butyrolactone</span>|Butyrolactone|&lt;= 3.0|&lt;= 3.0|
|<span id="cdps">CDPS</span>|tRNA-dependent cyclodipeptide synthases|5.0|8.0 beta|
|<span id="crocagin">crocagin</span>|Crocagin-like|7.0|7.0|
|<span id="cyanobactin">cyanobactin</span>|Cyanobactins like patellamide (AY986476)|&lt;= 3.0|7.0|
|<span id="cyclic-lactone-autoinducer">cyclic-lactone-autoinducer</span>|agrD-like cyclic lactone autoinducer peptides (AF001782)|6.0|6.0|
|<span id="cytokinin">cytokinin</span>|Adenine-type cytokinins, e.g. fusatin and trans-zeatin|8.0 beta|8.0 beta|
|<span id="darobactin">darobactin</span>|Darobactin-like compounds|7.1|7.1|
|<span id="deazapurine">deazapurine</span>|Deazapurine|8.0 beta|8.0 beta|
|<span id="ectoine">ectoine</span>|Ectoine|&lt;= 3.0|&lt;= 3.0|
|<span id="epipeptide">epipeptide</span>|D-amino-acid containing RiPPs such as yydF (D78193)|6.0|6.0|
|<span id="fungal_cdps">fungal_CDPS</span>|Fungal cyclodipeptide synthases|8.0 beta|8.0 beta|
|<span id="fungal-ripp">fungal-RiPP</span>|Fungal RiPP with POP or UstH peptidase types and a modification|5.0|5.0|
|<span id="furan">furan</span>|Furan|&lt;= 3.0|5.0|
|<span id="glycocin">glycocin</span>|Glycocin|&lt;= 3.0|&lt;= 3.0|
|<span id="guanidinotides">guanidinotides</span>|Pheganomycin-style protein ligase-containing cluster|4.0|6.0|
|<span id="hgle-ks"><span id="hgleks">hglE-KS</span></span>|Heterocyst glycolipid synthase-like PKS|5.0|5.0|
|<span id="HR-T2PKS">HR-T2PKS</span>|Highly reducing type II PKS like ishigamide and skyllamycin|7.1|7.1|
|<span id="hserlactone">hserlactone</span>|Homoserine lactone|&lt;= 3.0|&lt;= 3.0|
|<span id="hydrogen-cyanide">hydrogen-cyanide</span>|Hydrogen cyanide (AF208523, doi:10.1128/jb.182.24.6940-6949.20)|7.1|7.1|
|<span id="hydroxytropolone">hydroxy-tropolone</span>|7-hydroxytropolone-like cluster|8.0 beta|8.0 beta|
|<span id="HR-T2PKS">HR-T2PKS</span>|Highly reducing type II PKS like ishigamide and skyllamycin|7.1|7.1|
|<span id="indole">indole</span>|Indole|&lt;= 3.0|4.0|
|<span id="isocyanide">isocyanide</span>|Isocyanides (doi:10.1093/nar/gkad573)|7.1|7.1|
|<span id="isocyanide-nrp">NRP with isocyanide</span>|Isocyanides (doi:0.1128/mBio.00785-18)|7.1|7.1|
|<span id="ladderane">ladderane</span>|Ladderane|&lt;= 3.0|7.1|
|<span id="lanthipeptide-class-i">lanthipeptide class I</span>|Class I lanthipeptides like nisin|4.2|8.0 beta|
|<span id="lanthipeptide-class-ii">lanthipeptide class II</span>|Class II lanthipeptides like mutacin II (U40620)|4.2|6.0|
|<span id="lanthipeptide-class-iii">lanthipeptide class III</span>|Class III lanthipeptides like labyrinthopeptin (FN178622)|4.2|6.0|
|<span id="lanthipeptide-class-iv">lanthipeptide class IV</span>|Class IV lanthipeptides like venezuelin (HQ328852)|4.2|6.0|
|<span id="lanthipeptide-class-v">lanthipeptide class V</span>|Glycosylated lanthipeptide/linaridin hybrids like MT210103|5.1|6.0|
|<span id="lassopeptide">lassopeptide</span>|Lasso peptide|&lt;= 3.0|5.0|
|<span id="leupeptin">leupeptin</span>|leupeptin-like compounds|8.0 beta|8.0 beta|
|<span id="linaridin">linaridin</span>|Linear arid peptide such as cypemycin (HQ148718) and salinipeptin (MG788286)|&lt;= 3.0|&lt;= 3.0|
|<span id="lincosamides">lincosamides</span>|NRPS-adjacent biosynthesis of lincosamides|8.0 beta|8.0 beta|
|<span id="lipolanthine">lipolanthine</span>|Lanthipeptide class containing N-terminal fatty acids such as MG673929|5.0|5.0|
|<span id="melanin">melanin</span>|Melanin|&lt;= 3.0|&lt;= 3.0|
|<span id="methanobactin">methanobactin</span>|Copper-chelating/transporting peptides (doi:10.1126/science.aap9437)|7.0|7.0|
|<span id="microviridin">microviridin</span>|Microviridin|&lt;= 3.0|&lt;= 3.0|
|<span id="mycosporine">mycosporine</span>|Molecules containing mycosporine-like amino acid|7.0|8.0 beta|
|<span id="naggn">NAGGN</span>|N-acetylglutaminylglutamine amide|5.0|5.0|
|<span id="napaa">NAPAA</span>|Non-alpha poly-amino acids like e-Polylysin|6.0|7.0|
|<span id="ni-siderophore"><span id="nrp-independent-siderophore">NI-siderophore</span></span>|NRPS-independent, IucA/IucC-like siderophores (*siderophore* prior to 7.0)|7.0|8.0 beta|
|<span id="nitropropanoic-acid">nitropropanoic-acid</span>|3-Nitropropanoic acid (neurotoxin)|8.0 beta|8.0 beta|
|<span id="nrps">NRPS</span>|Non-ribosomal peptide synthetase|&lt;= 3.0|8.0 beta|
|<span id="nrp-metallophore">NRP-metallophore</span>|Non-ribosomal peptide metallophores|7.0|7.0|
|<span id="nucleoside">nucleoside</span>|Nucleoside|&lt;= 3.0|5.0|
|<span id="oligosaccharide">oligosaccharide</span>|Oligosaccharide|&lt;= 3.0|&lt;= 3.0|
|<span id="opine-like-metallophore">opine-like-metallophore</span>|Opine-like zincophores like staphylopine (doi:10.1128/mSystems.00554-20)|7.0|7.0|
|<span id="other">other</span>|Cluster containing a secondary metabolite-related protein that does not fit into any other category|4.0|5.0|
|<span id="pbde">PBDE</span>|Polybrominated diphenyl ether|4.1|4.1|
|<span id="phenazine">phenazine</span>|Phenazine|&lt;= 3.0|&lt;= 3.0|
|<span id="phosphoglycolipid">phosphoglycolipid</span>|Phosphoglycolipid|&lt;= 3.0|&lt;= 3.0|
|<span id="phosphonate">phosphonate</span>|Phosphonate|&lt;= 3.0|7.0|
|<span id="polyhalogenated-pyrrole">polyhalogenated-pyrrole</span>|Polyhalogenated pyrrole|8.0 beta|8.0 beta|
|<span id="polyyne">polyyne</span>|Polyyne|8.0 beta|8.0 beta|
|<span id="ppys-ks"><span id="ppysks">PpyS-KS</span></span>|PPY-like pyrone|4.2|4.2|
|<span id="prodigiosin">prodigiosin</span>|Serratia-type non-traditional PKS prodigiosin biosynthesis pathway|6.0|6.0|
|<span id="proteusin">proteusin</span>|Proteusin|&lt;= 3.0|&lt;= 3.0 |
|<span id="pufa">PUFA</span>|Polyunsaturated fatty acid|&lt;= 3.0|&lt;= 3.0|
|<span id="pyrrolidine">pyrrolidine</span>|Pyrrolidines like described in BGC0001510|6.0|6.0|
|<span id="ranthipeptide">ranthipeptide</span>|Cys-rich peptides (aka. SCIFF: six Cys in fourty-five) like in CP001581:3481278-3502939|6.0|6.0|
|<span id="ras-ripp">RaS-RiPP</span>|Streptide-like thioether-bond RiPPs|5.0|5.0|
|<span id="rcdps">RCDPS</span>|Fungal Arginine-containing cyclic dipeptides|7.0|7.0|
|<span id="redox-cofactor">redox-cofactor</span>|Redox-cofactors such as PQQ (NC_021985:1458906-1494876)|6.0|6.0|
|<span id="resorcinol">resorcinol</span>|Resorcinol|&lt;= 3.0|&lt;= 3.0|
|<span id="sactipeptide">sactipeptide</span>|Sactipeptide|&lt;= 3.0|6.0|
|<span id="spliceotide">spliceotide</span>|RiPPs containing plpX type spliceases (NZ_KB235920:17899-42115)|6.0|6.0|
|<span id="t1pks">T1PKS</span>|Type I PKS (Polyketide synthase)|&lt;= 3.0|&lt;= 3.0|
|<span id="t2pks">T2PKS</span>|Type II PKS|&lt;= 3.0|5.0|
|<span id="t3pks">T3PKS</span>|Type III PKS|&lt;= 3.0|&lt;= 3.0|
|<span id="terpene">terpene</span>|Terpene|&lt;= 3.0|8.0 beta|
|<span id="thioamitides">thioamitides</span>|Thioamitide RiPPs as found in JOBF01000011|5.1|6.0|
|<span id="thioamide-nrp">thioamide-NRP</span>|Thioamide-containing non-ribosomal peptide|5.0|5.0|
|<span id="transatpks"><span id="transat-pks">transAT-PKS</span></span>|Trans-AT PKS|&lt;= 3.0|8.0 beta|
|<span id="triceptide">triceptide</span>|Triceptides|7.1|7.1|
|<span id="tropodithietic-acid">tropodithietic-acid</span>|Tropodithietic acid|5.0|5.0|

### Relaxed
|Label|Description|Added|Last updated|
|-----|-----------|:---:|:----------:|
|<span id="fungal-ripp-like">fungal-RiPP-like</span>|Fungal RiPP-likes|7.0|7.0|
|<span id="nrps-like"><span id="nrpsfragment">NRPS-like</span></span>|NRPS-like fragment|5.0|8.0 beta|
|<span id="phosphonate-like">phosphonate-like</span>|Phosphonate-like (prior to 7.0 this was the phosphonate rule)|7.0|7.0|
|<span id="pks-like">PKS-like</span>|Other types of PKS|5.0|5.0|
|<span id="ripp-like">RiPP-like</span>|Other unspecified ribosomally synthesised and post-translationally modified peptide product (RiPP)|4.1|7.0|
|<span id="rre-containing">RRE-containing</span>|RRE-element containing cluster|6.0|6.0|
|<span id="terpene">terpene-precursor</span>|Compound likely used as a terpene precursor|8.0 beta|8.0 beta|
|<span id="transatpks-like"><span id="transat-pks-like">transAT-PKS-like</span></span>|Trans-AT PKS fragment, with trans-AT domain not found|&lt;= 5.0|5.0|

### Loose
|Label|Description|Added|Last updated|
|-----|-----------|:---:|:----------:|
|<span id="fatty_acid">fatty_acid</span>|Fatty acid (loose strictness, likely from primary metabolism)|&lt;= 3.0|4.2|
|<span id="halogenated">halogenated</span>|Halogenase-containing cluster, potentially generating a halogenated product|5.0|5.0|
|<span id="lysine">lysine</span>|Fungal lysine primary metabolism|8.0 beta|8.0 beta|
|<span id="saccharide">saccharide</span>|Saccharide (loose strictness, likely from primary metabolism)|&lt;= 3.0|&lt;= 3.0|


##Obsolete/Previous Types

|Label|Description|Added|Removed|Notes|
|-----|-----------|:---:|:-----:|-----|
|<span id="lap">LAP</span>|Linear azol(in)e-containing peptides|&lt;= 3.0|8.0 beta|Merged with *thiopeptide* into *azole-containing RiPP*|
|<span id="mycosporine-like">mycosporine-like</span>|Molecules containing mycosporine-like amino acid|7.0|8.0 beta|Renamed to *mycosporine*|
|<span id="thiopeptide">thiopeptide</span>|Thiopeptide|4.2|8.0 beta|Merged with *LAP* into *azole-containing RiPP*|
|<span id="siderophore">siderophore</span>|Siderophore|&lt;= 3.0|7.0|Renamed to *NI-siderophore*|
|<span id="bacteriocin">bacteriocin</span>|Bacteriocin or other unspecified ribosomally synthesised and <br/> post-translationally modified peptide product (RiPP)|4.1|6.0|Renamed to *RiPP-like*|
|<span id="fused">fused</span>|Pheganomycin-style protein ligase-containing cluster|4.0|6.0|Renamed to *guanidinotides*|
|<span id="head_to_tail">head_to_tail</span>|Head-to-tail cyclised RiPP (subtilosin-like)|&lt;= 3.0|6.0|Merged into *sactipeptides*|
|<span id="lanthidin">lanthidin</span>|Glycosylated lanthipeptide/linaridin hybrids like MT210103|5.1|6.0|Renamed to *lanthipeptide class V*|
|<span id="lantipeptide"><span id="lanthipeptide">lanthipeptide</span></span>|Lanthipeptides|4.2|6.0|Split into subclasses, e.g. *lanthipeptide class I*|
|<span id="tfua-related">TfuA-related</span>|TfuA-related RiPPs|5.1|6.0|Renamed to *thioamitides*|
|<span id="otherks">otherks</span>|Other types of PKS|&lt;= 3.0|5.0|Split into *PKS-like* and *hglE-KS*|
|<span id="microcin">microcin</span>|Microcin|&lt;= 3.0|5.0|Merged into *lasso peptide*|
|<span id="cf_saccharide">cf_saccharide</span>|Possible saccharide|&lt;= 3.0|5.0|Renamed to *saccharide*|
|<span id="cf_fatty_acid">cf_fatty_acid</span>|Possible fatty acid|&lt;= 3.0|5.0|Renamed to *fatty_acid*|
|<span id="cf_putative">cf_putative</span>|Putative cluster of unknown type identified <br/> with the ClusterFinder algorithm|&lt;= 3.0|5.0|These predictions are now subregions, not clusters|

