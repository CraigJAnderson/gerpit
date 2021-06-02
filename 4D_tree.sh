DIR=$1 #working directory
PCE=$2 #progressiveCactus Environment path
HAL=$3 # location of hal alignment

##to begin, you'll need to generate a neutral phylogeny of your hal alignment for use with gerp. If you don't care about accuracy, just generate a tree with the following hal command:
#halStats 1509_ca.hal --tree

##otherwise, below is a way to derive a tree from 4D sites across your phylogeny from example mouse gff from ensembl 100.  In short, I pull out the coding sequence from a specific flavour of gene model "ensembl_havana" and export as a bed. I'm getting gene models and resources for successful completion: -W 10:00 -M 15000
mkdir ${DIR}/src && cd ${DIR}/src

wget ftp://ftp.ensembl.org/pub/release-100/gff3/mus_musculus/Mus_musculus.GRCm38.100.gff3.gz
grep ensembl_havana Mus_musculus.GRCm38.100.gff3 > Mus_musculus.GRCm38.100_2.gff3
agat_convert_sp_gxf2gxf.pl -g Mus_musculus.GRCm38.100_2.gff3 -o Mus_musculus.GRCm38.100_3.gff3
awk '{if ($3 =="CDS") print $0}' Mus_musculus.GRCm38.100_3.gff3 > Mus_musculus.GRCm38.100_4.gff3
agat_convert_sp_gff2bed.pl --gff Mus_musculus.GRCm38.100_4.gff3 -o Mus_musculus.GRCm38.100_4.bed

##Get 4D sites with respect to the annotation bed just made. Takes hours to complete for --conserved, but is very quick for a single genome (-W 20:00 -M 15000)
source ${PCE}
hal4dExtract --conserved ${HAL} C57B6J Mus_musculus.GRCm38.100_4.bed Mus_musculus.GRCm38.100_4D.bed

##convert to the more common bed format
bed12ToBed6 -i Mus_musculus.GRCm38.100_4D.bed > Mus_musculus.GRCm38.100_5.bed

##I want a bed of each 4D site on only autosomes and the X chromosome
awk '{if ($1 != "Y") print $0}' Mus_musculus.GRCm38.100_5.bed > Mus_musculus.GRCm38.100_6.bed

##here I use phylofit to estimate the phylogenetic tree to provide gerp. Took a while: -W 50:00 -M 15000 and I need to supply a basic binary tree, which I manually edited down from the original hal alignment tree.
#(loxAfr3,(((canFam3,felCat8),(oviAri3,bosTau8)),(rheMac3,(oryCun2,(jacJac1,(micOch1,(Rattus,(CAROLI_EiJ,(CAST_EiJ,(C57B6J,C3H_HeJ))))))))));

hal2maf --onlyOrthologs --noAncestors --noDupes --unique --refGenome C57B6J --refTargets Mus_musculus.GRCm38.100_6.bed 1509_ca.hal mm10_4D.maf ; phyloFit --tree \"(loxAfr3,(((canFam3,felCat8),(oviAri3,bosTau8)),(rheMac3,(oryCun2,(jacJac1,(micOch1,(Rattus,(CAROLI_EiJ,(CAST_EiJ,(C57B6J,C3H_HeJ))))))))));\" --msa-format MAF --out-root mm10_4D.mod mm10_4D.maf

##output for neutral tree in #mm10_4D.mod
#ALPHABET: A C G T
#ORDER: 0
#SUBST_MOD: REV
#TRAINING_LNL: -23887879.554310
#BACKGROUND: 0.221286 0.277973 0.278516 0.222225
#RATE_MAT:
#  -1.081956    0.218350    0.692751    0.170855
#   0.173821   -0.937146    0.212741    0.550583
#   0.550401    0.212326   -0.936489    0.173761
#   0.170133    0.688703    0.217776   -1.076612
#TREE: (loxAfr3:0.0994233,(((canFam3:0.111048,felCat8:0.0971083):0.0516944,(oviAri3:0.0393319,bosTau8:0.0340846):0.151127):0.029594,(rheMac3:0.133016,(oryCun2:0.217621,(jacJac1:0.193543,(micOch1:0.129733,(Rattus:0.0890185,(CAROLI_EiJ:0.0217885,(CAST_EiJ:0.00448278,(C57B6J:0.00143669,C3H_HeJ:0.00146719):0.00389095):0.0186485):0.0645882):0.0666358):0.123377):0.0812895):0.0177262):0.0221412):0.0994233);

#store tree
tail -n 1 mm10_4D.mod | sed 's/TREE: //g' > tree.txt
