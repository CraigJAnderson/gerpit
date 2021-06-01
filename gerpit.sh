#!/bin/bash
#gerpit.sh
JOB=$1
BED=$2
SPE=$3
OUTDIR=$4
HANDLE=$5
HAL=$6
GERP=$7
##make space for job
mkdir ${OUTDIR}/${HANDLE}_${JOB}
cd ${OUTDIR}/${HANDLE}_${JOB}
sed -n ${JOB}p ${BED} > ${JOB}.bed
##get the maf from the hal file, specific to the bed coordinates
hal2maf --refGenome ${SPE} --refTargets ${JOB}.bed --noDupes --noAncestors --global --onlyOrthologs --unique ${HAL} ${JOB}.maf

##strip the header, then the tree for gerp++
sed -i 2d ${JOB}.maf

##the following will stitch together the individual alignment fragments from a maf produced by hal2maf. Care is taken to infer and fill in gaps not reported in the maf, which would otherwise lead to an inaccurate alignment.
egrep -nr "^$" ${JOB}.maf | sed 's/://g' | awk '{ print prev"\t"$0} { prev = $0 }'| sed 1d | awk '{print $1+2","$2-1}' > list.txt
for x in $(seq 1  1 $(wc -l list.txt | cut -f1 -d" ")) ; do TMP=$(sed -n ${x}p list.txt) ; sed -n ${TMP}p ${JOB}.maf | sed 's/\./\t/g' | awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > maf_part${x}.txt ; done
BREAKS=$(wc -l list.txt)
python maf2fa.py $BREAKS
sed -e 's/^/>/g' prep.txt | sed -e 's/\t/\n/g' > out.fa

##run gerp,append to bed coordinates and calculate neutral and RS averages for all_mean. Output only the first elem for all_el.gerp
${GERP}/gerpcol -a -f out.fa -t ${DIR}/src/tree.txt -e ${SPE} -j

##output scores for individual sites
while read line ; do bits=($line) ; for x in $(seq ${bits[1]} 1 ${bits[2]}) ; do echo $x | Y=${bits[0]} awk '{print ENVIRON["Y"]"\t"$1"\t"$1+1}' >> tmp.bed ; done ; sed '$d' tmp.bed | paste - out.fa.rates >> ${OUTDIR}/${HANDLE}.ind.bed ; done < ${JOB}.bed

##unhash if you want elements output
#${GERP}/gerpelem -f ./out.fa.rates
#REG=$(wc -l out.fa.rates.elems | awk -F" " '{print $1}')
#for x in $(seq 1 1 $REG) ; do cat ${JOB}.bed >> rates.bed ; done
#paste rates.bed out.fa.rates.elems | awk '{print $1"\t"$2+$6"\t"$2+$7"\t"$4"\t"$8"\t"$9}' | sort -k 1,1 -k2,2n  > elems.txt
#sed -n ${JOB}p ${BED} > ${JOB}_2.bed
#paste rates.bed out.fa.rates.elems | awk '{print $1"\t"$2+$6"\t"$2+$7"\t"$4"\t"$8"\t"$9}' | bedtools intersect -a stdin -b ${JOB}_2.bed -wa | sort -nk 6,6 | head -n 1 >> ${OUTDIR}/${HANDLE}.top_el.gerp
#cat elems.txt >> ${OUTDIR}/${HANDLE}.all_el.gerp

#tidy up
cd ..
rm -r ${HANDLE}_${JOB}
