If you want to calculate GERP scores (http://mendel.stanford.edu/SidowLab/downloads/gerp/) across a phylogeny for a range of genomic features, then included in this rough repo is a means of doing so. This was used to derive GERP scores for mouse CTCF, enhancers and promoters as part of the Liver Cancer Evolution consortium, many details of which can be found here:https://www.nature.com/articles/s41586-020-2435-1

WARNING: do not even attempt to use this if you don't have progressiveCactus running or if you're interested in whole genome measures as it'll take an extremely long time to run.

Finally, always check the output. Is it sensible? Is anything missing? Did your jobs run properly?

Supplied environment in .yml file with key package version, includes pandas, bedtools and agat, which is only necessary if you want your own 4D tree. 
<pre>conda activate gerpit</pre>

variable you'll need to have.
<pre>#INPUT_FEATURES=$1 #input bedfile to calculate gerp scores across
#PREFIX=$2 #opportunity to shorten file handle
#DIR=$3 #working directory
#BIN=$4 #location of these files
#GENOME_NAME=$5 # target genome name within hal alignment
#PCE=$6 #progressiveCactus Environment path
#HAL=$7 # location of hal alignment
#GERP=$8 #gerp bin directory
</pre>

Run 4D_tree.sh here. Note reduced input variables necessary
<pre>bash 4D_tree.sh $DIR $PCE $HAL
</pre>

=======

The following should really be run interactively so you can get a feel for times and how best to partition features into smaller chunks with respect to your resources.

<pre>##generate feature list. $1:number of features, $2:output file
mkdir ${DIR}/output

FEAT_COUNT=$(wc -l ${INPUT_FEATURES} | cut -d" " -f1)
</pre>

make index files to break up bed into manageable jobs. 50 features should only take an hour or so, change the 50 below if you are more or less patient.
<pre>paste <(seq 1 50 ${FEAT_COUNT}) <(seq 50 50 ${FEAT_COUNT}) | sed "$ s/$/${FEAT_COUNT}/" > ${DIR}/output/${PREFIX}_list.txt
while read line
 do bits=($line)
 source ${PCE} 
 cd ${DIR}
 for y in {${bits[0]}..${bits[1]}} 
  do ${BIN}/gerpit.sh \$y ${INPUT_FEATURES} ${GENOME_NAME} ${DIR}/output ${PREFIX}_${bits[0]} ${HAL} ${GERP}
 done
done < ${DIR}/output/${PREFIX}_list.txt
unset FEAT_COUNT
</pre>

to run on a cluster, this makes the most sense:
<pre>while read line
 do bits=($line)
 bsub -M 4000 -R "rusage[mem=4000]" -W 42:00 -o pp_${x}.o -e pp_${x}.e "source ${PCE} 
 cd ${DIR}
 for y in {${bits[0]}..${bits[1]}} 
  do ${BIN}/gerpit.sh \$y ${INPUT_FEATURES} ${GENOME_NAME} ${DIR}/output ${PREFIX}_${bits[0]} ${HAL} ${GERP}
 done "
done < ${DIR}/output/${PREFIX}_list.txt
</pre>

Once these are finished, deal with N's and calculate mean GERP scores for each feature
<pre>##this puts aportioned features together and then splits by chromomsome so sorting and calculations are faster: -M 40000 -W 10:00
cat ${DIR}/output/${PREFIX}*.ind.bed > ${DIR}/output/${PREFIX}.all.ind.bed

##split into chromosomes 1-19 & X. If missing data is from N's and are dealt with by refering to the original feature
for x in $(seq 1 1 19) X
 do cat ${DIR}/output/${PREFIX}.all.ind.bed | grep -P -e "^${x}[ \t]" | awk '{if (NF == 5) print $0}' | sort -k1,1V -k2,2g | uniq > ${DIR}/output/${PREFIX}.all.${x}.ind.bed
done

##split up original by chromosome
for x in $(seq 1 1 19) X
 do cat ${INPUT_FEATURES} | grep -P -e "^${x}[ \t]" > ${DIR}/output/original.${x}.bed
done

##calculate the mean gerp score for a feature, excluding any regions that contain N's
for x in $(seq 1 1 19) X 
 do grep -v NA ${DIR}/output/${PREFIX}.${x}.ind.bed > ${DIR}/output/${PREFIX}.${x}.ind.noNA.bed
 while read line ; do bits=($line) ; DIST=$((${bits[2]}-${bits[1]})) ; echo ${line} | sed 's/ /\t/g' | bedtools intersect -a stdin -b ${DIR}/output/${PREFIX}.${x}.ind.noNA.bed -wb | cut -f7- > ${DIR}/output/${PREFIX}.${x}.tmp
 #the following will derive the merged feature coordinates, reporting the size without N's and the extreme edges of the feature (so if it's split in two by any number of N's, we still report the new feature coordinates so that it includes the string of N's. The assumption is that N's in the middle will be very limited and it won't be worth splitting the feature, which is the correct way to deal with this if stretches of N's are more extensive. 
 bedtools merge -i ${DIR}/output/${PREFIX}.${x}.tmp > ${DIR}/output/${PREFIX}.${x}.merge.tmp
 NEW_FEAT_DIST=$(awk '{print $1"\t"$2"\t"$3"\t"$3-$2}' ${DIR}/output/${PREFIX}.${x}.merge.tmp | awk '{SUM+=$4}END{print SUM}')
 START=$(head -n 1 ${DIR}/output/${PREFIX}.${x}.merge.tmp | awk '{print $1" "$2}')
 END=$(tail -n 1 ${DIR}/output/${PREFIX}.${x}.merge.tmp | awk '{print $3}')
 NEW_FEAT=$(echo ${START} ${END})
 MEAN=$(awk '{SUM+=$5}END{print SUM}' ${DIR}/output/${PREFIX}.${x}.tmp | X=${NEW_FEAT_DIST} awk '{printf("%.6f\n", $1/ENVIRON["X"])}')
 echo ${NEW_FEAT} ${MEAN} ${NEW_FEAT_DIST} ${DIST} $(($DIST-$NEW_FEAT_DIST)) $line | sed 's/ /\t/g' >> ${DIR}/output/original.${x}.out.ind.bed
 rm ${DIR}/output/${PREFIX}.${x}.tmp
 rm ${DIR}/output/${PREFIX}.${x}.merge.tmp
 done < ${DIR}/output/original.${x}.bed
done

rm ${DIR}/output/${PREFIX}.${x}.ind.noNA.bed

#put together all parts of genome, which will be sorted
for x in $(seq 1 1 19) X
 do cat ${DIR}/output/${PREFIX}.${x}.out.ind.bed >> ${DIR}/output/${PREFIX}.mean.bed
 cat ${DIR}/output/${PREFIX}.${x}.ind.bed >> ${DIR}/output/${PREFIX}.ind.bed
done

#I've not removed intermediary files, as these can be interesting or necessary to debug, which I've put zero effort into sourcing.
</pre>

Example output on a per nucleotide level "ind.bed" file:
<pre>
#chr, start, end, neutral, GERP score
1       59170222        59170223        1.64    0.422
1       59170223        59170224        1.64    -1.35
1       59170224        59170225        1.64    1.64
1       59177592        59177593        0.4     0.4
1       59177593        59177594        0.4     -0.8
1       59177594        59177595        0.4     0.4
</pre>

n.b. N's will appear as NA in columns 4 and 5.

Example output on a per feature level "mean.bed" file:
<pre>
#chr, start(accounting for N's), end(accounting for N's), mean GERP score, length without N's, length original, original chr, original start, original end, metadata....
1       89229866        89230254        -0.132018       388     415     27      1       89229866        89230281        c3h_shadowCTCF_224      caroli_CTCF_754 shadowCTCF
1       89247203        89247638        -0.139228       435     435     0       1       89247203        89247638        c3h_shadowCTCF_225      caroli_CTCF_755 shadowCTCF
</pre>

