#!user/loci/bin/bash
#Three files are needed 1, EDTA's TEanno.gff3 2, RepeatsMask's masked.out 3, and the result format of reference sequence ref 4,SV is #Chr#Start#End#SVID#type#Leng 6 columns
# Example SV file [without header] HuBei_LG_01	13888	13941	DEL00000003	DEL	54

#EDTAresult
EDTA_result=${1}  #for example: lyc_gao_hube_clean_assembly.FINAL.fasta.mod.EDTA.TEanno.gff3
#RepeatsMaskresult
RepeatsMask_result=${2} #for example:  lyc_gao_hube_clean_assembly.FINAL.fasta.mod.MAKER.masked.out
#ref
ref=${3}
#SV file
SV_file=${4}


#TE result
cat ${EDTA_result} | grep -v '#' | awk 'BEGIN{FS="[\t;]";OFS="\t"}{print $1,$4,$5,$5-$4+1,$11}' > TMP_TE.bed

#Tandem repeats result
cat ${RepeatsMask_result} | sed '1,3d' | awk 'BEGIN{FS=" ";OFS="\t"}{print $5,$6,$7,$7-$6+1,$11}' > TMP_Simple_repeats.bed

bedtools subtract -a TMP_Simple_repeats.bed -b TMP_TE.bed -A -F 0.1 -f 0.1  > TMP-Simple_repeats2.bed

#Chrom leng
perl fasta_size_report.pl -f ${ref} | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2 }' > TMP-CHROM-leng.txt


cat ${SV_file} | sed '1d' > TMP-SV.txt

bedtools intersect -f 0.5 -wb -a TMP-SV.txt -b TMP-Simple_repeats2.bed > VNTRs.result.txt

bedtools intersect -f 0.5 -wb -a TMP-SV.txt -b TMP-Simple_repeats2.bed | awk 'BEGIN{OFS="\t";FS="\t"}{print $4}' | sort -n | uniq > TMP.unique.VNTRs.name.txt

cp TMP.unique.VNTRs.name.txt VNTRs.final.txt

#removing VNTRs
cp TMP-SV.txt TMP-SV-removeVNTRs.txt
cat TMP.unique.VNTRs.name.txt | while read y ; do sed -i "/$y/d" TMP-SV-removeVNTRs.txt ; done

#2#Then, SVs were classified as NAHR if both breakpoints were annotated as the same repeat class and the sequence identity between the two breakpoints is more than 80%. 
#NAHR
#In the first step, the flanking range of 200bp above and below the breakpoints on the left and right sides of the SV was taken
#100bp on the left flank, and the sequence was extracted

#bedtools flank -i TMP-SV-removeVNTRs.txt -g TMP-CHROM-leng.txt -l 180 -r 0 | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3+20,$4,$5,$6}' > TMP-LeftBreakpoint-SV-ragin.bed
bedtools flank -i TMP-SV-removeVNTRs.txt -g TMP-CHROM-leng.txt -l 100 -r 0 | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$4,$5,$6}' > TMP-LeftBreakpoint-SV-ragin.bed

bedtools getfasta -fi ${ref} -bed TMP-LeftBreakpoint-SV-ragin.bed -name -fo TMP-LeftBreakpoint-SV-ragin.fa
#100bp on the right flank, and the sequence was extracted

#bedtools flank -i TMP-SV-removeVNTRs.txt -g TMP-CHROM-leng.txt -l 0 -r 180 | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2-20,$3,$4,$5,$6}' > TMP-RightBreakpoint-SV-ragin.bed 
bedtools flank -i TMP-SV-removeVNTRs.txt -g TMP-CHROM-leng.txt -l 0 -r 100 | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$4,$5,$6}' > TMP-RightBreakpoint-SV-ragin.bed

bedtools getfasta -fi ${ref} -bed TMP-RightBreakpoint-SV-ragin.bed -name -fo TMP-RightBreakpoint-SV-ragin.fa

#3.#Blast all vs all
cat TMP-LeftBreakpoint-SV-ragin.fa TMP-RightBreakpoint-SV-ragin.fa > TMP-Bothbreakpoint.fa
#makeblastdb
makeblastdb -in TMP-Bothbreakpoint.fa -dbtype nucl -out TMPBothbreakpoint
#Blastn
blastn -db TMPBothbreakpoint -query TMP-Bothbreakpoint.fa -outfmt 6 -out TMP.all-vs-all.tsv -num_threads 100

paste TMP-LeftBreakpoint-SV-ragin.bed TMP-RightBreakpoint-SV-ragin.bed | awk 'BEGIN{OFS="\t";FS="\t"}{print $4"::"$1":"$2"-"$3,$4"::"$7":"$8"-"$9}' > TMP-two-breakpoint.txt
#NAHR
awk 'NR==FNR{a[$1,$2];next} ($1,$2) in a' TMP-two-breakpoint.txt TMP.all-vs-all.tsv | sed 's/::/\t/g' | awk 'BEGIN{OFS="\t";FS="\t"}{if($5>75){print $0}}' | uniq > NAHRs.result.txt
awk 'NR==FNR{a[$1,$2];next} ($1,$2) in a' TMP-two-breakpoint.txt TMP.all-vs-all.tsv | sed 's/::/\t/g' | awk 'BEGIN{OFS="\t";FS="\t"}{if($5>75){print $1}}' | uniq > NAHRs.final.txt

#Finally, the SVs overlapped with interspersed repetitive sequences were inferred into TE-mediated mechanisms.
#TE
#Removal NAHR
cat NAHRs.txt | while read y ; do sed -i "/$y/d" TMP-SV-removeVNTRs.txt ; done

mv TMP-SV-removeVNTRs.txt TMP-SV-removeVNTRs.NAHRs.txt


#All overlaps with TE
bedtools intersect -wa -wb -a TMP-SV-removeVNTRs.NAHRs.txt -b HB_TE.bed > TMP.ALL-TEinsertions.txt
#All sites that do not overlap with TE
bedtools intersect -wa -wb -a TMP-SV-removeVNTRs.NAHRs.txt -b HB_TE.bed | awk 'BEGIN{OFS="\t";FS="\t"}{print $4}' | sort -n | uniq > TMP.ALL.TEname.txt
bedtools intersect -wa -wb -a TMP-SV-removeVNTRs.NAHRs.txt -b HB_TE.bed | awk 'BEGIN{OFS="\t";FS="\t"}{print $0}' | sort -n | uniq > TE.result.txt
#This is the set of SVS associated with TE
#Determine single TE insertion（STEIs)
bedtools groupby -i TMP.ALL-TEinsertions.txt -g 1-4 -c 9 -o count | awk 'BEGIN{OFS="\t";FS="\t"}{if($5==1) {print $4}}' > STEIs.final.txt
#Determine multiple TE insertions（MTEIs)
bedtools groupby -i TMP.ALL-TEinsertions.txt -g 1-4 -c 9 -o count | awk 'BEGIN{OFS="\t";FS="\t"}{if($5!=1) {print $4}}' > MTEIs.final.txt

cat TMP.ALL.TEname.txt | while read y ; do sed -i "/$y/d" TMP-SV-removeVNTRs.NAHRs.txt ; done

cat TMP-SV-removeVNTRs.NAHRs.txt | awk 'BEGIN{OFS="\t";FS="\t"}{print $4}' > NHR.final.txt

rm TMP*
