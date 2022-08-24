#!/usr/bin/env bash
#perform variant calling on the vegfa gene chromatogram files
for sample in `ls *.ab1`
do
        tracy decompose -v -a homo_sapiens_hg19 -r Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz $sample -o $sample
done

# convert to vcf file

for i in `ls *bcftools`
do
	output=`echo $i | cut -f1 -d"."`
	bcftools view $i -o ${output}.vcf
done

# screen variants

for i in `ls *vcf`
do
out_file_prefix=`basename $i .vcf`
vcftools --vcf $i --snps variant_ids_list.txt --recode --out $out_file_prefix
done


# extract information of interest 

echo -e "LOCATION\tID\tALLELE_VARIATION" > threes.txt

for i in `ls 3*bcf`
do
        #bcftools view $i | grep -v "#" | awk '{print $2"\t"$3"\t"$4"/"$5}' >> threes.txt
        bcftools view -v snps $i | bcftools query -f '%POS\t%ID\t%REF/%ALT\n' >> threes.txt
done
#bcftools view 3002.bcf | grep -v "#" | awk '{print $2"\t"$3"\t"$4"/"$5}'

echo -e "LOCATION\tID\tALLELE_VARIATION" > twos.txt

for i in `ls 2*bcf`
do
       bcftools view -v snps $i | bcftools query -f '%POS\t%ID\t%REF/%ALT\n' >> twos.txt
        # bcftools view $i | grep -v "#" | awk '{print $2"\t"$3"\t"$4"/"$5}' >> twos.txt
done
