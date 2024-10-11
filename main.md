```bash
nano ~/.bashrc
/home/user/tools/affymetrix_power_tools/bin/
source ~/.bashrc

cd /home/user/tools/affymetrix_power_tools/bin/
chmod +x *

apt-format-result \
  --calls-file test_batch.txt \
  --annotation-file Axiom_PMDA.na36.r7.a8.annot.db \
  --export-vcf-file output.vcf

cat ah.txt | cut -f 1-806 > test_batch.txt
# Problem: probably need to remove annotation columns from the file.

awk 'NR >= 7 {for (i = 2; i <= 805; i++) print $i}' ah.txt | sort | uniq -c
```

8849114 A
207764927 AA
90911653 AB
   6406 AC
     20 AD
      1 AG
      2 AJ
11981132 B
299500857 BB
   8396 BC
    161 BD
      2 BE
    671 C
 202298 CC
     73 CD
      2 CE
    434 D
   4237 DD
15344830 NoCall
 527062 NoCall_1
 153398 ZeroCN

 ```bash
awk 'NR < 7 {print; next} {for (i = 2; i <= 805; i++) if ($i != "AA" && $i != "AB" && $i != "BB" && $i != "NoCall" && $i != "NoCall_1") next} 1' ah.txt > ah1.txt

# changing A and B to 0 and 1:
awk 'NR >= 7 {
    genotypes = ""
    # Process columns 2 to 805
    for (i = 2; i <= 805; i++) {
        call = $i
        if (call == "AA") {
            gt = "0/0"
        } else if (call == "AB") {
            gt = "0/1"
        } else if (call == "BB") {
            gt = "1/1"
        } else if (call == "NoCall" || call == "NoCall_1") {
            gt = "./."
        } else {
            gt = "./."  # Convert any other values to missing genotype "./."
        }
        genotypes = genotypes gt "\t"
    }

    # Print the first column (unchanged), followed by the converted genotype string, and the rest of the columns
    print $1 "\t" genotypes substr($0, index($0, $806))  # $806 and beyond (rest of the row) unchanged
}' ah1.txt > ah2.txt

cat ah1.txt | head -n 6 | tail -n 1 > samples_header.tsv

# counting empty cells in 811 dbSNP
awk -F'\t' 'NF < 811 || $811 == "" {count++} END {print count}' ah2.txt
# 49171

awk -F'\t' 'NF < 836 || $836 == "" {count++} END {print count}' ah2.txt
# 7348
# Let's use extended rsID but first let'd delete all rows without rsIDs
awk -F'\t' '$836 != ""' ah2.txt > ah3.txt

# changing it into a vcf format:
paste <(awk -F'\t' '{print $807, $808, $836, $819, $820, ".", ".", ".", "GT"}' OFS='\t' ah3.txt) <(cut -f 2-805 ah3.txt) > ah4.txt
awk 'BEGIN {ORS=""; print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"} {ORS="\n"; print}' sample_ids.tsv > header_vcf.txt
cat header_vcf.txt ah4.txt > ah.vcf
plink --vcf ah.vcf --make-bed --out ah

# NO WAY! IT ACTUALLY WORKED! Let's run add sex and phenotypes and run the analysis
plink --bfile ah --update-sex "head_ah.tsv - sex.tsv"  --pheno "head_ah.tsv - pheno.tsv" --make-bed --out ah1
awk '$6 == -9' ah1.fam | awk '{print $1"\t" $2}' > missing_phenotype.tsv
plink --bfile ah1 --remove missing_phenotype.tsv --allow-no-sex --make-bed --out ah2
plink --bfile ah2 --geno 0.02 --make-bed --out ah3
plink --bfile ah3 --mind 0.02 --make-bed --out ah4
plink --bfile ah4 --maf 0.001 --make-bed --out ah5
plink --bfile ah5 --genome --min 0.2 --out pihat_min0.2
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
plink2 --bfile ah5 --pca 10 --out ah_pca

plink --bfile ah5 --covar ah_pca.eigenvec --logistic --hide-covar --thread-num 8 --out simple_logistic
plink --bfile ah5 --covar ah_pca.eigenvec --logistic --dominant --hide-covar --out dominant_results
plink --bfile ah5 --covar ah_pca.eigenvec --logistic --recessive --hide-covar --out recessive_results
plink --bfile ah5 --assoc --out assoc_results
plink --bfile ah5 --fisher
plink --bfile ah5 --model
plink2 --bfile ah5 --glm --covar ah_pca.eigenvec

cat dominant_results.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat recessive_results.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat simple_logistic.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat assoc_results.assoc | awk '$9 != "NA"' | sort -gk 9,9 | head
cat plink.assoc.fisher | awk '$8 != "NA"' | sort -gk 8,8 | head
cat plink.model | awk '$10 != "NA"' | sort -gk 10,10 | head
cat plink2.PHENO1.glm.logistic.hybrid | awk '$15 != "NA"' | sort -gk 15,15 | head
```
