nano ~/.bashrc
/home/user/tools/affymetrix_power_tools/bin/
source ~/.bashrc

cd /home/user/tools/affymetrix_power_tools/bin/
chmod +x *

apt-format-result \
  --calls-file ah.txt \
  --annotation-file Axiom_PMDA.na36.r7.a8.annot.db \
  --export-vcf-file output.vcf

Problem: probably need to remove annotation columns from the file.
