# Holostei
Detecting purifying selection in Holostei vegfc genes

All this is just a script to prepare 10 fasta files for conservation analysis by comparing non-synonymous to synonymous codon changes in the mRNA. The final analysis is not done here but by HyPhy (https://github.com/veg/hyphy):

    hyphy slac --pvalue 0.05 --alignment holostei_vegfcd_mRNA.nxs --tree holostei_vegfcd.phy_phyml_tree.txt
