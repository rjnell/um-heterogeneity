# Load data
lumc_data = readxl::read_xlsx("data/Supplementary Tables.xlsx", sheet=1, skip=1)
ids = lumc_data$Tumor_ID
sample = ids[11]

output_path = "/mnt/f/mutations/res/"
output_path_win = "F:/mutations/res/"

# Create dirs
system(paste0("wsl mkdir ", output_path, sample))

system(paste0('wsl freebayes -f "/mnt/f/38.fa" "/mnt/f/LUMC/',sample,'.bam" > "/mnt/f/mutations/res/',sample,'/',sample, '.vcf"'))

SS_SR_VCF_path = paste0(output_path, sample, "/", sample, ".vcf")
SS_SR_VCF_path_win = paste0(output_path_win, sample, "/", sample, ".vcf")
SS_SR_AV_path = paste0(output_path, sample, "/", sample,  ".av")
SS_SR_CSV_path = paste0(output_path, sample, "/", sample)
SS_SR_CSV_path_win = paste0(output_path_win, sample, "/", sample)
SS_SR_OUTPUT_path_win = paste0(output_path_win, sample, "/", sample, ".csv")

# Convert2annovar
annovar1 = paste0("wsl perl /mnt/f/mutations/bin/annovar/convert2annovar.pl -format vcf4 ",
                  SS_SR_VCF_path,
                  " > ",
                  SS_SR_AV_path)
system(annovar1)

# Table_annovar
annovar2 = paste0("wsl perl /mnt/f/mutations/bin/annovar/table_annovar.pl ",
                  SS_SR_AV_path,
                  " /mnt/f/mutations/bin/annovar/humandb -buildver hg38 -remove -protocol refGene -operation gx -nastring . -csvout -polish -xref /mnt/f/mutations/bin/annovar/example/gene_xref.txt",
                  " -out ", SS_SR_CSV_path)
system(annovar2)

# Open vcf_file for gene and sample
vcf_file = vcfR::read.vcfR(SS_SR_VCF_path_win)      

cols = c("GT","ABQ","AD","ADF","ADR","DP","FREQ","GQ","PVAL","RBQ","RD","RDF","RDR","SDP","RO","QR","AO","QA","GL")
vcf_table = matrix(nrow = nrow(vcf_file@fix), ncol = 4+length(cols))
vcf_table[,1:4] = vcf_file@fix[,c(1,2,4,5)]
colnames(vcf_table) = c("CHR","POS","REF","ALT",cols)
for (index in row(vcf_file@fix)[,1]) {
  formats = strsplit(vcf_file@gt[index, 1], split=":")[[1]]  
  values = strsplit(vcf_file@gt[index, 2], split=":")[[1]]  
  vcf_table[index,formats] = values
}

# Open annovar_file
annovar_file = read.csv(paste0(SS_SR_CSV_path_win, ".hg38_multianno.csv"), sep=",", stringsAsFactors = F)
as = NULL
for (i in 1:nrow(annovar_file)) {
  as = c(as,paste(annovar_file[i,1:2], collapse = "-"))
}

res = NULL
for (i in 1:nrow(vcf_table)) {
  w = which(as == paste(vcf_table[i,1:2], collapse = "-"))
  if (length(w) == 1) {
    res = rbind(res, c(vcf_table[i,],annovar_file[w,]))  
  }
  else {
    res = rbind(res, c(vcf_table[i,],rep(NA,ncol(annovar_file))))  
  }
  
}

write.table(res[,c(1,2,3,4,7,30,31,32,33)],paste0(SS_SR_CSV_path_win, "-res.tsv"), sep="\t")
v = read.csv(paste0(SS_SR_CSV_path_win, "-res.tsv"), sep="\t", stringsAsFactors = F)
View(v)
