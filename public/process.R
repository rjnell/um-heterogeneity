# Script to anonymise data

# LUMC-data
lumc_data = read_xlsx("data/Supplementary Tables.xlsx", sheet=1, skip=1)

# Load LUMC data
lumc_normalised = as.matrix(read.csv("data/normalised-80.tsv",sep="\t", stringsAsFactors = F, check.names=F, row.names = 1))
colnames(lumc_normalised) = lumc_data$LUMC_ID[match(colnames(lumc_normalised), lumc_data$Tumor_ID)]

write.table(lumc_normalised, "public/gene-expression-data", sep="\t")


# Purity
purity = read_xlsx("data/Supplementary Tables.xlsx", sheet=2)
ID = lumc_data$LUMC_ID[match(purity$Tumor_ID, lumc_data$Tumor_ID)]
ID[which(nchar(purity$Tumor_ID) == 12)] = "LUMC-34 (2nd)"
write.table(cbind(ID, purity[,3:7]), "public/purity", sep="\t")

View(cbind(ID, purity[,3:7]))


# Add
add = read_xlsx("data/Supplementary Tables.xlsx", sheet=3)
ID = lumc_data$LUMC_ID[match(add$Tumor_ID, lumc_data$Tumor_ID)]
write.table(cbind(ID, add[,2:6]), "public/add", sep="\t")

View(cbind(ID, add[,2:6]))


# Add
chr = read_xlsx("data/Supplementary Tables.xlsx", sheet=5)
ID = lumc_data$LUMC_ID[match(chr$V1...1, lumc_data$Tumor_ID)]
res = cbind(ID, chr[2:6])
colnames(res) = c("ID", "Approach", "Normalised result", "99%-CI_low", "99%-CI_high", "Conclusion")
write.table(res, "public/chr3p", sep="\t")

res = cbind(ID, chr[7:11])
colnames(res) = c("ID", "Approach", "Normalised result", "99%-CI_low", "99%-CI_high", "Conclusion")
write.table(res, "public/chr8q", sep="\t")

