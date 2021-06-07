library(dplyr)


dir <- '/home/arthur/Bureau/Armillaria-main'
data_all <- read.csv(paste(dir, '/consensus_seq.csv', sep=''))
species <- c("mellea", "ostoyae", "cepistipes", "borealis", "gallica", "tabescens")
data <- droplevels(data_all[-c(1:6,23:27),])

filter_snp_inter=function(specie, specie_set){
  snp_species <- data.frame(data[1])[1]
  for (i in seq(3,length(specie_set))){
    if (dim(table(droplevels(specie_set[i])))==1){
      out_specie_set = data[which((data$species %in% specie)==FALSE),][names(specie_set)[i]]
      if (length(table(intersect(names(table(droplevels(specie_set[i]))), names(table(droplevels(out_specie_set)))))) == 0){
        snp_species[,names(specie_set[i])] <- data[i]
      }
    }
  }
  return(snp_species)
}


filter_snp_intra=function(specie, specie_comb){
  snp_species=data.frame(data[1])[1]
  for (i in seq(3,length(specie_comb))){
    if (dim(table(droplevels(specie_comb[i])))>1){
      for (j in seq(length(specie))){
        specie_set <- filter(specie_comb, specie_comb$species == specie[j])[i]
        out_specie_set <- filter(specie_comb, specie_comb$species == specie[-j])[i]
        if (length(table(intersect(names(table(droplevels(specie_set))),names(table(droplevels(out_specie_set)))))) == 0){
          snp_species[,names(specie_set)] <- data[i]
        }
      }
    }
  }
  return(snp_species)
}


ind_bore_osto <- filter(data, data$species %in% c("borealis", "ostoyae"))
ind_cepi_gall <- filter(data, data$species %in% c("cepistipes", "gallica"))

snp_inter_bore_osto <- filter_snp_inter(c("borealis", "ostoyae"), ind_bore_osto)
  write.csv(snp_inter_bore_osto, paste(dir, '/snp/snp_inter_bore_osto.csv', sep = ''))
snp_inter_cepi_gall <- filter_snp_inter(c("cepistipes", "gallica"), ind_cepi_gall)
  write.csv(snp_inter_cepi_gall, paste(dir, '/snp/snp_inter_cepi_gall.csv', sep = ''))

snp_intra_bore_osto <- filter_snp_intra(c("borealis", "ostoyae"), ind_bore_osto)
  write.csv(snp_intra_bore_osto, paste(dir, '/snp/snp_intra_bore_osto.csv', sep = ''))
snp_intra_cepi_gall <- filter_snp_intra(c("cepistipes", "gallica"), ind_cepi_gall)
  write.csv(snp_intra_cepi_gall, paste(dir, '/snp/snp_intra_cepi_gall.csv', sep = ''))

colnames(snp_inter_bore_osto) <-paste("BO_vs_CG", gsub('X', '', colnames(snp_inter_bore_osto)), sep = "_")
colnames(snp_inter_cepi_gall) <-paste("CG_vs_BO", gsub('X', '', colnames(snp_inter_cepi_gall)), sep = "_")
colnames(snp_intra_bore_osto) <-paste("bore_osto", gsub('X', '', colnames(snp_intra_bore_osto)), sep = "_")
colnames(snp_intra_cepi_gall) <-paste("cepi_gall", gsub('X', '', colnames(snp_intra_cepi_gall)), sep = "_")

a <- cbind(data$X, snp_inter_bore_osto[,-1], snp_inter_cepi_gall[,-1], snp_intra_bore_osto[,-1], snp_intra_cepi_gall[,-1])

# write.csv(snp_inter_bore_osto, paste(dir, '/snp/snp_inter_bore_osto.csv', sep = ''))
# write_xlsx(cbind(snp_mellea, snp_ostoyae, snp_cepistipes, snp_borealis, snp_gallica, snp_tabescens), paste(dir, '/snp/snp_species.xlsx', sep = ''))
