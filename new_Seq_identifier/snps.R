library(dplyr)
library(writexl)


dir = 'C:/Users/bapt0/Desktop/Seq_identifier'
data = read.csv(paste(dir, '/consensus_seq.csv', sep=''))
species=c("mellea", "ostoyae", "cepistipes", "borealis", "gallica", "tabescens")
data=filter(data, (data$species %in% c("mellea", "tabescens"))==FALSE)

filter_snp_inter=function(specie, specie_set){
  snp_species=data.frame(table(data$X))[1]
  for (i in seq(3,length(specie_set))){
    if (dim(table(specie_set[i]))==1){
      out_specie_set = filter(data, (data$species %in% specie)==FALSE)[names(specie_set)[i]]
      if ((names(table(specie_set[i])) %in% names(table(out_specie_set)))==FALSE){
        snp_species[,names(specie_set[i])]=data[i]
      }
    }
  }
  return(snp_species)
}


filter_snp_intra=function(specie, specie_comb){
  snp_species=data.frame(specie_comb$X)[1]
  for (i in seq(3,length(specie_comb))){
    if (dim(table(specie_comb[i]))>1){
      for (j in seq(length(specie))){
        specie_set = filter(specie_comb, specie_comb$species == specie[j])[i]
        out_specie_set = filter(specie_comb, specie_comb$species == specie[-j])[i]
        if (length(table(intersect(specie_set,out_specie_set)))==0){
        # if ((names(table(specie_set)) %in% names(table(out_specie_set)))==FALSE){
          snp_species[,names(specie_set)]=specie_comb[i]
        }
      }
    }
  }
  return(snp_species)
}


ind_bore_osto = filter(data, data$species %in% c("borealis", "ostoyae"))
ind_cepi_gall = filter(data, data$species %in% c("cepistipes", "gallica"))

snp_inter_bore_osto = filter_snp_inter(c("borealis", "ostoyae"), ind_bore_osto)
  write.csv(snp_inter_bore_osto, paste(dir, '/snp/snp_inter_bore_osto.csv', sep = ''))
snp_inter_cepi_gall = filter_snp_inter(c("cepistipes", "gallica"), ind_cepi_gall)
  write.csv(snp_inter_cepi_gall, paste(dir, '/snp/snp_inter_cepi_gall.csv', sep = ''))

snp_intra_bore_osto = filter_snp_intra(c("borealis", "ostoyae"), ind_bore_osto)
  write.csv(snp_intra_bore_osto, paste(dir, '/snp/snp_intra_bore_osto.csv', sep = ''))
snp_intra_cepi_gall = filter_snp_intra(c("cepistipes", "gallica"), ind_cepi_gall)
  write.csv(snp_intra_cepi_gall, paste(dir, '/snp/snp_intra_cepi_gall.csv', sep = ''))


# write.csv(snp_inter_bore_osto, paste(dir, '/snp/snp_inter_bore_osto.csv', sep = ''))
# write_xlsx(cbind(snp_mellea, snp_ostoyae, snp_cepistipes, snp_borealis, snp_gallica, snp_tabescens), paste(dir, '/snp/snp_species.xlsx', sep = ''))
