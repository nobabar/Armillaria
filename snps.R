library(dplyr)
library(writexl)
library(rjson)

# dir <- commandArgs
# print(dir)
# dir <- 'C:/Users/bapt0/Desktop/Seq_identifier'

# group_1 <- c("borealis", "ostoyae")
# group_2 <- c("cepistipes", "gallica")
# out_group <- c("tabescens", "mellea")
# species <- names(table(cbind(group_1, group_2, out_group)))

# group <- list(species = species,
#                groups = list(group_1 = group_1,
#                              group_2 = group_2), 
#                out_group = out_group)
# groups_json <- toJSON(group)
# write(groups_json, 'C:/Users/bapt0/Desktop/Seq_identifier/species.JSON')

# groups <- fromJSON(file='C:/Users/bapt0/Desktop/Seq_identifier/species.JSON')
# 
# data_all <- read.csv(paste(dir, '/consensus_seq.csv', sep=''))
# data_groups <- filter(data_all, (data_all$species %in% groups$out_group)==FALSE)
# data_out <- filter(data_all, data_all$species %in% out_group)


filter_snp_inter=function(specie, specie_set, data_set){
  snp_species <- data_set[2]
  for (i in seq(3,length(specie_set))){
    if (dim(table(specie_set[i]))==1){
      out_specie_set <- filter(data_set, (data_set$species %in% specie)==FALSE)[names(specie_set)[i]]
      if ((names(table(specie_set[i])) %in% names(table(out_specie_set)))==FALSE){
        snp_species[,names(specie_set[i])] <- data_set[i]
      }
    }
  }
  return(snp_species)
}


filter_snp_intra=function(specie, specie_comb, data_set){
  snp_species <- data_set[2]
  for (i in seq(3,length(specie_comb))){
    if (dim(table(specie_comb[i]))>1){
      for (j in seq(length(specie))){
        specie_set <- filter(specie_comb, specie_comb$species == specie[j])[i]
        out_specie_set <- filter(specie_comb, specie_comb$species == specie[-j])[i]
        if (length(table(intersect(specie_set,out_specie_set)))==0){
        # if ((names(table(specie_set)) %in% names(table(out_specie_set)))==FALSE){
          snp_species[,names(specie_set)] <- data_set[i]
        }
      }
    }
  }
  return(snp_species)
}

by_groups=function(group, data_set){
  if (deparse(substitute(group)) == "out_group"){
    snp_group <- data_set[,c(1, 2)]
    for (i in seq(length(group))){
      ind_group <- filter(data_set, data_set$species == group[i])
      
      snp_group_i <- filter_snp_inter(group[i], ind_group, data_set)
      write.csv(snp_group_i, paste(dir, '/snp/snp_', group[i], '.csv', sep = ''))
      colnames(snp_group_i) <- paste(group[i], gsub('X', '', colnames(snp_group_i)), sep = "_")
      
      snp_group <- cbind(snp_group, snp_group_i[,-1])
    }
  }else{
    ind_group <- filter(data_set, data_set$species %in% group)
    
    snp_group_inter <- filter_snp_inter(group, ind_group, data_set)
    write.csv(snp_group_inter, paste(dir, '/snp/snp_inter', paste(group, collapse='_'), '.csv', sep = ''))
    colnames(snp_group_inter) <- paste(paste(group, collapse=''), 'inter', gsub('X', '', colnames(snp_group_inter)), sep = "_")
    
    snp_group_intra <- filter_snp_intra(group, ind_group, data_set)
    write.csv(snp_group_intra, paste(dir, '/snp/snp_intra', paste(group, collapse='_'), '.csv', sep = ''))
    colnames(snp_group_intra) <- paste(paste(group, collapse=''), '_intra_', gsub('X', '', colnames(snp_group_intra)), sep = "_")
    
    snp_group <- cbind(snp_group_inter, snp_group_intra[,-1])
  }
  return (snp_group)
}

in_out=function(dir){
  groups <- fromJSON(file='C:/Users/bapt0/Desktop/Seq_identifier/species.JSON')
  
  data_all <- read.csv(paste(dir, '/consensus_seq.csv', sep=''))
  data_groups <- filter(data_all, (data_all$species %in% groups$out_group)==FALSE)
  
  snps <- data_groups[,c(1,2)]
  for (group in names(groups)){
    if (group == "out_group"){
      snps_out_group <- by_groups(groups$out_group, data_all)
    }else{
      snp_group <- by_groups(groups$group, data_groups)
      snps <- cbind(snps, snp_group[,-1])
    }
  }
  write.csv(snps, paste(dir, '/snp/snps.csv', sep = ''))
  write_xlsx(snps, paste(dir, '/snp/snps.xlsx', sep = ''))
  write.csv(snps_out_group, paste(dir, '/snp/snps_out_group.csv', sep = ''))
  write_xlsx(snps_out_group, paste(dir, '/snp/snps_out_group.xlsx', sep = ''))
}
  
# snp_group_1 <- groups(group_1, data_groups)
# snp_group_2 <- groups(group_2, data_groups)
# snp_group <- data.frame(table(snp_group))
# snp_group <- filter(snp_group, snp_group$Freq != 0)[,-length(snp_group)]
# snps <- cbind(data_groups[,c(1, 2)], snp_group_1[,-1], snp_group_2[,-1])

# write.csv(snps, paste(dir, '/snp/snps.csv', sep = ''))
# write_xlsx(snps, paste(dir, '/snp/snps.xlsx', sep = ''))
# write.csv(snps_out_group, paste(dir, '/snp/snps_out_group.csv', sep = ''))
# write_xlsx(snps_out_group, paste(dir, '/snp/snps_out_group.xlsx', sep = ''))

# snps_out_group <- filter(snps_out_group, snps_out_group$species %in% out_group)
# for (specie in out_group){
#   sub_snp_out <- filter(snps_out_group, snps_out_group$species == specie)
#   sub_snp_out <- cbind(snps_out_group[2], select(snps_out_group, contains(specie)))
#   sub_out_group <- data.frame(table(sub_snp_out))
#   sub_out_group <- filter(sub_out_group, sub_out_group$Freq != 0)
# }
