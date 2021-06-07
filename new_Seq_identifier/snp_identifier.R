library(dplyr)
library(writexl)


dir = 'C:/Users/bapt0/Desktop/Seq_identifier'
data = read.csv(paste(dir, '/consensus_seq.csv', sep=''))
species=c("mellea", "ostoyae", "cepistipes", "borealis", "gallica", "tabescens")


filter_snp=function(specie, specie_set){
  snp_species=setNames(data.frame(table(data$species))[1], paste('snp_', specie, sep=''))
  for (i in seq(3,length(specie_set))){
    if (dim(table(specie_set[i]))==1){
      tag=0
      for (out_specie in species){
        if ((specie==out_specie)==FALSE){
          if(names(table(specie_set[i])) %in% names(table(filter(data, data$species==out_specie)[names(specie_set)[i]]))){
            tag=tag+1
          }
        }
      }
      if (tag < 2){
        snp_species[,names(specie_set)[i]]=rbind(paste(names(table(ind_borealis[names(specie_set)[i]])), collapse=""),
                                                 paste(names(table(ind_cepistipes[names(specie_set)[i]])), collapse=""),
                                                 paste(names(table(ind_gallica[names(specie_set)[i]])), collapse=""),
                                                 paste(names(table(ind_mellea[names(specie_set)[i]])), collapse=""),
                                                 paste(names(table(ind_ostoyae[names(specie_set)[i]])), collapse=""),
                                                 paste(names(table(ind_tabescens[names(specie_set)[i]])), collapse=""))
      }
    }
  }
  return(snp_species)
}


for (specie in species){
  eval(parse(text = paste('ind_', specie, ' = filter(data, data$species==specie)', sep = "")))
}


snp_mellea = filter_snp(specie = "mellea", specie_set=ind_mellea)
  write.csv(snp_mellea, paste(dir, '/snp_mellea.csv', sep = ''))
snp_ostoyae = filter_snp(specie = "ostoyae", specie_set=ind_ostoyae)
  write.csv(snp_ostoyae, paste(dir, '/snp_ostoyae.csv', sep = ''))
snp_cepistipes = filter_snp(specie = "cepistipes", specie_set=ind_cepistipes)
  write.csv(snp_cepistipes, paste(dir, '/snp_cepistipes.csv', sep = ''))
snp_borealis = filter_snp(specie = "borealis", specie_set=ind_borealis)
  write.csv(snp_borealis, paste(dir, '/snp_borealis.csv', sep = ''))
snp_gallica = filter_snp(specie = "gallica", specie_set=ind_gallica)
  write.csv(snp_gallica, paste(dir, '/snp_gallica.csv', sep = ''))
snp_tabescens = filter_snp(specie = "tabescens", specie_set=ind_tabescens)
  write.csv(snp_tabescens, paste(dir, '/snp_tabescens.csv', sep = ''))

# write_xlsx(cbind(snp_mellea, snp_ostoyae, snp_cepistipes, snp_borealis, snp_gallica, snp_tabescens), paste(dir, '/snp_species.xlsx', sep = ''))
