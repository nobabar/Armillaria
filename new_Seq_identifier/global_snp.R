library(dplyr)


data = read.csv('C:/Users/bapt0/Desktop/Seq_identifier/consensus_seq.csv')
species=c("mellea", "ostoyae", "cepistipes", "borealis", "gallica", "tabescens")


for (specie in species){
  eval(parse(text = paste('ind_', specie, ' = filter(data, data$species==specie)', sep = "")))
}
snp_species=data.frame(table(data$species))


for (s in species){
  specie_set=eval(parse(text=paste('ind_', s, sep='')))
  i=3
  while (i <= length(specie_set)){
    if (dim(table(specie_set[i]))==1){
      tag=0
      for (out_specie in species){
        if ((s==out_specie)==FALSE){
          if(names(table(specie_set[i])) %in% names(table(filter(data, data$species==out_specie)[names(specie_set)[i]]))){
            tag=tag+1
          }
        }
      }
      if (tag < 2){
        snp_species[,names(specie_set)[i]] = rbind(paste(names(table(ind_mellea[names(specie_set)[i]])), collapse=""),
                                                   paste(names(table(ind_ostoyae[names(specie_set)[i]])), collapse=""),
                                                   paste(names(table(ind_cepistipes[names(specie_set)[i]])), collapse=""),
                                                   paste(names(table(ind_borealis[names(specie_set)[i]])), collapse=""),
                                                   paste(names(table(ind_gallica[names(specie_set)[i]])), collapse=""),
                                                   paste(names(table(ind_tabescens[names(specie_set)[i]])), collapse=""))
      }
    }
    i=i+1
  }
}
