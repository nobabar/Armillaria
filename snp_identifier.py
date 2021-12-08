import pandas as pd
import json
import glob
import os


files_dir = r"C:\Users\bapt0\Desktop\Stage_cyril\ITSarmillaire\test"
script_dir = r'C:\Users\bapt0\Desktop\Seq_identifier'


def load_to_dict(file):
    with open(file, "r") as j:
        x = json.loads(j.read())
    return x


def identifier_intra(sequence, snp_dict, group):
    id_dict = {}
    for specie in group:
        id_dict[specie] = {}
        id_dict[specie]["score"] = 0
        id_dict[specie]["N"] = 0
        for snp in snp_dict[specie].keys():
            nucl = sequence[snp-1]
            if nucl in snp_dict[specie][snp]:
                id_dict[specie]["score"] = id_dict[specie]["score"] + 1
            if nucl == "N":
                id_dict[specie]["N"] = id_dict[specie]["N"] + 1
    # print(f"{id_dict=} - {id_dict.keys()=}: {id_dict.values()=}\n\n")
    return id_dict


def identifier_inter(sequence, snp_dict, group):
    id_dict = {}
    id_dict["score"] = 0
    id_dict["N"] = 0
    for snp in snp_dict[group].keys():
        nucl = sequence[snp-1]
        if nucl in snp_dict[group][snp]:
            id_dict["score"] = id_dict["score"] + 1
        if nucl == "N":
            id_dict["N"] = id_dict["N"] + 1
    # print(f"{id_dict=} - {id_dict.keys()=}: {id_dict.values()=}\n\n")
    return id_dict


def group_branch(sequence, snp_dict, phylo):
    specie_ind = {}
    out_member = False
    if "out_group" in phylo.keys():
        out_group_specie = identifier_intra(sequence, snp_dict, phylo["out_group"])
        for specie in out_group_specie.keys():
            specie_ind = {specie: [out_group_specie[specie]["score"] , out_group_specie[specie]["N"]]}
            
            
            if out_group_specie[specie]["score"] >= int(0.8*len(snp_dict[specie])-out_group_specie[specie]["N"]):
                specie_ind["specie"] = specie
                out_member = True
    for groups in phylo["groups"]:
        group_inter = identifier_inter(sequence, snp_dict, groups)
        specie_ind["_".join(phylo["groups"][groups])] = [group_inter["score"] , group_inter["N"]]
        
        
        group_specie = identifier_intra(sequence, snp_dict, group["groups"][groups])
        for specie in group_specie.keys():
            specie_ind[specie] = [group_specie[specie]["score"] , group_specie[specie]["N"]]
            
            
            if out_member == False and group_inter["score"] >= len(snp_dict[groups])-1 and group_specie[specie]["score"] >= len(snp_dict[specie])-1:
                if "specie" in specie_ind.keys():
                    specie_ind["specie"] = specie_ind["specie"].append(specie)
                else:
                    specie_ind["specie"] = specie
    if "specie" not in specie_ind.keys():
        specie_ind["specie"] = "undefined"
    return specie_ind

group = load_to_dict(script_dir+'\species.JSON')
species = group.pop("species")
# locals().update(species)

snps_groups = pd.read_csv(script_dir + '\snp\snps.csv')
snps_out_group = pd.read_csv(script_dir + '\snp\snps_out_group.csv')

snps = {}
for key, value in group.items():
    if key == "out_group":
        for specie in value:
            snp = snps_out_group[snps_out_group["species"]==specie]
            snp = snp.filter(like=specie, axis=1) 
            snp = snp.rename(columns = lambda x: x.rsplit('_', 1)[1])
            snps[specie] = {}
            for col in snp:
                snps[specie][int(col)] = set((snp[col]))
    elif key == "groups":
        for sub_group in value:
            snp = snps_groups[snps_groups["species"].isin(value[sub_group])]
            snp = snp.filter(like="inter", axis=1)
            snp = snp.filter(like=value[sub_group][0], axis=1)
            snp = snp.rename(columns = lambda x: x.rsplit('_', 1)[1])
            snps[sub_group] = {}
            for col in snp:
                snps[sub_group][int(col)] = set((snp[col]))
            for specie in value[sub_group]:
                snp = snps_groups[snps_groups["species"]==specie]
                snp = snp.filter(like=specie, axis=1)
                snp = snp.filter(like="intra", axis=1)
                snp = snp.rename(columns = lambda x: x.rsplit('_', 1)[1])
                snps[specie] = {}
                for col in snp:
                    snps[specie][int(col)] = set((snp[col]))

files_directory = []
files_names = []
for file_name in glob.glob(files_dir+"/**/*.fas", recursive=True):
    files_directory.append(file_name)
    head, tail = os.path.split(file_name)
    files_names.append(tail.split('.', 1)[0])

for file in files_directory:
    # the directory to the alingned file
    aligned_dir = file.replace(".fas", "_aligned_cleaned.fasta")
    JSON = file.replace(".fas", "_species.JSON")
    csv = file.replace(".fas", "_species.csv")
    
    # copy the content of the aligned file
    with open(aligned_dir, "r") as align:
        cont = align.read().split('>')[1:]
    
    names, seq = [], []
    for x in cont:
        names.append(x.split('\n', 1)[0])
        seq.append(x.split('\n', 1)[1].replace('\n', ''))
    
    specie = {}
    for i in range(len(names)):
        specie[names[i]] = group_branch(seq[i], snps, group)
    
    with open(JSON, 'w') as result:
        result.write(json.dumps(specie, indent=2))
