import subprocess as sp
import json

from consensus_sequence import sites

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Directories
general_dir = 'C:/Users/bapt0/Desktop/Seq_identifier'
ref_file_dir = 'C:/Users/bapt0/Desktop/Seq_identifier/ref_seq.fas'

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Groups and species
groups = {}
species = []
out = input("Is there an out group of species? [y/n]: ")
if out == "y":
    specie_list = input("----Enter the list of species for your out group (seperated by a simple comma): ").replace(" ", "").split(",")
    groups["out_group"] = specie_list
    species = specie_list
i = 1
groups["groups"]={}
while(True):
    group = input(f"Enter the list of species for group {i} (enter s to stop): ")
    if group == "s":
        break
    groups["groups"][f"group_{i}"] = group.replace(" ", "").split(",")
    species = species + group.replace(" ", "").split(",")
    i += 1

with open(general_dir+'/species.JSON', 'w+') as JSON:
    JSON.write(json.dumps(groups))

# out_group = ["tabescens", "mellea"]
# group_1 = ["borealis", "ostoyae"]
# group_2 = ["cepistipes", "gallica"]
# groups={out_group: out_group, group_1: group_1, group_2: group_2}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Cleaning the set of reference sequence and
# Determination of possible snp sites
while(True):
    excel = input("Would you like an extra excel file as output? [y/n]: ")
    if excel != "y" and excel != "n":
        print("-----Wrong input!-----")
    else:
        break
snp_sites = sites(ref_file_dir, excel)
if excel == "y":
    print("Polymorphic sites localized. Csv and Excel files created in the directory of the reference file")
else:
    print("Polymorphic sites localized. Csv file created in the same directory as the reference file")

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Call R script
# Determine the snps from the isolated sites for each specie
# cmd = ["in_out", general_dir+'/new_Seq_identifier/snps.R'] + [general_dir]
cmd = ["in_out", 'C:/Users/bapt0/Desktop/Seq_identifier/new_Seq_identifier/snps.R'] + [general_dir, excel]
sp.check_output(cmd)



