from Bio import Entrez
import json


def load_to_dict(file):
    with open(file, "r") as j:
        x = json.loads(j.read())
    return x


general_dir = r'C:\Users\bapt0\Desktop\Stage_cyril\ITSarmillaire\GenBank'
script_dir = r'C:\Users\bapt0\Desktop\Seq_identifier'

group = load_to_dict(script_dir+'\species.JSON')
specie_list = group.pop("species")


def genbank(species, general_dir):
    Entrez.email = "bapt09.rousseau@gmail.com"
    Entrez.tool = "Armillaria_SNPidentifier"
    
    for specie in species:
        file = general_dir+r'\genbank_A'+specie+'.fasta'
        open(file, 'w+').close()
        
        handle_1 = Entrez.esearch(db="nucleotide",
                                  term="armillaria "+specie+" internal transcribed spacer",
                                  retmax=40,
                                  retmode="fasta")
        id_list = ",".join(Entrez.read(handle_1)["IdList"])
        handle_1.close()
        
        handle_2 = Entrez.epost(db = 'nucleotide', id = id_list)
        record = Entrez.read(handle_2)
        webEnv = record['WebEnv']
        queryKey = record['QueryKey']
        handle_2.close()
        
        handle_3 = Entrez.efetch(db = 'nucleotide', id = id_list,
                               rettype = 'fasta', retmod = 'text')
        with open(file, 'a') as file:
            for line in handle_3.readlines():
                file.write(line)
