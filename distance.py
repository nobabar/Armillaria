import DNA_class
import re

def create_dna(file):
    # name = []
    with open(file, "r") as seq:
        lines = seq.read().strip()
        nameseq = re.compile("(^>.+)\s((?:[ATCG-]+\s*)+)", re.MULTILINE)
        seqs = nameseq.findall(lines)
             # name.append(line[0])
             # content=line[1]
        print(seqs)
        for i in seqs:
            dna = DNA_class.Dna(i[0],i[1])
            print(dna.name, dna.seq, "\n\n")
