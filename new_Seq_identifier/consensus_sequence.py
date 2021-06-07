import pandas as pd
import os


ref_file_dir = r"C:\Users\bapt0\Desktop\Seq_identifier\ref_seq.fas"


with open(ref_file_dir, "r") as file:
    cont = file.read().split('>')[1:]


name, seq = [], []
for x in cont:
    name.append(x.split('\n', 1)[0])
    seq.append(x.split('\n', 1)[1].replace('\n', ''))

species = list(i.rsplit(".", 1)[1] for i in name)


loci = []
gaps = []
for i in range(len(seq[0])):
    locus = ''
    gp = ''
    for x in seq:
        locus += x[i]
        if x[i] == '-':
            gap = '-'
    loci.append(locus)
    gaps.append(gap)


# truncate start and end
while loci[0].count("-") > 0.2*len(cont):
    loci.pop(0)
while loci[-1].count("-") > 0.2*len(cont):
    loci.pop(-1)


snp = {}
seq_final = []
for i in range(len(loci)):
    if len(set(loci[i])) >= 2:
        snp_loc = []
        for j in range(len(loci[i])):
            snp_loc.append(list(loci[i])[j])
        tag = True
        for k in set(loci[i]):
            if loci[i].count(k) == 1:
                tag = False
        if len(set(loci[i])) == 2 and "-" in set(loci[i]):
            tag = False
        if tag:
            if i+1 in snp:
                snp[i+1].append(snp_loc)
            else:
                snp[i+1] = snp_loc
    seq_final.append(loci[i])


with open(ref_file_dir.replace(".fas", "_cleaned.fasta"), "w") as cln:
    for i in range(len(seq_final[0])):
        tag = ''
        for x in seq_final:
            tag = tag + x[i]
        cln.write('>' + name[i] + '\n' + tag + '\n')


d = {}
d["species"] = species
for i in snp.keys():
    a = list(snp[i])
    d[i] = a
df = pd.DataFrame(data=d, columns=d.keys(), index=name)
head = os.path.split(ref_file_dir)[0]
df.to_csv(head+'\consensus_seq.csv')
df.to_excel(head+'\consensus_seq.xlsx')
