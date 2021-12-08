from Bio import SeqIO
import termcolor
import glob
import os


general_dir = r"C:\Users\bapt0\Desktop\Stage_cyril\FOL662\FOL66269_SEQ_BEC87"
folders = []
for folder in os.walk(general_dir):
    if "ab1" in folder[0]:
        folders.append(folder[0])

for folder in folders:
    files_directory, files_names = [], []
    for file_name in glob.iglob(folder+"/**/*.ab1", recursive=True):
        files_directory.append(file_name)
        head, tail = os.path.split(file_name)
        files_names.append(tail.split('.', 1)[0])

    seqName, seqStr, seqQual = [], [], []
    for files in files_directory:
        record = SeqIO.parse(files, "abi")
        # count = SeqIO.write(record, files.replace(".ab1", ".fastq"), "fastq")
        for seq in record:
            seqName.append(seq.id)
            seqStr.append(str(seq.seq))
            seqQual.append(seq.letter_annotations["phred_quality"])

    i = 0
    while i < len(seqQual):
        for j in range(len(seqQual[i])):
            if seqQual[i][j] < 20:
                if j == 0:
                    seqStr[i] = seqStr[i][1:]
                else:
                    seqStr[i] = seqStr[i][:j-1]+'N'+seqStr[i][j-1:]
        if seqStr[i].count("N") > 0.15*len(seqStr[i]):
            seqName.pop(i)
            seqStr.pop(i)
            seqQual.pop(i)
        else:
            while seqStr[i][:5].count("N") > 4:
                seqStr[i] = seqStr[i][5:]
            while seqStr[i][-5:].count("N") > 4:
                seqStr[i] = seqStr[i][:-5]
            print(seqName[i] + ":", termcolor.colored(seqStr[i].count("N"), "cyan"),
                  "N for a length of", termcolor.colored(len(seqStr[i]), "cyan"),
                  "thus a ratio of", termcolor.colored(seqStr[i].count("N")/len(seqStr[i])*100, "cyan"))
            i += 1
    if len(seqName) > 0:
        open(folder.rsplit("\\", 1)[0]+"\\"+folder.rsplit("\\", 1)[1]+".fas", "w+").close()
        with open(folder.rsplit("\\", 1)[0]+"\\"+folder.rsplit("\\", 1)[1]+".fas", "a") as file:
            for i in range(len(seqName)):
                file.write(">"+seqName[i]+"\n"+seqStr[i]+"\n")
