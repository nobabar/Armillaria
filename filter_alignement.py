import glob
import os


general_dir = r"C:\Users\bapt0\Desktop\Stage_cyril\ITSarmillaire\test"

# directory to the muscle command file
muscle_dir = r"C:\Users\bapt0\Desktop\Seq_identifier\muscle3.8.31_i86win32.exe"

# the directory to the file containing the reference sequences, must be .fas
ref_file_dir = r"C:\Users\bapt0\Desktop\Seq_identifier\ref_seq_cleaned.fasta"

files_directory = []
files_names = []
for file_name in glob.glob(general_dir+"/**/*.fas", recursive=True):
    files_directory.append(file_name)
    head, tail = os.path.split(file_name)
    files_names.append(tail.split('.', 1)[0])

# pick the first sequence of the reference file
with open(ref_file_dir, "r") as ref:
    ref_cont = ref.read().split('>')[1:]
    ref_name = ref_cont[0].split('\n', 1)[0]
    ref_seq = ref_cont[0].split('\n', 1)[1].replace('\n', '')

ref_names, ref_seqs = [], []
for x in ref_cont:
    ref_names.append(x.split('\n', 1)[0])
    ref_seqs.append(x.split('\n', 1)[1])

for file in files_directory:
    # file=files_directory[1]
    # the directory to the aligned file
    aligned_dir = file.replace(".fas", "_aligned.fasta")
    
    # retrieve that sequence in the aligned file
    with open(aligned_dir, "r") as align:
        aligned_ref = align.read().split('>'+ref_name)[1].split('>')[0].replace('\n', '')
    
    # copy the content of the aligned file
    with open(aligned_dir, "r") as align:
        cont = align.read().split('>')[1:]
    
    name, seq, ref_seq_aligned = [], [], []
    for x in cont:
        name.append(x.split('\n', 1)[0])
        if name[-1] in ref_names:
            ref_seq_aligned.append(x.split('\n', 1)[1].replace('\n', ''))
            name.pop(-1)
        else:
            seq.append(x.split('\n', 1)[1].replace('\n', ''))
    
    i=0
    while i<len(seq)-1:
        score_max = 0
        for sequence in ref_seq_aligned:
            score = 0
            for j in range(len(seq[i])):
                if seq[i][j] == sequence[j]:
                    score += 1
            if score > score_max:
                score_max = score
        if score_max < 0.5*len(seq[i]):
            print(name[i], score_max, len(seq[i]))
            seq.pop(i)
            name.pop(i)
        else:
            i += 1
    
    loci = []
    for i in range(len(seq[0])):
        locus = ''
        for x in seq:
            locus += x[i]
        loci.append(locus)
    
    aligned_ref = list(aligned_ref)
    i = 0
    while i < len(ref_seq)-1:
        if ref_seq[i] != aligned_ref[i]:
            aligned_ref.pop(i)
            loci.pop(i)
        else:
            i += 1
        # print(f"{i}: {ref_seq[i]=} != {aligned_ref[i]=}")
    loci = loci[:len(ref_seq)]
    
    with open(aligned_dir.replace(".fasta", "_cleaned.fasta"), "w") as cln:
        for i in range(len(loci[0])):
            tag = ''
            for x in loci:
                tag = tag + x[i]
            if name[i] not in ref_names:
                cln.write('>' + name[i] + '\n' + tag + '\n')
