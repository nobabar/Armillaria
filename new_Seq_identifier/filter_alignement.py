import glob
import os


general_dir=r"C:\Users\bapt0\Desktop\Stage_cyril\ITSarmillaire"

# the directory to the file containing the reference sequences, must be .fas
reference_file_dir = r"C:\Users\bapt0\Desktop\Seq_identifier\ref_seq_cleaned.fasta"

files_directory=[]
files_names=[]
for file_name in glob.glob(general_dir+"/**/*.fas", recursive=True):
    files_directory.append(file_name)
    head, tail = os.path.split(file_name)
    files_names.append(tail.split('.', 1)[0])


# the directory to the alingned and unaligned files, making sure they are empty
aligned_dir = files_directory[0].rsplit("\\", 1)[0]+r"\aligned.fasta"


with open(reference_file_dir, "r") as ref:
    ref_cont = ref.read().split('>')[1:2]
    ref_name = ref_cont[0].split('\n',1)[0]
    ref_seq = ref_cont[0].split('\n',1)[1].replace('\n','')
with open(aligned_dir, "r") as align:
    align_cont = align.read().split('>'+ref_name)[1].split('>')[0]
    align_seq = align_cont.split('\n',1)[1].replace('\n','')


with open(aligned_dir, "r") as align:
    cont = align.read().split('>')[1:]


name, seq = [], []
for x in cont:
    name.append(x.split('\n',1)[0])
    seq.append(x.split('\n',1)[1].replace('\n',''))

    
loci = []
for i in range(len(seq[0])):
    locus = ''
    for x in seq:
        locus += x[i]
    loci.append(locus)

i=0
while i<len(ref_seq):
   if list(ref_seq)[i]!=list(align_seq)[i]:
       align_seq=align_seq[:i-1]+align_seq[i+1:]
       loci.pop(i)
   else:
       i+=1

