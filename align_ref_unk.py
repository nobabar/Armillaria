from Bio.Align.Applications import MuscleCommandline
import glob
import os


general_dir = r"C:\Users\bapt0\Desktop\Stage_cyril\ITSarmillaire\test"

# the directory to the file containing the reference sequences, must be .fas
reference_file_dir = r"C:\Users\bapt0\Desktop\Seq_identifier\ref_seq_cleaned.fasta"


# the directory to the muscle script that will perform the alignement
# can be download from here https://www.drive5.com/muscle/downloads.htm
muscle_dir = r"C:\Users\bapt0\Desktop\Seq_identifier\muscle3.8.31_i86win32.exe"
# If you don't know where is your muscle.exe
# for root, dir, files in os.walk("C:\\Users"):    #you can partially complete the path
#     if "muscle3.8.31_i86win32.exe" in files:
#         muscle_dir=os.path.join(root, "muscle3.8.31_i86win32.exe")


files_directory = []
files_names = []
for file_name in glob.glob(general_dir+"/**/*.fas", recursive=True):
    files_directory.append(file_name)
    head, tail = os.path.split(file_name)
    files_names.append(tail.split('.', 1)[0])


for file in files_directory:
    # the directory to the alingned and unaligned files, making sure they are empty
    unaligned_dir = file.replace(".fas", "_aligned.fasta")
    open(unaligned_dir, "w+").close()
    aligned_dir = file.replace(".fas", "_aligned.fasta")
    # open(aligned_dir, "w+").close()
    
    with open(unaligned_dir, "a") as align:
        with open(reference_file_dir, "r") as ref:
            for line in ref:
                align.write(line)
        with open(file, "r") as unk:
            for line in unk:
                align.write(line)
    
    muscle_cline = MuscleCommandline(muscle_dir, input=unaligned_dir, out=aligned_dir)
    stdout, stderr = muscle_cline()
