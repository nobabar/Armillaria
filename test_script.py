from Bio.Align.Applications import MuscleCommandline
import numpy as np
import glob
import os
import re


import jsonloader
import identifier

# Your working environment - need modifications!!

# the directory to the file you want to identify, must be a folder
identify_general_dir = r"C:\Users\bapt0\Desktop\Seq_identifier\Seq"
# can be seq or fas files
files_type = "seq"
# the directory to the file containing the reference sequences, must be .fas
reference_file_dir = r"C:\Users\bapt0\Desktop\Seq_identifier\ref_seq.fas"
reference_json = r".\data\index_ref.JSON"

# the directory to the alingned and unaligned files, making sure they are empty
unaligned_dir = r"C:\Users\bapt0\Desktop\Seq_identifier\Seq\unaligned.fasta"
open(unaligned_dir, "w+").close()
aligned_dir = r"C:\Users\bapt0\Desktop\Seq_identifier\Seq\aligned.fasta"
open(aligned_dir, "w+").close()

# the directory to the muscle script that will perform the alignement
# can be download from here https://www.drive5.com/muscle/downloads.htm
muscle_dir = r"C:\Users\bapt0\Desktop\Seq_identifier\muscle3.8.31_i86win32.exe"
# If you don't know where is your muscle.exe
# for root, dir, files in os.walk("C:\\Users"):    #you can partially complete the path
#     if "muscle3.8.31_i86win32.exe" in files:
#         muscle_dir=os.path.join(root, "muscle3.8.31_i86win32.exe")


# Checking the input settings
assert os.path.exists(
    identify_general_dir), f"Le chemin d'accès files_dir \'{identify_general_dir}\' n'est pas valide"
assert os.path.exists(
    reference_file_dir), f"Le chemin d'accès ref_dir \'{reference_file_dir}\' n'est pas valide"
assert os.path.exists(
    muscle_dir), f"Le chemin d'accès msucle_dir \'{muscle_dir}\' n'est pas valide"
identify_files = []
identify_files_directory = []
files_names = []
n = 0
for n, file_name in enumerate(glob.iglob(f"{identify_general_dir}\*.{files_type}", recursive=True)):
    identify_files_directory.append(file_name)
    head, tail = os.path.split(file_name)
    files_names.append(">"+tail.split('.', 1)[0])
print(f"We found {n+1} files(s) in .{files_type}", *files_names, sep='\n')


dict_ref=jsonloader.load_to_dict(reference_json)


m = 0
ostoyae = 0
cepistipes = 0
borealis = 0
mellea = 0
gallica = 0
tabescens = 0
with open(reference_file_dir, "r") as ref:
    for line in ref:
        if re.match("^>\w+", line):
            m += 1  # there are m reference sequences and n files to identify
            if "ostoyae" in line:
                ostoyae += 1
            if "cepistipes" in line:
                cepistipes += 1
            if "borealis" in line:
                borealis += 1
            if "mellea" in line:
                mellea += 1
            if "gallica" in line:
                gallica += 1
            if "tabescens" in line:
                tabescens += 1


def concat(concat_dir, id_files, ref_dir, f_type):
    '''
    Parameters
    ----------
    concat_dir : string (file directory)
        adir is the directory to the file which will contain both the reference
        files and the file to identify.
    id_dir : sring (folder directory)
        fdir is the directory to the file we need to identify.
    rdir : string (file directory)
        rdir is the directory to the reference file which contain all the
        reference sequences needed for the identification.
    ftype : string
        ftype is the type of file we will need to handle, it can either be
        .seq or .fas files.

    Returns
    -------
    The aligned.fas file is created and contain both the reference sequences
    and the sequences to analyse

    '''
    with open(concat_dir, "a") as align:
        with open(ref_dir, "r") as ref:
            for line in ref:
                align.write(line)
        for files_name in id_files:
            with open(files_name, "r") as seq:
                if f_type == "seq":
                    head, tail = os.path.split(files_name)
                    align.write(f"\n>{tail.split('.', 1)[0]}\n")
                for lines in seq:
                    align.write(lines.strip())


def align(muscle_exe, infile, outfile):
    '''
    Parameters
    ----------
    file : string (file directory)
        adir is the directory to the file which contain both the reference
        files and the file to identify, it will be used as both in and out file
        for the muscle alignment.
    muscle_exe : string (file directory)
        muscle_exe is the directory to the muscle3.8.31_i86win32.exe file that
        will perform the alignement of the sequences.
    Returns
    -------
    The sequences contained in the aligned.fas are aligned by muscle multi
    sequences alignement.
    '''
    muscle_cline = MuscleCommandline(muscle_exe, input=infile, out=outfile)
    stdout, stderr = muscle_cline()


def divide_chunks(l, n):
    '''
    Parameters
    ----------
    l : list
        l is the list you want to divide into chunks.
    n : integer
        n is the size of the chunks the list will be divideed into.
    Yields
    ------
    list of lists
        the original list is now divided into chunks (other lists) of size n.
    '''
    for i in range(0, len(l), n):
        yield l[i:i + n]


def hamming_distance(file, names):
    '''
    Parameters
    ----------
    file : string (file directory)
        file is the directory to the aligned file contaning all sequences.
    names : list
        names is the list of all sequences to identify.
    Returns
    -------
    distance : list
        distance is a list containing the hamming distances between each
        sequences to identify and each reference sequences.
    name : list
        name is the list of the references sequences in the correct order.
    '''
    name = []
    dna_ref = []
    dna_id = []
    with open(file, "r") as seq:
        line = seq.readline().strip()
        while re.match("^>\w+", line):
            name.append(line)
            line = seq.readline().strip()
            lines = []
            while re.match("[ATCG-]+", line):
                lines.append(line)
                line = seq.readline().strip()
            if name[-1] in names:
                dna_id.append(lines)
            else:
                dna_ref.append(lines)
                name.pop()
    distance = []
    for i in dna_id:
        end = re.compile("-{1,}$")
        i="".join(i)
        end.sub("", i)
        for j in dna_ref:
            j="".join(j)
            x = 0
            switch = True
            for (k, l) in zip(i, j):
                if switch == True:
                    if k != "-" and l != "-":
                        switch = False
                if switch == False:
                    if k=="-" or l=="-":
                        continue
                    elif k != l:
                        x += 1
            distance.append(x)
    return distance, name


def identify(list_distance, n_ref, no, nc, nb, nm, ng, nt):
    '''
    Parameters
    ----------
    list_distance : list
        this list is the one containg the hamming distances of all sequences.
    n_ref : integer
        n_ref is the number of reference sequences.
    no, nc, nb, nm, ng, nt : integer
        these integers are the number of each species contained among the
        references sequences.
    Returns
    -------
    specie : list
        this list contain the most probable specie identified for each file.
    '''
    list_distance = list(divide_chunks(list_distance, n_ref))
    specie_1 = []
    specie_2 = []
    for i in list_distance:
        index_1 = i.index(sorted(i)[0])
        if index_1 < no:
            specie_1.append(f"ostoyae {index_1}")
        if no <= index_1 < no+nc:
            specie_1.append(f"cepistipes {index_1}")
        if no+nc <= index_1 < no+nc+nb:
            specie_1.append(f"borealis {index_1}")
        if no+nc+nb <= index_1 < no+nc+nb+nm:
            specie_1.append(f"mellea {index_1}")
        if no+nc+nb+nm <= index_1 < no+nc+nb+nm+ng:
            specie_1.append(f"gallica {index_1}")
        if no+nc+nb+nm+ng <= index_1:
            specie_1.append(f"tabescens {index_1}")
        index_2 = i.index(sorted(i)[1])
        if index_2 < no:
            specie_2.append(f"ostoyae {index_2}")
        if no <= index_2 < no+nc:
            specie_2.append(f"cepistipes {index_2}")
        if no+nc <= index_2 < no+nc+nb:
            specie_2.append(f"borealis {index_2}")
        if no+nc+nb <= index_2 < no+nc+nb+nm:
            specie_2.append(f"mellea {index_2}")
        if no+nc+nb+nm <= index_2 < no+nc+nb+nm+ng:
            specie_2.append(f"gallica {index_2}")
        if no+nc+nb+nm+ng <= index_2:
            specie_2.append(f"tabescens {index_2}")
    return specie_1, specie_2


concat(unaligned_dir, identify_files_directory, reference_file_dir, files_type)
align(muscle_dir, unaligned_dir, aligned_dir)
distances, names = hamming_distance(aligned_dir, files_names)
species_1, species_2 = identify(distances, m, ostoyae, cepistipes,
                                borealis, mellea, gallica, tabescens)

d = dict(zip(list(range(1, len(species_1)+1)), zip(species_1, species_2)))
print("""
 ------------------------------------------------------------------------
|     unknown sample     || 1st armillary specie || 2nd armillary specie |
 ------------------------------------------------------------------------""")
for i in range(len(d)):
    row1 = (f"{names[i]}").center(24)
    if (d[i+1][0]).split()[0] == (d[i+1][1]).split()[0]:        
        row2 = (f"{d[i+1][0].split()[0]} {min(d[i+1][0].split()[1], d[i+1][1].split()[1])}").center(22)
        row3 = ("/").center(22)
    else:
        row2 = (f"{d[i+1][0]}").center(22)
        row3 = (f"{d[i+1][1]}").center(22)  
    print(f"|{row1}||{row2}||{row3}|")
print(""" ------------------------------------------------------------------------ """)
