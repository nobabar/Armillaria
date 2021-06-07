from pandas_ods_reader import read_ods
import glob
import os


path=r'C:\Users\bapt0\Desktop\Stage_cyril\ITSarmillaire'


n=0
files_names = []
for n, file_name in enumerate(glob.iglob(path+"/**/*.seq", recursive=True)):
    head, tail = os.path.split(file_name)
    files_names.append(tail.split('_', 1)[0].split('-', 1)[0])
# print(*files_names, sep='\n')


path_ods = r'C:\Users\bapt0\Desktop\Seq_identifier\ITsaPasser_Ete2019.ods'
df = read_ods(path_ods, 1)
names=[]
for i in range(len(df)):
    if df.iloc[i]['ITSaPasser'] == 1:
        names.append(df.iloc[i]['ind'])


for i in names:
    if i in files_names:
        print(f"{i} has already been sequenced")
