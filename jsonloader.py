import json


def load_to_dict(file):
    with open(file, "r") as j:
        x = json.loads(j.read())
        x = key_int(x)
    return x

def key_int(dic):
    d = {int(k):str(l) for k,l in dic.items()}
    return d