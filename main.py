import numpy as np
import pandas as pd


def read_text():
    with open("./in/GSE83523_CE_blastomeres_rawcounts.txt", "r") as f:
        text = f.read()
    return text

def make_df_from_text():
    text = read_text().split("\n")
    items = [line.split("\t") for line in text]
    columnnames = [line[0] for line in items]
    df = pd.DataFrame(columns=columnnames)
    for line in items:
        if len(line) > 2:
            df[line[0]] = np.array(line[1:])
    pass

if __name__ == '__main__':
    a = pd.read_csv("./in/GSE83523_CE_blastomeres_rawcounts.txt", sep="\t")
    print(a)
