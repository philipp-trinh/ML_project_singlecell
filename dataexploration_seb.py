import os

from dataexploration import make_adata
print("cwd: ", os.getcwd())
data = make_adata()
print(data)