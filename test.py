x=[["a","b"],["c","d"],["1","2"],["X","Y"]]
import itertools

for i in itertools.product(*x):
    print(i)