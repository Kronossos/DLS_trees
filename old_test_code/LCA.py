from Tree import *
import random
from Tarjan import *
from Random_tree import *



labels=[]
labels.append('KING')
labels.append('SHEEP')
labels.append('DRAGON')

spec = BinTree(BinNode(Leaf(labels[2]),BinNode(Leaf(labels[0]),Leaf(labels[1]))))
gen = BinTree(BinNode(Leaf(labels[2]),BinNode( BinNode(Leaf(labels[2]),Leaf(labels[1])), BinNode(Leaf(labels[0]),Leaf(labels[1])))))


gen.create_set_label()
son_set=gen.son_set()




print(gen)
# TarjanOLCA(gen,xd)
print()


spec.create_set_label()
spec.label_to_node()
print(spec)
print()
TarjanOLCA(spec.root(),son_set,label_dict=spec.label_dict)







# ala=gen.son("R").son("R")
# for key in ala.__dict__:
#     try: print(key,": ",ala.__dict__[key].label())
#     except: print(key,": ",ala.__dict__[key])