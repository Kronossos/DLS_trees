from Tree import *
import random
import copy
from Tarjan import *

random.seed(1)

def generate_tree(leafs,multi=True,delete_chances=0.8):
    names=[Leaf(i) for i in range(leafs)]

    while len(names)>1:
        left = random.choice(names)
        right = random.choice(names)

        if not multi and left == right:
            continue
        if left == right and random.random() < delete_chances:
            continue

        names.append(BinNode(left, copy.deepcopy(right)))

        if left != right:
            if right.is_leaf() and random.random() > delete_chances and multi:
                pass
            else:
                names.remove(right)

        if left.is_leaf() and random.random() > delete_chances and multi:
            pass
        else:
            names.remove(left)



    return names[0]

gen=generate_tree(6,True,0.75)

spec=generate_tree(3,False)


print(gen)
print(spec)

count=0
for x in gen:
    if not x.is_leaf():
        x.set_label("G_"+str(count))
        count += 1

count=0
for x in spec:
    if not x.is_leaf():
        x.set_label("S_"+str(count))
        count += 1
print("\n")
# print(gen)
# print(spec)


