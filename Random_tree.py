from Tree import *
import random


def generate_tree(leafs,multi=True,delete_chances=0.8):
    names=[Leaf(i) for i in range(leafs)]

    while len(names)>1:
        left = random.choice(names)
        right = random.choice(names)

        if not multi and left == right:
            continue
        if left == right and random.random() < delete_chances:
            continue

        names.append(BinNode(left, right))

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

print(generate_tree(3,True,0.7))

print(generate_tree(3,False))

