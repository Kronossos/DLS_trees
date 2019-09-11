from Tree import *
import random
import copy


def generate_tree(leafs=[1,2,3,4],multi=True,delete_chances=0.8):
    names=[Leaf(i) for i in leafs]

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


def tree_to_text(node):
    text="({},{})"
    left=""
    right=""

    if node.left.is_leaf():
        left="{}".format(node.left.label())
    else:
        left=tree_to_text(node.left)

    if node.right.is_leaf():
        right = "{}".format(node.right.label())
    else:
        right = tree_to_text(node.right)

    return text.format(left,right)


def create_random_set(species_tree_size=15,data_size={3:25,4:8,5:6,6:4,8:3,9:2}):
    codes = []
    for x in range(97, 123):
        codes.append(chr(x))
    spec_labels=codes[:species_tree_size]
    data_file=str(sum(data_size.values()))
    # spec_labels = [x for x in range(species_tree_size)]
    spec = generate_tree(spec_labels, False)
    for tree_size, subset_size in data_size.items():
        for x in range(subset_size):
            subset_labels = random.sample(spec_labels,tree_size)
            subset_tree = generate_tree(subset_labels, False)
            data_file += "\n{}".format(tree_to_text(subset_tree))
    data_file += "\n#S.T.\n{}".format(tree_to_text(spec))
    return data_file


print(create_random_set())


