# http://telliott99.blogspot.com/2010/03/clustering-with-upgma.html

from Wunsch import Alignment
from Tree import *
import random

def reconstruct_ancestors(self):
    if self.left.label():
        if self.right.label():

            container = Alignment(self.left.label(), self.right.label())[1]

            new_label = ""
            for i in range(len(container[0])):
                new_label += random.choice(container)[i]

            self.set_label(new_label.replace("_", ""))
        else:
            self.right.reconstruct_ancestors()
            self.reconstruct_ancestors()
    else:
        self.left.reconstruct_ancestors()
        self.reconstruct_ancestors()
def history_cost(self):
    tree_cost = 0
    for child in iter(self):
        try:
            if not child.is_leaf():
                tree_cost += Alignment(child.left.label(), child.right.label())[0]
        except TypeError:
            return None
    return tree_cost


BinNode.reconstruct_ancestors = reconstruct_ancestors
BinNode.history_cost = history_cost


def Reconstruct_History(sequences, cost=Alignment):
    F = [[float('Inf') for col in range(len(sequences))] for row in range(len(sequences))]

    for row in range(len(sequences)):
        col = 0
        while col < row:
            F[row][col] = cost(sequences[row], sequences[col])[0]
            col += 1

    sequences_list = []
    for seq in sequences:
        sequences_list.append(Leaf(seq))

    while len(F) > 1:

        col_min_val = float("inf")
        for row in range(len(F)):
            col = min(F[row])
            if col < col_min_val:
                col_min_val = col
                col_min = F[row].index(min(F[row]))
                row_min = row

        container = cost(sequences_list[row_min].label(), sequences_list[col_min].label())[1]
        new_label = ""
        for i in range(len(container[0])):
            new_label += random.choice(container)[i]
        new_label = new_label.replace("_", "")

        new_node = BinNode(sequences_list[row_min], sequences_list[col_min])
        new_node.set_label(new_label)
        sequences_list.append(new_node)

        for index in sorted([col_min, row_min], reverse=True):
            del sequences_list[index]
            del F[index]
            for row in range(len(F)):
                del F[row][index]

        for row_number in range(len(F)):
            F[row_number].append(float("Inf"))

        new_row = []
        for seq_number in range(len(sequences_list) - 1):
            new_row.append(cost(new_label, sequences_list[seq_number].label())[0])
        new_row.append(float("Inf"))
        F.append(new_row)

    return BinTree(sequences_list[0])


sequences = []
sequences.append('ATAAGGCCAT')
sequences.append('ACGTTAAGCGT')
sequences.append('ATAAAGCCAT')
sequences.append('ACGTAAGCTT')
sequences.append('ACGTAAAGCGT')

t = BinTree(BinNode(BinNode(Leaf(sequences[0]), Leaf(sequences[1])),
                    BinNode(Leaf(sequences[2]), BinNode(Leaf(sequences[3]), Leaf(sequences[4])))))

print(t)
print("\n------------------------\n")
t.reconstruct_ancestors()
print(t)
print("\n------------------------\n")
t = Reconstruct_History(sequences)
print(t)
print(t.history_cost())
