from LCA import BinNode, BinTree


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
