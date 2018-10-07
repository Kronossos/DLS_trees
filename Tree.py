# -*- coding: utf-8 -*-
import random

class Leaf:
    def __init__(self, value):
        self.value = value

    def __iter__(self):
        yield self

    def __str__(self, level=0, blank=False, sign="  "):
        if not blank:
            blank = []
        line = list(" " * level + sign + str(self.value))
        for fill in blank:
            try:
                if line[fill] == " ":
                    line[fill] = "│"
            except:
                pass
        return "".join(line)

    def is_leaf(self):
        return True

    def label(self):
        return self.value

    def height(self):
        return 0

class BinNode:
    def __init__(self, left=None, right=None,value=None):
        # for node in [left, right]:
            # if node.__class__.__name__ not in ["BinNode", "Leaf"]:
            #     exit("Given nodes are not possible part of BinNode.")
        self.left = left
        self.right = right
        self.value = value

    def __iter__(self):
        for node in self.left:
            yield node
        yield self
        for node in self.right:
            yield node

    def __str__(self, level=0, blank=False, sign="  "):
        # Python jest o tyle upierdliwy, że pamięta odnośnik do blanka.
        # Jeśli się wywoła "czystą" wersję to lista zawierać będzie jakieś elementy.
        # Taki sposób to omija.
        if not blank: blank = []
        switch = False
        line = list(" " * level + sign + "█>")

        if self.label():
            line += self.value()

        # Troszkę sztuczka, ponieważ python nawet jeśli wewnąrze pętli zacznie kasować elementy z końca
        # to sam for i tak po nich przejdzie. Więc super przydatne w tym momencie.
        for fill in blank[:]:
            if switch or "".join(line[fill]) == "└":
                switch = True
                blank.pop(-1)
            elif line[fill] == " ":
                line[fill] = "│"

        tree = "".join(line)

        if not self.left.is_leaf(): blank.append(level + 2)

        tree += "\n" + self.left.__str__(level + 2, blank, "├-")
        tree += "\n" + self.right.__str__(level + 2, blank, "└-")
        return tree

    def is_leaf(self):
        return False

    def son(self, which):
        if which == "L":
            return self.left
        elif which == "R":
            return self.right
        else:
            print("Unknow son!")

    def set_label(self, value):
        self.value = value

    def label(self):
        return self.value

class BinTree:
    def __getattr__(self, name):
        print(self.node.locals())
        return self.node.locals()[name]()

    def __init__(self, node):
        if node.__class__.__name__ in ["BinNode", "Leaf"]:
            self.node = node
        else:
            exit("Given root is not possible part of BinTree.")

    def __str__(self):
        if self.root: return self.node.__str__()

    def __iter__(self):
        return iter(self.node)

    def root(self):
        return self.node

