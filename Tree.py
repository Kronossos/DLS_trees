# -*- coding: utf-8 -*-
import random


class Leaf:
    def __init__(self, value):
        self.value = value
        self.color = "white"

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
        return 1

    def set_label(self, value):
        self.value = value


class BinNode:
    def __init__(self, left, right, value=None):
        self.father = None

        self.left = left
        left.father = self

        self.right = right
        right.father = self

        self.value = value
        self.color = "white"

    def __iter__(self):
        for node in self.left:
            yield node
        for node in self.right:
            yield node
        yield self

    def __str__(self, level=0, blank=False, sign="  "):
        if not blank: blank = []
        switch = False
        line = list(" " * level + sign + "█>")

        # if self.label():
        #     line += self.label()

        for fill in blank[:]:
            if switch or "".join(line[fill]) == "└":
                switch = True
                blank.pop(-1)
            elif line[fill] == " ":
                line[fill] = "│"

        tree = "".join(line)

        if self.label():
            tree += str(self.label())

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

    def __init__(self, node):
        if node.__class__.__name__ in ["BinNode", "Leaf"]:
            self.node = node
            self.label_dict = {}
        else:
            exit("Given root is not possible part of BinTree.")

    def __str__(self):
        if self.root: return self.node.__str__()

    def __iter__(self):
        return iter(self.node)

    def root(self):
        return self.node

    def create_set_label(self):
        for x in self:
            if not x.is_leaf():
                left_labels = set(x.son("L").label())
                right_labels = set(x.son("R").label())

                x.set_label(tuple(right_labels.union(left_labels)))

            else:
                new_label = (x.label(),)
                x.value = new_label

    def create_string_label(self, sign="N_"):
        count = 0
        for x in self:
            if not x.is_leaf():
                x.set_label(sign + str(count))
                count += 1

    def create_list_label(self):
        for x in self:
            if not x.is_leaf():
                left_labels = x.son("L").label()
                right_labels = x.son("R").label()

                x.set_label(right_labels + left_labels)

            else:
                new_label = [x.label()]
                x.value = new_label

    def label_to_node(self, label_function=create_set_label):

        label_dict = {}

        for node in self:
            label_dict[node.label()] = node

        self.label_dict = label_dict

    def label_dict(self):
        return self.label_dict

    def son_set(self):
        son_set = []
        for v in self:
            if not v.is_leaf():
                son_set.append((v.son("L"), v.son("R")))
        return son_set