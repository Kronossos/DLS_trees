#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys

files = sys.argv[1:]


class TreeNode:
    def __init__(self, id: int,change_id: tuple,change_type: tuple, duplication_prefix: tuple):
        self.id = id
        self.duplication_prefix = duplication_prefix

        connections = {}

        # Tmove=1, Clost=0
        # {1: [list of trees], 0: [list of trees]}
        # or {tree: connection_type,  tree_2: connection_type, ...} <---- current
        for i in range(len(change_id)):
            connections[change_id[i]]=change_type[i]






class Scenario:
    def __init__(self, filename):
        nodes = {}

        with open(filename, "r") as scenario_file:
            for line in scenario_file:
                id, change, change_type, duplication_prefix = line.split("\t")
                new_tree = TreeNode(int(id), tuple(change), tuple(change_type), tuple(duplication_prefix))
                nodes[id]=new_tree
