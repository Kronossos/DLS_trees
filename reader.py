#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import collections


class TreeNode:
    def __init__(self, id: int, level: int, change_id: list, change_type: list, duplication_prefix: list):
        self.id = id
        self.duplication_prefix = duplication_prefix
        self.level = level

        connections = {}

        # Tmove=1, Clost=0
        # {1: [list of trees], 0: [list of trees]}
        # or {tree: connection_type,  tree_2: connection_type, ...} <---- current
        # Store ID or object?
        for i in range(len(change_id)):
            connections[change_id[i]] = change_type[i]
        self.connections = connections

    def __str__(self):
        return str(self.id) + ":" + str(self.duplication_prefix)


class Scenario:
    def __init__(self, filename):

        # Should scenario store dict with levels (easier print then)?
        nodes = collections.OrderedDict()

        with open(filename, "r") as scenario_file:
            for line in scenario_file:
                id, level, duplication_prefix, change, change_type = line.strip().split(" ")  # tab,space or comma?
                new_tree = TreeNode(int(id), int(level), eval(change), eval(change_type), eval(duplication_prefix))
                nodes[id] = new_tree
        self.nodes = nodes

    def __iter__(self):
        for tree_id in reversed(self.nodes):
            yield self.nodes[tree_id]

    def __str__(self):
        pass


class AllScenarios:
    def __init__(self, files_list):
        scenarios = []
        for file in files_list:
            scenarios.append(Scenario(file))
        self.scenarios = scenarios

    def __iter__(self):
        for scenario in self.scenarios:
            yield scenario

    def __str__(self):
        pass


def test():
    # files = sys.argv[1:]
    # files = ["scenario_1", "scenario_2"]
    files = ["pg_sc_1", "pg_sc_2"]
    a = AllScenarios(files)
    for sc in a:
        print()
        for tree in sc:
            print(tree)


test()
