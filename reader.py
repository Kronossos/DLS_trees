#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import collections
import random


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
        self.name=filename

    def __iter__(self):
        for tree_id in reversed(self.nodes):
            yield self.nodes[tree_id]

    def __str__(self):
        return str(self.name)


class AllScenarios:
    def __init__(self, files_list):
        scenarios = {}
        for file in files_list:
            scenarios[file] = Scenario(file)
        self.scenarios = scenarios

    def __iter__(self):
        for scenario in self.scenarios:
            yield self.scenarios[scenario]

    def random_scenario(self):
        random_trees={}
        for scenario in self:
            chosen_tree=random.choice(list(scenario.nodes.values()))
            random_trees[(scenario.name,chosen_tree.id)] = chosen_tree.duplication_prefix
        return random_trees


    def rate_scenario(self,chosen_trees):
        all_trees = list(chosen_trees.values())
        max_of_dup = [max(x) for x in zip(*all_trees)]
        return max_of_dup

    def select_scenarios(self,chose_fun = random_scenario,iter=100000):

        min_cost=float("inf")
        min_scenatio=[]

        for i in range(iter):
            chosen_scenario = chose_fun(self)
            rated_scenatio = self.rate_scenario(chosen_scenario)

            current_cost=sum(rated_scenatio)
            if current_cost<min_cost:
                min_cost = current_cost
                min_scenatio=rated_scenatio
        return min_scenatio,min_cost

    def __str__(self):
        pass


def test():
    files = sys.argv[1:]
    a = AllScenarios(files)
    # for sc in a:
    #     print(sc)
    #     for tree in sc:
    #         print(tree)

    print(a.select_scenarios())

test()
