#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import time
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

    def __iter__(self):
        yield self.duplication_prefix


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
        self.name = filename

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
        random_trees = {}
        for scenario in self:
            chosen_tree = random.choice(list(scenario.nodes.values()))
            random_trees[(scenario.name, chosen_tree.id)] = chosen_tree.duplication_prefix
        return random_trees

    def decrease_vector(self, select_type="start"):
        max_trees = []
        for scenario in self:
            all_dup_pref = [tree.duplication_prefix for tree in scenario]
            max_trees.append(self.rate_scenario(all_dup_pref))
        max_tree = self.rate_scenario(max_trees)
        print(max_tree)

        if select_type == "start":
            index = 0
            while index < len(max_tree):
                if max_tree[index] == 0:
                    index += 1
                    continue

                max_tree_temp = max_tree[:]
                max_tree_temp[index] -= 1

                for scenario in self:
                    for tree in scenario:
                        for i in range(len(tree.duplication_prefix) - 1, index - 1, -1):
                            if max_tree_temp[i] - tree.duplication_prefix[i] < 0:
                                break
                        else:
                            break
                    else:
                        index += 1
                        break
                else:
                    max_tree = max_tree_temp
            return max_tree, sum(max_tree)

        if select_type == "end":

            index = len(max_tree) - 1
            while index > 0:
                if max_tree[index] == 0:
                    index -= 1
                    continue

                max_tree_temp = max_tree[:]
                max_tree_temp[index] -= 1

                for scenario in self:
                    for tree in scenario:
                        for i in range(index, -1, -1):
                            if max_tree_temp[i] - tree.duplication_prefix[i] < 0:
                                break
                        else:
                            break
                    else:
                        index -= 1
                        break
                else:
                    max_tree = max_tree_temp
            return max_tree, sum(max_tree)

        if select_type == "random":

            index_list = [x for x in range(len(max_tree)) if x != 0]

            while index_list:
                index_list_position = random.randint(0, len(index_list) - 1)
                index = index_list[index_list_position]

                max_tree_temp = max_tree[:]
                max_tree_temp[index] -= 1

                for scenario in self:
                    for tree in scenario:
                        for i in index_list:
                            if max_tree_temp[i] - tree.duplication_prefix[i] < 0:
                                break
                        else:
                            break
                    else:
                        index_list.pop(index_list_position)
                        break
                else:
                    max_tree = max_tree_temp
            return max_tree, sum(max_tree)

    def rate_scenario(self, chosen_trees):
        all_trees = list(chosen_trees)
        max_of_dup = [max(x) for x in zip(*all_trees)]
        return max_of_dup

    def select_scenarios(self, chose_fun=random_scenario, iter_num=1000):

        min_cost = float("inf")
        min_scenatio = []

        for i in range(iter_num):
            chosen_scenario = chose_fun(self)
            rated_scenatio = self.rate_scenario(chosen_scenario.values())

            current_cost = sum(rated_scenatio)
            if current_cost < min_cost:
                min_cost = current_cost
                min_scenatio = rated_scenatio
        return min_scenatio, min_cost

    def __str__(self):
        pass


def test():
    files = sys.argv[1:]
    a = AllScenarios(files)
    # for sc in a:
    #     print(sc)
    #     for tree in sc:
    #         print(tree)

    print("total random")
    start_time = time.time()
    print(a.select_scenarios())
    print("Done in {} .".format(time.time() - start_time))

    print("start")
    start_time = time.time()
    print(a.decrease_vector())
    print("Done in {} .".format(time.time() - start_time))

    print("end")
    start_time = time.time()
    print(a.decrease_vector("end"))
    print("Done in {} .".format(time.time() - start_time))

    print("index random")
    start_time = time.time()
    print(a.decrease_vector("random"))
    print("Done in {} .".format(time.time() - start_time))


test()
