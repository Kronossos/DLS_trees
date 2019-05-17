#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#amor@mimuw Kronossos alek1
import time
import sys
import collections
import random
import itertools

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
                new_tree.owner = filename
                nodes[id] = new_tree
        if not nodes:
            raise IndexError('Cannot create scenario from an empty file!')
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
        all_files=len(files_list)
        counter=0
        for file in files_list:
            try:
                scenarios[file] = Scenario(file)
            except:
                pass
                counter += 1
        self.scenarios = scenarios

        corrupted_percent=counter*100/all_files
        print("Data loaded. {}% of the data was corrupted: {} of {}.".format(corrupted_percent,counter,all_files))

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
        #Może zjeść pamięć dla dużych zbiorów. DOTESTOWAĆ!
        for scenario in self:
            all_dup_pref = [tree.duplication_prefix for tree in scenario]
            max_trees.append(self.rate_scenario(all_dup_pref))
        max_tree = self.rate_scenario(max_trees)

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
                        for i in range(len(tree.duplication_prefix)):
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
                        for i in range(len(tree.duplication_prefix)):
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
                        for i in range(len(tree.duplication_prefix)):
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

    def subset_calc(self,buffer_size=10):
        # TO CHANGE! NAIVE WAY OF DEALING WITH GETTING TREES!

        current_set = [[tree.duplication_prefix for tree in scenario ] for scenario in self]

        while current_set == 1:
            current_set = 1 #split by buffer size



    def check_all(self,trees=False):
        if not trees:
            # TO CHANGE! NAIVE WAY OF DEALING WITH GETTING TREES!
            trees = [[tree.duplication_prefix for tree in scenario] for scenario in self]

        min_cost = float("inf")
        min_scenatio = []

        size=1
        for scenario in trees:
            size*=len(scenario)
        print(size)

        counter=0
        for chosen_scenario in itertools.product(*trees):
            rated_scenatio = self.rate_scenario(chosen_scenario)
            current_cost = sum(rated_scenatio)
            if current_cost < min_cost:
                min_cost = current_cost
                min_scenatio = rated_scenatio
            counter+=1
            i = counter * 100 / size
            sys.stdout.write("\rCompleted {} percent.".format(i))
            sys.stdout.flush()

        return min_scenatio, min_cost

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


    # a.subset_calc()



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

    # print("FULL")
    # start_time = time.time()
    # print(a.check_all())
    # print("Done in {} .".format(time.time() - start_time))




test()
