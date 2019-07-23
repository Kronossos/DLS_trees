#!/usr/bin/env bash

echo "Creating directories..."
rm -rf test_sets
rm -rf test_res
rm -rf test_scenario

mkdir test_sets
mkdir test_res
mkdir test_scenario

touch test_res/rme_PG.txt
touch test_res/rme_FHS.txt
touch test_res/heu_PG.txt
touch test_res/heu_FHS.txt

for i in {1..1000}
do
   echo $i
   rm -rf test_scenario
   mkdir test_scenario
   mkdir test_scenario/FHS
   mkdir test_scenario/PG

   python3 Random_tree.py >test_sets/$i"_full.txt"
   python2 ../rme/rme.py -P test_sets/$i"_full.txt" | sed "s/MEscore/$i/g" >> test_res/rme_PG.txt
   python2 ../rme/rme.py -F test_sets/$i"_full.txt" |tail -n 1 | sed "s/MEscore/$i/g" >> test_res/rme_FHS.txt

   spec=`cat test_sets/$i"_full.txt"|tail -n 1`
   genes=`cat test_sets/$i"_full.txt"|sed '2,49!d'`
   name=0

   for g in $genes
   do
      ((name+=1))
      ../dlsgen/dlsgen2 -gs -da -pilaEv  $g $spec >test_scenario/PG/PG_$name.txt
      ../dlsgen/dlsgen2 -ga -da -pilaEv  $g $spec >test_scenario/FHS/FHS_$name.txt
   done

   python3 ../reader.py `ls -d /home/alek/PycharmProjects/DLS_trees/test_module/test_scenario/FHS/*`|tail -n 1 | sed "s/MEscore/$i/g" >> test_res/heu_FHS.txt
   python3 ../reader.py `ls -d /home/alek/PycharmProjects/DLS_trees/test_module/test_scenario/PG/*` |tail -n 1 | sed "s/MEscore/$i/g" >> test_res/heu_PG.txt
done

python3 draw_test.py
