import sys
import re
file_name = sys.argv[1]


word_set=set()
with open(file_name) as file:
    for line in file:
        line = re.findall("\w+",line)
        word_set.update(line)

# print(word_set)
# print(len(word_set))

codes=[]
for x in range(97,123):
    codes.append(chr(x))
for x in range(48,58):
    codes.append(chr(x))

codes=codes[:len(word_set)]

translate=zip(word_set,codes)

#sed 's/acoe/a/g;s/agna/b/g;s/amph/c/g;s/anne/d/g;s/arth/e/g;s/aves/f/g;s/chlo/g/g;s/chon/h/g;s/echi/i/g;s/embr/j/g;s/fung/k/g;s/mamm/l/g;s/moll/m/g;s/oste/n/g;s/prot/o/g;s/rept/p/g'

command="sed \'"

for key,val in translate:
    # print(key,val)
    command += "s/{}/{}/g;".format(key,val)
command += "\'"
# command += "\' {}".format(file_name)

print(command)

#sed 's/CAEEL/a/g;s/ARATH/b/g;s/BRARE/c/g;s/ORYSA/d/g;s/SCHMA/e/g;s/MOUSE/f/g;s/TETNG/g/g;s/SCHPO/h/g;s/APIME/i/g;s/MONDO/j/g;s/FUGRU/k/g;s/CAEBR/l/g;s/HUMAN/m/g;s/MACMU/n/g;s/CIOIN/o/g;s/DROME/p/g;s/ANOGA/q/g;s/DROPS/r/g;s/CANFA/s/g;s/GASAC/t/g;s/YEAST/u/g;s/AEDAE/v/g;s/ORYLA/w/g;s/PANTR/x/g;s/CHICK/y/g;s/XENTR/z/g;s/BOVIN/0/g;s/RAT/1/g;' sets/treefam/s.txt

