# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

heu_FHS={}
with open("test_res/heu_FHS.txt") as file:
    for line in file:
        id,value=line.split()
        heu_FHS[id]=int(value)
heu_FHS = pd.DataFrame(heu_FHS.items(), columns=['ID', 'MEscore'])

heu_PG={}
with open("test_res/heu_PG.txt") as file:
    for line in file:
        id,value=line.split()
        heu_PG[id]=int(value)
heu_PG = pd.DataFrame(heu_PG.items(), columns=['ID', 'MEscore'])

rme_FHS={}
with open("test_res/rme_FHS.txt") as file:
    for line in file:
        id,value=line.split()
        rme_FHS[id]=int(value)
rme_FHS = pd.DataFrame(rme_FHS.items(), columns=['ID', 'MEscore'])


rme_PG={}
with open("test_res/rme_PG.txt") as file:
    for line in file:
        id,value=line.split()
        rme_PG[id]=int(value)
rme_PG = pd.DataFrame(rme_PG.items(), columns=['ID', 'MEscore'])


plt.close('all')
PG_score = heu_PG['MEscore'] - rme_PG['MEscore']
PG_score=PG_score.value_counts()
print(PG_score)
PG_score.plot(kind='bar',title="Model PG")
plt.xlabel('Różnica')
plt.ylabel('Ilość prób')
plt.savefig('PG.png', dpi=300, format='png', bbox_inches='tight')


plt.close('all')
FHS_score = heu_FHS['MEscore'] - rme_FHS['MEscore']
FHS_score=FHS_score.value_counts()
print(FHS_score)
FHS_score.plot(kind='bar',title="Model PG")
plt.xlabel('Różnica')
plt.ylabel('Ilość prób')
plt.savefig('FHS.png', dpi=300, format='png', bbox_inches='tight')