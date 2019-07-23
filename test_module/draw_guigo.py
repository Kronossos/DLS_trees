# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np




guigo_1_PG={"Heurystyka": 6,"RME": 6}
guigo_1_PG = pd.DataFrame(guigo_1_PG.items(), columns=['Program', 'MEscore'])

guigo_2_PG={"Heurystyka": 5,"RME": 5}
guigo_2_PG = pd.DataFrame(guigo_2_PG.items(), columns=['Program', 'MEscore'])

guigo_1_FHS={"Heurystyka": 4,"RME": 4}
guigo_1_FHS = pd.DataFrame(guigo_1_FHS.items(), columns=['Program', 'MEscore'])

guigo_2_FHS={"Heurystyka": 4,"RME": 4}
guigo_2_FHS = pd.DataFrame(guigo_2_FHS.items(), columns=['Program', 'MEscore'])


plt.close('all')
guigo_1_PG.plot.bar(x='Program',y='MEscore',legend=False,title="Guigo 1: Model PG",ylim=(0,9))
plt.xlabel('')
plt.ylabel('MEscore')
plt.savefig('G1_PG.png', dpi=300, format='png', bbox_inches='tight')

plt.close('all')
guigo_2_PG.plot.bar(x='Program',y='MEscore',legend=False,title="Guigo 2: Model PG",ylim=(0,9))
plt.xlabel('')
plt.ylabel('MEscore')
plt.savefig('G2_PG.png', dpi=300, format='png', bbox_inches='tight')

plt.close('all')
guigo_1_FHS.plot.bar(x='Program',y='MEscore',legend=False,title="Guigo 1: Model FHS",ylim=(0,9))
plt.xlabel('')
plt.ylabel('MEscore')
plt.savefig('G1_FHS.png', dpi=300, format='png', bbox_inches='tight')

plt.close('all')
guigo_2_FHS.plot.bar(x='Program',y='MEscore',legend=False,title="Guigo 2: Model FHS",ylim=(0,9))
plt.xlabel('')
plt.ylabel('MEscore')
plt.savefig('G2_FHS.png', dpi=300, format='png', bbox_inches='tight')
