import numpy as np
import math
import matplotlib.pyplot as plt
fr=open('feature_importance.txt','r')
line = fr.readlines()

imp=np.zeros((509))
for eachline in line:
    line=eachline.strip().split(',')
    site=int(line[0])
    value=int(line[1])
    imp[site]=value
fr.close()

importance=imp
plt.bar(range(0,204),importance[0:204],1,color='#FFCE00',label='4-bits Binary')
plt.bar(range(204,408),importance[204:408],1,color='#0375B4',label='OPF')
plt.bar(range(408,488),importance[408:488],1,color='#007849',label='k-mers')
plt.bar(range(488,490),importance[488:490],1,color='#FF0033',label='Site Location')
plt.bar(range(490,497),importance[490:497],1,color='#CCCC00',label='Information Theory')
plt.bar(range(497,509),importance[497:509],1,color='#262228',label='SNP')

plt.xlabel('feature')
plt.ylabel('importance-score')
plt.xlim(0,510)
plt.ylim(0,1200)
plt.legend(loc='upper left')
plt.show()

for i in importance:
    print(i)
