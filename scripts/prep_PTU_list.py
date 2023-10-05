import pandas as pd  
import re  

plasmids = pd.read_csv('../data/table_50PTUs_RS200.csv', sep='\t')

top_50_ptus = []
with open('../data/PTUs.txt', 'r') as f:
	for line in f.readlines():
		top_50_ptus.append(line.strip('\n'))

print(top_50_ptus)

for ptu in top_50_ptus:
	ids = plasmids[plasmids['PTU']==ptu]['ID']
	ptu_renamed = re.sub("\\/", "", ptu)
	with open('../data/'+ptu_renamed+'_accs.txt', 'w') as f:
		for id in ids:
			f.write('%s\n' % id)