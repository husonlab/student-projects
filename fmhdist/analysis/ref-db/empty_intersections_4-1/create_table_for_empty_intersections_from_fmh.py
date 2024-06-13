# hacky script to create markdown tables for the empty intersection analysis
import pandas as pd

genomes = ['GCA_022750515.1', 'GCA_008974285.1', 'GCA_014858625.1', 'GCA_018394375.1', 'GCA_001314365.1', 'GCF_000142945.1']

seed10="""GCA_001314365.1 vs GCA_018394375.1 is empty
GCA_008974285.1 vs GCA_014858625.1 is empty
GCA_008974285.1 vs GCA_022750515.1 is empty
GCA_014858625.1 vs GCA_018394375.1 is empty
GCA_018394375.1 vs GCA_022750515.1 is empty
GCA_018394375.1 vs GCF_000142945.1 is empty"""

seed20="""GCA_001314365.1 vs GCA_008974285.1 is empty
GCA_008974285.1 vs GCA_014858625.1 is empty
GCA_008974285.1 vs GCA_018394375.1 is empty
GCA_008974285.1 vs GCA_022750515.1 is empty
GCA_008974285.1 vs GCF_000142945.1 is empty
GCA_018394375.1 vs GCA_022750515.1 is empty"""

seed30="""GCA_001314365.1 vs GCA_008974285.1 is empty
GCA_001314365.1 vs GCA_018394375.1 is empty
GCA_008974285.1 vs GCA_014858625.1 is empty
GCA_008974285.1 vs GCA_022750515.1 is empty
GCA_014858625.1 vs GCA_018394375.1 is empty
GCA_018394375.1 vs GCA_022750515.1 is empty
GCA_018394375.1 vs GCF_000142945.1 is empty"""

seed40="""GCA_001314365.1 vs GCA_008974285.1 is empty
GCA_001314365.1 vs GCA_018394375.1 is empty
GCA_008974285.1 vs GCA_014858625.1 is empty
GCA_008974285.1 vs GCA_018394375.1 is empty
GCA_008974285.1 vs GCF_000142945.1 is empty
GCA_014858625.1 vs GCA_018394375.1 is empty
GCA_018394375.1 vs GCF_000142945.1 is empty"""

seed50="""GCA_001314365.1 vs GCA_008974285.1 is empty
GCA_001314365.1 vs GCA_018394375.1 is empty
GCA_008974285.1 vs GCA_014858625.1 is empty
GCA_008974285.1 vs GCA_018394375.1 is empty
GCA_008974285.1 vs GCA_022750515.1 is empty
GCA_008974285.1 vs GCF_000142945.1 is empty
GCA_014858625.1 vs GCA_018394375.1 is empty
GCA_018394375.1 vs GCA_022750515.1 is empty
GCA_018394375.1 vs GCF_000142945.1 is empty"""

def create_df(input):
    obj = {}
    for g in genomes:
        obj[g] = [0 for _ in range(len(genomes))]
    df = pd.DataFrame(data=obj)
    df.index = genomes    
    for i in input:
        df[i[0]][i[1]] = 1
        df[i[1]][i[0]] = 1
    return df

def split_input(input):
    return [(i.split()[0], i.split()[2]) for i in input.split('\n')]

t10 = create_df(split_input(seed10))
t20 = create_df(split_input(seed20))
t30 = create_df(split_input(seed30))
t40 = create_df(split_input(seed40))
t50 = create_df(split_input(seed50))

print(t10.to_markdown())
print(t20.to_markdown())
print(t30.to_markdown())
print(t40.to_markdown())
print(t50.to_markdown())
