import pandas as pd
import glob
import tqdm
maps = glob.glob('*.gmap')

for map in tqdm.tqdm(maps):
    map_df = pd.read_csv(map, sep='\t')
    map_df['cM'] = map_df['cM'].sort_values().values
    map_df['pos'] = map_df['pos'].sort_values().values
    map_df['CHM13v2.chr'] = (map_df['cM'].diff()/map_df.pos.diff()*1000000).fillna(0)
    map_df = map_df.rename(columns={'CHM13v2.chr':'cM/Mb'})
    map_df.to_csv(map + '.resorted.gmap.gz', sep='\t', index=False)