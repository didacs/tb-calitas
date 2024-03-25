#!/usr/bin/env python

import argparse
import os
import pandas as pd
import psycopg2
from dotenv import load_dotenv

parser = argparse.ArgumentParser(description='')
parser.add_argument('--spacers', help='TXT file containing a list of spacer ids (SPXXX), one per line, no header',
                    required=True)
parser.add_argument('--outfile', help='Output CSV file with spacer metadata', required=True)
args = parser.parse_args()

sp_list = tuple(pd.read_csv(args.spacers, header=None)[0])
if len(sp_list) == 1:
    sp_list = str(sp_list).replace(',','')

query = f"""
select 
    sp.file_registry_id$ as sp_id,
    tg.name$ as gene,
    sp.chromosome,
    sp.start,
    sp.end,
    sp.strand,
    sp.cut_position as cut_pos,
    dna.bases as bases,
    sp.Cas_PAM as PAM
from spacer$raw as sp
join dna_oligo$raw as dna on sp._sync_key = dna._sync_key
join target_gene$raw as tg on tg.id = sp.target_gene
where sp.is_registered$ = true and sp.archived$ = false and sp.file_registry_id$ in {sp_list}
"""

load_dotenv()
username = os.getenv('WAREHOUSE_USERNAME')
password = os.getenv('WAREHOUSE_PASSWORD')
url = os.getenv('WAREHOUSE_URL')
assert username and password and url, 'WAREHOUSE credentials are missing'
conn = psycopg2.connect(f"dbname=warehouse user={username} password={password} port=5432 host={url}")
cur = conn.cursor()
cur.execute(query)
sp = cur.fetchall()
names = [x[0] for x in cur.description]
df = pd.DataFrame(sp, columns=names)
sp_list = tuple(pd.read_csv(args.spacers, header=None)[0])
missing = [item for item in list(sp_list) if item not in list(df['sp_id'])]
assert (len(df) == len(sp_list)), f'Not all spacers were retrieved from benchling {missing}'
df.to_csv(args.outfile, index=0)




