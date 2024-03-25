#!/usr/bin/env python

import argparse
import re

import pandas as pd
import pybedtools


def parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--sites', required=True, help='calitas output file is the input')
    parser.add_argument('--gtf', required=True, help='gene set annotation in gtf format')
    parser.add_argument('--output', required=True, help='excel file with multiple tabs')
    return parser.parse_args()






def get_counts(x):
    if 'count' in x:
        return int(x.split('count:')[1])
    else:
        return '*'

def overlap_cosmic_cgc(df):
    # load cosmic genes
    # TODO: make sure genome build used was GRCh38/hg38
    url = 's3://tomebfx-data/COSMIC/v99/cosmic/GRCh38/cancer_gene_census.csv'
    cgc = pd.read_csv(url)
    cgc['strand'] = '.'
    cgc[['chrom', 'start', 'end']] = cgc['Genome Location'].str.split(':|-', expand=True)
    cgc['chrom'] = 'chr' + cgc['chrom']
    cols = ['chrom', 'start', 'end', 'Gene Symbol', 'Tier', 'strand']
    cgc_for_bed = cgc[cols]
    # some genes have not defined coordinates, cannot overlap coordinates here so let's drop them
    cgc_for_bed = cgc_for_bed.replace('', pd.NA)  # Convert empty strings to pd.NA (optional)
    cgc_for_bed = cgc_for_bed.dropna()
    # create bedtool object for CGC genes
    cgc_bed = pybedtools.BedTool.from_dataframe(cgc_for_bed)

    # create bedtool object for calitas sites
    cols = 'chromosome coordinate_start coordinate_end name score strand'.split()
    bed = pybedtools.BedTool.from_dataframe(df[cols])

    # intersect calitas sites CGC genes, type loj (left outer join)
    names = 'chromosome coordinate_start coordinate_end name score strand _chr _start _end gene_symbol cosmic_tier _strand'.split()
    cgc_sites = bed.intersect(cgc_bed, wa=True, wb=True, loj=True).to_dataframe(names=names)
    cgc_cols = 'name gene_symbol cosmic_tier'.split()
    df = df.merge(cgc_sites[cgc_cols], how="left")
    df['cosmic_tier'] = df['cosmic_tier'].astype(int)
    df['cosmic_tier'] = df['cosmic_tier'].replace(-1,'.')
    return df


def overlap_gene_set(df, gtf):
    """
    Computes the overlaps between calitas sites and the gene annotation provided as GTF.

    Returns the input dataframe with the following added columns:
        threat_tier: 4 tier category for cryptic sites based on the overlapping gene feature
            Tier I: coding regions only
            Tier II: non-coding regions (exonic UTR or intronic) of coding genes
            Tier III: exonic or intronic regions of non-coding genes
            Tier IV: all regions not in Tier I-III
        overlapping_gene: information about overlapping genes in the format of
            <gene_id>,<strand>,<gene_biotype>;
        overlapping_feature: specific gene feature overlapped by the cryptic site
        same_strand: whether the cryptic site and the gene are in the same strand,
            coulb be "both" if the site overlaps with gene on both strands

    Runs pybedtools.intersect
    """

    df['name'] = (df["chromosome"] +':'+
                  df['coordinate_start'].astype(str) +':'+
                  df['coordinate_end'].astype(str) +':'+
                  df['strand'])
    original_df = df.copy()

    # make sure name is unique
    duplicated_rows = df[df["name"].duplicated()]
    assert duplicated_rows.empty, f"Duplicated sites:\n{duplicated_rows}"

    cols = 'chromosome coordinate_start coordinate_end name score strand'.split()
    bed = pybedtools.BedTool.from_dataframe(df[cols])

    names = 'chromosome coordinate_start coordinate_end name score strand seqname source feature _start _end _score _strand frame group'.split()
    intersected_df = bed.sort().intersect(b=gtf, sorted=True, wb=True, loj=True).to_dataframe(names=names)
    # group by guide name and collapse columns feature and group as lists
    # each row now corresponds to a unique site, and contains the list of overlapping features
    by_fields = 'chromosome coordinate_start coordinate_end name strand'.split()
    grouped_data = intersected_df.groupby(by_fields).agg({'feature': list, 'group': list, '_strand': list}).reset_index()
    # make sure all sites are here
    assert len(bed) == len(grouped_data)

    # get thread tier for each site
    # the tier is based on the overlaps returned by BedTool.intersect. Briefly,
    # if a site overlaps with CDS it returns Tier I
    # else if overlaps with a coding gene, returns Tier II
    # else if overlaps with a gene whose biotype is other than protein_coding, returns Tier III
    # else if site overlaps no features in the gtf file, returns IV
    grouped_data['threat_tier'] = grouped_data.apply(get_threat_tier, axis=1)
    # get overlapping gene and the corresponding feature (CDS,intron,UTR,non_coding_exon)
    grouped_data['overlapping_gene'] = grouped_data.apply(overlapping_gene, axis=1)
    grouped_data['overlapping_feature'] = grouped_data.apply(overlapping_feature, axis=1)
    grouped_data = grouped_data.drop(columns=["feature", "group", "_strand"])

    # merge with original df
    df = pd.merge(original_df, grouped_data)
    return df

def get_threat_tier(row):
    """"""
    assert len(row.feature) == len(row.group)

    if row.feature == row.group == ['.']:
        return 'IV'
    elif 'CDS' in row.feature: # and site_overlaps_protein_coding(row):
        return 'I'
    elif site_overlaps_protein_coding(row):
        return 'II'
    elif site_overlaps_non_coding(row):
        return 'III'
    else:
        raise

def site_overlaps_protein_coding(row):
    """"""
    if not 'gene' in row.feature: raise
    for feature,group in zip(row.feature,row.group):
        if feature == 'gene':
            result = parse_group(group)
            if result['gene_biotype'] == 'protein_coding':
                return True

def site_overlaps_non_coding(row):
    """"""
    if not 'gene' in row.feature: raise
    biotypes = []
    for feature,group in zip(row.feature,row.group):
        if feature == 'gene':
            result = parse_group(group)
            biotypes.append(result['gene_biotype'])
    if not 'protein_coding' in biotypes:
        return True

def parse_group(group):
    """"""
    # Regex pattern to match key-value pairs
    pattern = r'(\w+)\s+"([^"]*)"\s*;'
    # Find all matches using the pattern
    matches = re.findall(pattern, group)
    # Create a dictionary of key-value pairs
    return dict(matches)

def overlapping_gene(row):
    if not 'gene' in row.feature:
        return '-'
    df = pd.DataFrame({'feature': row.feature, 'group': row.group, 'strand': row._strand})
    genes = []
    for i,d in df.iterrows():
        if d.feature == 'gene':
            result = parse_group(d.group)
            gene_id = result['gene_id']
            gene_biotype = result['gene_biotype']
            genes.append(f'{gene_id},{d.strand},{gene_biotype}')
    return '; '.join(genes)

def overlapping_feature(row):
    """"""
    if not 'gene' in row.feature:
        # no gene overlap
        return '-'
    elif not 'exon' in row.feature:
        # gene overlap but no exon, must be intron
        return 'intron'
    elif 'CDS' in row.feature:
        # CDS overlap
        return 'CDS'
    for feature,group in zip(row.feature,row.group):
        result = parse_group(group)
        if feature == 'exon' and result['transcript_biotype'] == "mRNA":
            # no CDS, exon overlap in protein_coding must be UTR
            return 'UTR'
    if 'exon' in row.feature:
        # exon in transcript_biotype other than mRNA
        return 'non_coding_exon'
    else:
        raise


if __name__ == "__main__":
    args = parse_args()
    df = pd.read_csv(args.sites, sep='\t')
    # add pam from target
    df['pam_target'] = df['padded_target'].str[-3:]
    # get threat tier and overlapping genes and features
    df = overlap_gene_set(df, gtf=args.gtf)
    # overlap with cosmic genes
    df = overlap_cosmic_cgc(df)
    first_cols = [
        'guide_id',
        'total_mm_plus_gaps',
        'overlapping_gene',
        'threat_tier',
        'cosmic_tier',
        'overlapping_feature',
        'pam_target'
    ]
    
    new_order = first_cols + [col for col in df.columns if col not in first_cols]
    df = df[new_order]
    df = df.sort_values('score', ascending=False)
    df.to_string(args.output, index=0)
    pybedtools.helpers.cleanup(verbose=False, remove_all=False)

