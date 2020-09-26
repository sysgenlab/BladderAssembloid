## gene annotation
import pandas as pd
def geneID2uniprot():
    output = {} # { gene ID : uniprot ID }
    df = pd.read_table('../data/UNIPROT_GENESYMBOL.tab')
    for i in range(len(df)):
        uniprot, geneList = df['Entry'][i], df['Gene names'][i]
        if pd.isnull(geneList) == False:
            geneList = geneList.split()
            for gene in geneList:
                output[gene.upper()] = uniprot
    return output

def uniprot2geneID():
    output = {} # { uniprot ID : gene ID }
    df = pd.read_table('../data/UNIPROT_GENESYMBOL.tab')
    for i in range(len(df)):
        uniprot, geneList = df['Entry'][i], df['Gene names'][i]
        if pd.isnull(geneList) == False:
            geneList = geneList.split()
            gene = geneList[0]
            output[uniprot] = gene.upper()
    return output
