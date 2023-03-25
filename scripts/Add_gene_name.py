import pandas as pd

df = pd.read_csv('featurecounts_comb.txt', sep = '\t', comment='#')


df.columns = df.columns.str.replace('/storage/marwe/LisaRNA_DNA/scripts/nextflow_RNA/../../data/RNA/Bowtie_samtools/', '')
df.columns = df.columns.str.replace('_R1_001.fastq.sorted.bam', '')

# Creating an empty list to store gene data
gene_data = []

# Path to the gtf file
gtf = 'ncbi_dataset/data/GCF_000005845.2/genomic.gtf'


with open(gtf, 'r') as f:
    for line in f:
        if 'GeneID' in line:
            # Extract the gene_id and gene_name information from the line
            gene_id = line.split('gene_id "')[-1].split('";')[0]
            gene_name = line.split('; gene "')[-1].split('";')[0]
            # Add the extracted data to the gene_data list
            gene_data.append((gene_id, gene_name))

# Convert the gene_data list to a pandas dataframe and set the gene_id as the index
gene_dict = pd.DataFrame(gene_data, columns=["gene_id", "gene_name"]).set_index("gene_id")

# Convert the gene_dict dataframe to a dictionary, with gene_id as the key and gene_name as the value
gene_dict = gene_dict.to_dict()["gene_name"]

# Map the gene_id in the input dataframe to the corresponding gene_name using the gene_dict dictionary
df["gene_name"] = df["Geneid"].map(gene_dict)


df.to_csv('featurecounts_comb.txt', sep = '\t', index=False)