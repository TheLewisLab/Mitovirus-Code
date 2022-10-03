import getopt
import sys
import pandas as pd
from Bio import SeqIO
from Bio import SeqUtils
import numpy as np
import os

def fileIO(argv):
    inputfile = ''
    trans_table =''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:t:o:", ["ifile=", "table=", "ofile="])
    except getopt.GetoptError:
        print('CodonFrequency.py -i <inputfile> -t <translationtable> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('CodonFrequency.py -i <inputfile> -t <translationtable> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-t", "--trans"):
            trans_table = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    print('Input file is "', inputfile)
    print('Output file is "', outputfile)
    return inputfile, trans_table, outputfile


def main(argv):
    inputFile, trans_table, outputFile = fileIO(argv)
    

    frequency = {
    'TTT' : 0, 'TCT' : 0, 'TAT' : 0, 'TGT': 0,
    'TTC' : 0, 'TCC' : 0, 'TAC' : 0, 'TGC' : 0,
    'TTA' : 0, 'TCA' : 0, 'TAA' : 0, 'TGA' : 0,
    'TTG' : 0, 'TCG' : 0, 'TAG' : 0, 'TGG' : 0,

    'CTT' : 0, 'CCT' : 0,  'CAT' : 0, 'CGT' : 0,
    'CTC' : 0, 'CCC' : 0, 'CAC' : 0, 'CGC' : 0,
    'CTA' : 0, 'CCA' : 0, 'CAA' : 0, 'CGA' : 0,
    'CTG' : 0, 'CCG' : 0, 'CAG' : 0, 'CGG' : 0,

    'ATT' : 0, 'ACT' : 0, 'AAT' : 0, 'AGT' : 0,
    'ATC' : 0, 'ACC' : 0, 'AAC' : 0, 'AGC' : 0,
    'ATA' : 0, 'ACA' : 0, 'AAA' : 0, 'AGA' : 0,
    'ATG' : 0, 'ACG' : 0, 'AAG' : 0, 'AGG' : 0,

    'GTT' : 0, 'GCT' : 0, 'GAT' : 0, 'GGT' : 0,
    'GTC' : 0, 'GCC' : 0, 'GAC' : 0, 'GGC' : 0,
    'GTA' : 0, 'GCA' : 0, 'GAA' : 0, 'GGA' : 0,
    'GTG' : 0, 'GCG' : 0, 'GAG' : 0, 'GGG' : 0
    }
    
    dirname = os.path.dirname(__file__)

    fungi_mito_df = pd.read_csv(dirname + '/codon_references/fungi_mito_codon_frequency.csv')
    mammalia_mito_df = pd.read_csv(dirname + '/codon_references/mammalia_mito_codon_frequency.csv')
    viridiplantae_mito_df = pd.read_csv(dirname + '/codon_references/viridiplantae_mito_codon_frequency.csv')
    metazoan_mito_df = pd.read_csv(dirname + '/codon_references/metazoan_mito_codon_frequency.csv')

    fungi_genome_df = pd.read_csv(dirname + '/codon_references/fungi_genome_codon_frequency.csv')
    mammalia_genome_df = pd.read_csv(dirname + '/codon_references/mammalia_genome_codon_frequency.csv')
    viridiplantae_genome_df = pd.read_csv(dirname + '/codon_references/viridiplantae_genome_codon_frequency.csv')
    metazoan_genome_df = pd.read_csv(dirname + '/codon_references/metazoan_genome_codon_frequency.csv')
    bacteria_df = pd.read_csv(dirname + '/codon_references/bacteria_codon_frequency.csv')
    yeast_df = pd.read_csv(dirname + '/codon_references/S_Cerevisiae_frequency.csv')
    
    sra_name = str(os.getcwd()).split('/')[-5]
    
    with open('summary.csv', 'w') as f:
        f.write('sra,node,codon_table,sequence,cds,protein,length,start,stop,fungi_mito_r2,fungi_genome_r2,mammalia_mito_r2,mammalia_genome_r2,'
                'viridiplantae_mito_r2,viridiplantae_genome_r2,metazoan_mito_r2,metazoan_genome_r2,bacteria_r2,yeast_r2,at_content\n')
        
        node_name = 'NA'
        sequence_path = 'NA'

        for item in os.listdir('./'):
            if str(item)[:4] != 'node' and str(item).split('.')[-1] == 'fa':
                sequence_path = str(item)
                temp_list = str(item).split('.')
                node_name = temp_list[0]
                print(node_name)
        
        full_sequence = SeqIO.read(sequence_path, 'fasta')
        full_sequence = full_sequence.seq
        at_content = 100 - SeqUtils.GC(full_sequence)
        
        for cdc in SeqIO.parse(inputFile, 'fasta'):
            sequ = str(cdc.seq)
            ambig_counter = 0

            for i in range(0, len(sequ), 3):
                codon = sequ[i:i + 3]
                if len(codon) > 0 and len(codon) < 3:
                    codon = codon + 'A' * (3 - len(codon))

                if frequency.get(codon) == None:
                    ambig_counter += 1
                else:
                    frequency[codon] = 1 + frequency.get(codon)

            cdc_df = pd.DataFrame.from_dict(frequency, orient='index', columns=['Count'])
            cdc_df['Frequency'] = 1000 * (cdc_df['Count'] / sum(frequency.values()))
            cdc_df = cdc_df.sort_index()
            cdc_df.to_csv(outputFile)
            gene = cdc.seq.translate(table=trans_table)
            cdc_range = str(cdc.id).split(':')[1]
            start = cdc_range.split('-')[0]
            end = cdc_range.split('-')[1]
            fungi_mito_r2 = np.corrcoef(cdc_df['Frequency'], fungi_mito_df['Frequency'])[0,1]**2
            fungi_genome_r2 = np.corrcoef(cdc_df['Frequency'], fungi_genome_df['Frequency'])[0,1]**2
            mammalia_mito_r2 = np.corrcoef(cdc_df['Frequency'], mammalia_mito_df['Frequency'])[0,1]**2
            mammalia_genome_r2 = np.corrcoef(cdc_df['Frequency'], mammalia_genome_df['Frequency'])[0,1]**2
            viridiplantae_mito_r2 = np.corrcoef(cdc_df['Frequency'], viridiplantae_mito_df['Frequency'])[0,1]**2
            viridiplantae_genome_r2 = np.corrcoef(cdc_df['Frequency'], viridiplantae_genome_df['Frequency'])[0,1]**2
            metazoan_mito_r2 = np.corrcoef(cdc_df['Frequency'], metazoan_mito_df['Frequency'])[0,1]**2
            metazoan_genome_r2 = np.corrcoef(cdc_df['Frequency'], metazoan_genome_df['Frequency'])[0,1]**2
            bacteria_r2 = np.corrcoef(cdc_df['Frequency'], bacteria_df['Frequency'])[0,1]**2
            yeast_r2 = np.corrcoef(cdc_df['Frequency'], yeast_df['Frequency'])[0,1]**2
            f.write('%s,%s,%d,%s,%s,%s,%d,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n' % (sra_name, node_name, int(trans_table), full_sequence, cdc.seq, gene, len(gene), start, end, fungi_mito_r2, fungi_genome_r2, mammalia_mito_r2, mammalia_genome_r2, viridiplantae_mito_r2, viridiplantae_genome_r2, metazoan_mito_r2, metazoan_genome_r2, bacteria_r2, yeast_r2, at_content))

if __name__ == "__main__":
    main(sys.argv[1:])
