#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def read_gwas_file(file_path):
    # Load the file, inferring column names
    df = pd.read_csv(file_path, delim_whitespace=True, header=0)

    # Identify the chromosome column
    possible_chr_cols = ['#CHROM', 'CHR']
    chr_col = next((col for col in possible_chr_cols if col in df.columns), None)

    # Identify the position column
    possible_pos_cols = ['BP', 'POS']
    pos_col = next((col for col in possible_pos_cols if col in df.columns), None)

    # Identify the p-value column
    possible_pval_cols = ['P', 'p', 'PVAL', 'pval']
    pval_col = next((col for col in possible_pval_cols if col in df.columns), None)

    # Identify the rsID column
    possible_snp_cols = ['SNP', 'ID', 'rsID']
    snp_col = next((col for col in possible_snp_cols if col in df.columns), None)

    # Check if required columns are present
    if not chr_col or not pos_col or not pval_col or not snp_col:
        raise ValueError("Required columns ('CHR'/'#CHROM', 'BP'/'POS', 'P'/'p', or 'SNP'/'ID'/'rsID') not found. Check the file format.")

    return df[[chr_col, pos_col, snp_col, pval_col]].rename(columns={chr_col: 'CHR', pos_col: 'BP', snp_col: 'SNP', pval_col: 'P'})

def plot_manhattan(df, output_file):
    df['-log10(P)'] = -np.log10(df['P'])
    
    # Remove infinite values
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['-log10(P)'])

    # Manhattan plot setup
    df['CHR'] = df['CHR'].astype('category')
    df['ind'] = range(len(df))
    df_grouped = df.groupby(('CHR'))
    
    # Set the figure size wider for better spacing of x-axis ticks
    fig = plt.figure(figsize=(20, 6))
    ax = fig.add_subplot(111)

    colors = ['b', 'r']
    x_labels = []
    x_labels_pos = []

    # Variables to keep track of label positions for avoiding overlaps
    label_positions = []

    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='-log10(P)', color=colors[num % len(colors)], ax=ax, s=10, alpha=0.75)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] + group['ind'].iloc[0]) / 2)

        # Annotate SNPs with P < 1e-4
        significant_snps = group[group['P'] < 1e-4]
        for _, row in significant_snps.iterrows():
            x, y = row['ind'], row['-log10(P)']
            label = row['SNP']

            # Adjust label position to avoid overlap
            adjusted_y = y
            while any(abs(adjusted_y - prev_y) < 0.2 and abs(x - prev_x) < 100 for prev_x, prev_y in label_positions):
                adjusted_y += 0.2  # Shift the label upwards

            # Add label to the plot
            ax.text(x, adjusted_y, label, fontsize=4, rotation=0, ha='right', color='black')

            # Store the label position
            label_positions.append((x, adjusted_y))

    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('Manhattan Plot')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)  # Save plot at 300 dpi
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate a Manhattan plot from PLINK GWAS results.")
    parser.add_argument('-i', '--input', required=True, help="Input file containing GWAS results")
    
    args = parser.parse_args()

    # Check if the input file exists
    if not os.path.isfile(args.input):
        raise FileNotFoundError(f"Input file {args.input} not found.")

    # Generate output file name based on input file name
    input_file_name = os.path.basename(args.input)
    output_file_name = os.path.splitext(input_file_name)[0] + "_manhattan.png"

    # Read and plot the data
    df = read_gwas_file(args.input)
    plot_manhattan(df, output_file_name)
    print(f"Manhattan plot saved to {output_file_name}")

if __name__ == "__main__":
    main()

