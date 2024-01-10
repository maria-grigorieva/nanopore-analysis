import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def calculate_probabilities(matrix):
    return pd.DataFrame({letter: np.sum(matrix == letter, axis=0) / matrix.shape[0]
                         for letter in ['A', 'C', 'G', 'T']}).T

def plot_probabilities(df):
    # Set the colors for each letter
    colors = {'A': 'red', 'C': 'green', 'G': 'blue', 'T': 'orange'}

    # Plot the probabilities for each letter
    for letter in df.index:
        plt.plot(df.loc[letter], label=letter, color=colors[letter])

    # Increase the width of the plot to 20
    fig = plt.gcf()
    fig.set_size_inches(20, fig.get_figheight())

    # Set the legend, x-label, and y-label
    plt.legend()
    plt.xlabel('Position')
    plt.ylabel('Probability')

    # Add labels with the titles of the letters with the highest probability
    for column in df.columns:
        max_letter, max_prob = df[column].idxmax(), df[column].max()
        plt.text(int(column), max_prob, max_letter, ha='center', va='bottom', fontsize=10)

    plt.savefig('ATCCATTCCCTCCGATAGATGAAACCAGCAC.png')
    plt.clf()
    plt.cla()
    plt.close()

# Read the file and convert it to a pandas DataFrame
data = []
with open('../input_data/ATCCATTCCCTCCGATAGATGAAACCAGCAC.seq', 'r') as file:
    for line in file:
        data.append(list(line.strip()))

df_seqs = pd.DataFrame(data)

df_prob = calculate_probabilities(df_seqs)

# Save the DataFrames to an Excel file
with pd.ExcelWriter('ATCCATTCCCTCCGATAGATGAAACCAGCAC.xlsx') as writer:
    df_prob.to_excel(writer, sheet_name='Frequencies')
    df_seqs.to_excel(writer, sheet_name='Sequences')

plot_probabilities(df_prob)
