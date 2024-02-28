import matplotlib.pyplot as plt
from PIL import Image, ImageDraw
def plot_probabilities(reference, plots_dir, title=None):
    # Set the colors for each letter
    colors = {'A': 'red', 'C': 'green', 'G': 'blue', 'T': 'orange'}

    # Plot the probabilities for each letter
    for letter in reference['freq'].index:
        plt.plot(reference['freq'].loc[letter], label=letter, color=colors[letter])

    # Increase the width of the plot to 20
    fig = plt.gcf()
    fig.set_size_inches(20, fig.get_figheight())

    # Set the legend, x-label, and y-label
    plt.legend()
    plt.xlabel('Position')
    plt.ylabel('Probability')

    # Add labels with the titles of the letters with the highest probability
    for column in reference['freq'].columns:
        max_letter, max_prob = reference['freq'][column].idxmax(), reference['freq'][column].max()
        plt.text(int(column), max_prob, max_letter, ha='center', va='bottom', fontsize=10)

    ref_seq = reference['seq']
    ref_pos = reference['start_pos']
    figname = f'{plots_dir}/{ref_pos}:{ref_seq}.png' if title is None else f'{plots_dir}/{ref_seq}-{title}.png'
    plt.savefig(figname)
    plt.clf()
    plt.cla()
    plt.close()

def create_image_with_colored_sequence(records, output_file, limit=200):
    width = 20  # Width of each character box
    height = 20  # Height of each character box
    padding = 5  # Padding between character boxes
    max_score = 40  # Maximum score

    total_width = limit * (width + padding)
    total_height = height * len(records)

    image = Image.new("RGB", (total_width, total_height), "white")
    draw = ImageDraw.Draw(image)

    y = 0

    for record in records:
        sequence = record['sequence']
        scores = record['score']

        x = 0

        for i, score in enumerate(scores):
            character = sequence[i]

            # Calculate the color based on the score
            normalized_score = score / max_score  # Normalize the score between 0 and 1
            red = int(255 * (1 - normalized_score))
            green = int(255 * normalized_score)
            color = (red, green, 0)

            # Draw the character box with the corresponding color
            draw.rectangle([x, y, x + width, y + height], fill=color)
            draw.text((x+8, y), character, fill="black")

            x += width + padding

        y += height

    image.save(output_file)