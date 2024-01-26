from Bio import SeqIO
import sys
import argparse
from pathlib import Path
from PIL import Image, ImageDraw

def main(args):
    parser = argparse.ArgumentParser(description='Visualization of phred score for fastq files')
    args_info = [
        ['-i', '--input', str, 'Path to the fastq input file', 'input_file'],
        ['-l', '--limit', int, 'Limit number of sequences', 'limit', 1000],
        ['-sl', '--seq_length', int, 'Sequence length', None]
    ]

    for arg_info in args_info:
        parser.add_argument(arg_info[0], arg_info[1], type=arg_info[2], help=arg_info[3], default=arg_info[4])

    args = parser.parse_args()

    filename = args.input
    limit = args.limit
    seq_len = args.seq_length

    create_image_with_colored_sequence(filename, limit, seq_len=seq_len)


def create_image_with_colored_sequence(filename, limit=200, seq_len=None):
    records = list(SeqIO.parse(filename, "fastq"))[:limit]
    width = 20  # Width of each character box
    height = 20  # Height of each character box
    padding = 5  # Padding between character boxes
    max_score = 40  # Maximum score

    # Calculate the total width and height of the image
    total_width = limit * (width + padding)
    total_height = height * len(records)

    image = Image.new("RGB", (total_width, total_height), "white")
    draw = ImageDraw.Draw(image)

    y = 0

    for record in records:
        sequence = str(record.seq) if seq_len is None else str(record.seq)[:seq_len]
        scores = record.letter_annotations["phred_quality"] if seq_len is None else \
                 record.letter_annotations["phred_quality"][:seq_len]

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

    image.save(f'{Path(filename).stem}_phred_score_detailed.png')


if __name__ == "__main__":
    main(sys.argv[1:])