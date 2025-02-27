# nanopore-analysis

## Aptamer Search

This code searches for an aptamer among low-quality sequences with a known length and primers. It utilizes probabilistic analysis and slicing window techniques to identify potential aptamer sequences. The code takes in input data containing the DNA sequences and outputs the results, including plots and excel files with the analysis steps.

### Installation

To run this code, you need to have the following dependencies installed:

- pandas
- matplotlib
- numpy
- Biopython

You can install these dependencies by running the following command:

```
pip install pandas matplotlib numpy biopython
```

### Usage

1. Prepare the input data:
   - Create a directory and place the input data containing the DNA sequences in it.
   - Each sequence should be stored as a separate file.
   
2. Run the code:
   - Open a terminal and navigate to the directory where the code is saved.
   - Execute the following command:
   
     ```
     python aptamer_search.py -i <input_directory>
     ```
     
     Replace `<input_directory>` with the path to the directory containing the input data.

3. Arguments:
   - `-al, --alen`, int, 'Length of an aptamer', (default: 31)
   - `-i, --input`, str, 'Path to the input fastq file', (default: 'input_data')
   - `-o, --output`, str, 'Directory with output data', (default :'../results')
   - `-pl, --left_primer`, str, 'Left Primer'
   - `-pr, --right_primer`, str, 'Right Primer'
   - `-r, --ref`, str, 'Initial reference sequence'
   - `-p, --pos`, int, 'Start position of the reference sequence', (default: -1)
   - `-f, --fuzzy`, bool, 'Add fuzzy search', (default: False)
   - `-s, --save`, bool, 'Save to excel (True/False)', (default: False)
   - `-c, --complement`, bool, 'Add complementary primer', (default: False)
   - `-w, --weights`, bool, 'Take into account nucleotide phred scores', (default: False)
   - `-ph, --cutoff`, int, 'Phred score cut off', (default: 15)
   - `-flex, --flexible_shifts`, bool, 'Flexible slicing window (improves performance, but decreases accuracy)'

4. Output:
   - The code will output information about the analysis steps and the selected aptamer candidates.
   - It will also generate plots with frequencies and excel files containing the analysis steps (if `-s` argument is set to `True`).
   - The results will be saved in the `output` directory within the input directory.
   - The plots will be saved in the `plots/references` subdirectory of the output directory.

## Contributing

If you want to contribute to this project, feel free to submit a pull request with your suggestions or improvements. Please make sure to follow the existing code style and guidelines.

## License

This code is licensed under the MIT License. You can find the details in the [LICENSE](LICENSE) file.

## Author

This code was developed by magsend@gmail.com and contributors.

For any questions or inquiries, please contact [Your Email].