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

3. Optional arguments:
   - `-al, --alen`: Length of the aptamer (default: 31)
   - `-pl, --plen`: Length of the primer (default: 20)
   - `-i, --input`: Directory with input data (default: 'input_data')
   - `-rl, --reflen`: Initial reference length (default: 9)
   - `-r, --ref`: Initial reference sequence (default: 'auto')
   - `-p, --pos`: Start position of the reference sequence (default: -1)
   - `-s, --save`: Save results to excel (default: False)

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