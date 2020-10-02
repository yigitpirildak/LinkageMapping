import sys
from genetic_mapper import Mapper
import data_loader

if __name__ == "__main__":
    if (len(sys.argv) <= 1):
        print("Please enter path to frequency file!")
    else:
        freq_file = sys.argv[1]
        phenotypeFreq = data_loader.get_data(freq_file)
        print(phenotypeFreq)
        Mapper(phenotypeFreq).solve().print_solution()
