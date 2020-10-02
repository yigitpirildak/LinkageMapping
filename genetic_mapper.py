from enum import Enum


class PhenotypeClass(Enum):
    UNKNOWN = 1
    PARENTAL = 2
    DCO = 3
    SCO = 4

class PhenotypeFrequency:

    def __init__(self, phenotype, frequency):
        self.phenotype = phenotype
        self.frequency = frequency
        self.phenotypeClass = PhenotypeClass.UNKNOWN

    def set_phenotype_class(self, phenotypeClass):
        self.phenotypeClass = phenotypeClass

    def __lt__(self, other):
        return self.frequency < other.frequency

    def __str__(self):
        return f"{self.phenotype} {self.frequency} {self.phenotypeClass}"

    def __add__(self, other):
        return self.phenotype + other.phenotype

    def difference(self, otherPhenotype):
        differenceList = []
        for characterIndex in range(0, len(self.phenotype)):
            if self.phenotype[characterIndex] != otherPhenotype.phenotype[characterIndex]:
                differenceList.append(self.phenotype[characterIndex].upper())
        return differenceList

class Mapper:

    MIN_SCO_INDEX = 2
    MAX_SCO_INDEX = 5

    def __init__(self, frequencyMap):
        self.frequencyList = []
        for key, value in frequencyMap.items():
            self.frequencyList.append(PhenotypeFrequency(key, value))
        self.frequencyList.sort()
        self.__check_gene_names_and_consistency()
        self.__determine_phenotype_classes()
        self.totalFrequency = sum(freqObj.frequency for freqObj in self.frequencyList)

    def solve(self):
        # Compare PARENTALS to DCO to determine the gene in the middle
        parentalDcoDifference = self.frequencyList[self.MIN_SCO_INDEX - 2].difference(self.frequencyList[self.MAX_SCO_INDEX + 2])
        if (len(parentalDcoDifference) == 1):
            self.__rearrange_gene_order(parentalDcoDifference[0])
        else:
            parentalDcoDifference = self.frequencyList[self.MIN_SCO_INDEX - 1].difference(self.frequencyList[self.MAX_SCO_INDEX + 2])
            self.__rearrange_gene_order(parentalDcoDifference[0])

        dcoFrequency = (self.frequencyList[self.MIN_SCO_INDEX - 2].frequency + self.frequencyList[self.MIN_SCO_INDEX - 1].frequency) / 1000
        diff1 = self.__calculate_distance((self.geneOrder[0], self.geneOrder[1]))
        diff2 = self.__calculate_distance((self.geneOrder[1], self.geneOrder[2]))
        diff3 = diff1 + diff2
        coc = (self.frequencyList[self.MIN_SCO_INDEX - 2].frequency + self.frequencyList[self.MIN_SCO_INDEX - 1].frequency) / (diff1 * diff2 * self.totalFrequency)

        return MappingSolution(self.frequencyList,
                              self.geneOrder,
                              (diff1, diff2, diff3),
                              (diff1, diff2, diff3 - 2 * dcoFrequency),
                              coc)

    def __calculate_distance(self, genes):
        recombinationFrequency = 0
        for i in range(0, self.MAX_SCO_INDEX + 1):
            parentalDiffs = (self.frequencyList[self.MAX_SCO_INDEX + 2].difference(self.frequencyList[i]), self.frequencyList[self.MAX_SCO_INDEX + 1].difference(self.frequencyList[i]))

            for diff in parentalDiffs:
                for gene in genes:
                    if len(diff) == 1 and gene.upper() in diff:
                        recombinationFrequency += self.frequencyList[i].frequency
                        break

        recombinationFrequency /= self.totalFrequency
        return recombinationFrequency

    """
        Since we are attempting to map three genes, we know that lowest two are DCO,
        highest two frequencies are parental and the rest are SCO.
    """
    def __determine_phenotype_classes(self):
        for i in range(0, len(self.frequencyList)):
            if i < self.MIN_SCO_INDEX:
                phenotypeClass = PhenotypeClass.DCO
            elif i <= self.MAX_SCO_INDEX:
                phenotypeClass = PhenotypeClass.SCO
            else:
                phenotypeClass = PhenotypeClass.PARENTAL

            self.frequencyList[i].set_phenotype_class(phenotypeClass)

    def __check_gene_names_and_consistency(self):
        genes = self.frequencyList[0].phenotype.upper()
        self.geneOrder = list(genes)

        for phenotypeFreq in self.frequencyList:
            if phenotypeFreq.phenotype.upper() != genes:
                raise Exception("Inconsistent gene types!")

    def __rearrange_gene_order(self, midGene):
        midGene = midGene.upper()
        while (self.geneOrder[1] != midGene):
            self.geneOrder.insert(0, self.geneOrder[len(self.geneOrder) - 1])
            del self.geneOrder[len(self.geneOrder) - 1]


class MappingSolution:
    def __init__(self, frequencyList, genes, recombinationFreqsNoDco, recombinationFreqsWithDco, coc):
        self.frequencyList = frequencyList
        self.genes = genes
        self.recombinationFreqsNoDco = recombinationFreqsNoDco
        self.recombinationFreqsWithDco = recombinationFreqsWithDco
        self.coc = coc

    # Solutions are converted to centiMorgans
    def print_solution(self):
        print("================")
        self.__print_frequency_list()
        print("================")
        print(f"Gene Order : {self.genes}")
        print(f"Coefficient of coincidence(coc) : {self.coc}")
        print(f"Interference : {1 - self.coc}")
        print("================")
        print("Gene Distances without counting double crossovers")
        self.__print_gene_distance((self.genes[0], self.genes[1]), self.recombinationFreqsNoDco[0])
        self.__print_gene_distance((self.genes[1], self.genes[2]), self.recombinationFreqsNoDco[1])
        self.__print_gene_distance((self.genes[0], self.genes[2]), self.recombinationFreqsNoDco[2])
        print("================")
        print("Gene Distances counting double crossovers")
        self.__print_gene_distance((self.genes[0], self.genes[1]), self.recombinationFreqsWithDco[0])
        self.__print_gene_distance((self.genes[1], self.genes[2]), self.recombinationFreqsWithDco[1])
        self.__print_gene_distance((self.genes[0], self.genes[2]), self.recombinationFreqsWithDco[2])
        print("================")

    def __print_gene_distance(self, genes, recombinationFreq):
        print(f"{genes[0]}-{genes[1]} : {recombinationFreq * 100}cM")

    def __print_frequency_list(self):
        for phenotypeFreq in self.frequencyList:
            print(phenotypeFreq)
