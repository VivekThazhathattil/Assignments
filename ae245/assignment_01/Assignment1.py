import math
import sys
import os
import matplotlib
import xml.etree.ElementTree as ET

# Let the H2-air reaction be:
# H2 + 1/2 (O2 + 3.76 N2) ----> (a H2O) + (b OH) + (c O2) + (d H2) + (e O) + (f H) + ( g N2)

species_list = ['H2O', 'OH', 'O2', 'H2', 'O', 'H','N2']

# Formation equilibria reactions:
reactions = [
        "H2 + 0.5 O2 <----> H2O",
        "0.5 O2 <----> O",
        "0.5 H2 <----> H",
        "0.5 O2 + 0.5 H2 <----> OH"
        ] 

# Parse BURCAT_truncated.xml file
# Note: a truncated version of the original BURCAT.xml is used.
# It has only the species relevant to the current problem.
def BURCAT_parser(sp_list: list, filename: str = "BURCAT_truncated.xml") -> dict:
    """
    input: * sp_list which contains all the species relevant to the problem,
           * filename, which can be explicitly passed. Else "BURCAT_trucated.xml" file is used.

    output: coeffs of temperature fit data required for Kp, in the form of a dict with keys being species fed
            and values being a list of two lists each - first for low T and second for high T
    """
    if not os.path.exists(filename):
        print("File: " + filename + " not found. Exiting...")
        sys.exit(0)

    tree = ET.parse(filename)
    root = tree.getroot()
    T_data = {}
    all_species = root.findall('specie')
    for sp in sp_list:
        T_range1 = []
        T_range2 = []
        T_overall = []
        for specie in all_species:
            if specie.find('formula_name_structure').find('formula_name_structure_1').text.split(' ')[0] == sp:
                coefs = specie.find('phase').find('coefficients').find('range_1000_to_Tmax').findall('coef')
                for coef in coefs:
                    T_range1.append(float(coef.text))
                coefs = specie.find('phase').find('coefficients').find('range_Tmin_to_1000').findall('coef')
                for coef in coefs:
                    T_range2.append(float(coef.text))
                T_overall.append(T_range1) 
                T_overall.append(T_range2)
                break
        T_data[sp] = T_overall
    return T_data

def reaction_to_dicts (reaction: str):
    """
    input: reaction string from the 'reactions' list

    output: two dicts, LHS and RHS, with keys being the species and their
    number of moles as the values
    """
    def helper_func(_list):
        _dict = {}
        _list = [species.strip() for species in  _list.split('+')]
        # extract the number of moles and species into the _dict
        for species in _list:
            temp_species = species.split(' ')
            if len(temp_species) == 1:
                _dict[ temp_species[0].strip() ] = 1.0
            else:
                _dict[ temp_species[1].strip() ] = float(temp_species[0].strip())
        return _dict

    temp_list = reaction.split('<---->')     
    if len(temp_list) != 2:
        print("Error in processing equilibria reaction string: " + reaction + ".\n Exiting...")
        sys.exit(0)

    LHS = helper_func(temp_list[0])
    RHS = helper_func(temp_list[1])

    return LHS, RHS


def get_species_with_given_atom (atom: str, species_list: list) -> dict:
    """ 
    input: atom symbol. Eg: 'H', 'O', 'N'

    output: 'species_dict' dictionary whose keys are species containing 
    the given atom and key is the number of given atoms for that species

    """

    if len(atom) > 2 or atom[0].islower() or not atom.isalpha():
        print("Invalid atomic symbol. Exiting...")
        sys.exit(0)

    species_dict = {}

    # caveat: 1) cannot accurately count the number of some atoms in some species.
    # Eg: H atoms in CH3COOH species 
    # Here n(H) = 3, which should instead be 4. Similary n(C) = 1, but instead 
    # should be 2
    # 2) cannot specify phases for the same species

    for species in species_list:
        test_result = species.find(atom)
        next_index = len(atom) + test_result
        if test_result >= 0: # match found
            if next_index == len(species): # next_index exceeds last index
                species_dict[species] = 1.0
            elif next_index  == len(species) - 1: # next_index == last index
                if species[next_index].isnumeric():
                    species_dict[species] = float(species[next_index])
                else:
                    species_dict[species] = 1.0
            elif len(atom) + test_result < len(species) - 1 and len(atom) >= 0:
                if species[next_index].isnumeric():
                    species_dict[species] = float(species[next_index])
                else:
                    species_dict[species] = 1.0
            else:       # i.e for eg, next_index > len(species):
                print("Error in finding number of atoms. Exiting...")
                sys.exit(0)
    return species_dict

def get_atoms_for_species( species: str ) -> list:
    """ 
    input: species. Eg: 'CH4', 'H2O', 'HNO3'
    output: atoms present in the input species 
        CH4 -> ['C', 'H']
        H2O -> ['H', 'O']
        HNO3 -> ['H', 'N', 'O']
    """
    atoms = []
    for i in range(len(species)):
        if species[i].isupper():
            if i + 1 == len(species):
                if species[i] not in atoms:
                    atoms.append(species[i])
                return atoms
            elif species[i+1].islower() and "" + species[i] + species[i+1] not in atoms:
                atoms.append("" + species[i] + species[i+1])
            elif species[i] not in atoms:
                atoms.append(species[i])
    return atoms

#def get_Kp( reaction: str) -> float:
#    """
#    input: equilibrium reaction string
#
#    output: Kp for that reaction
#    """
#    LHS, RHS = reaction_to_dicts(reaction)

#### TEST FUNCTIONS ##################

## checks for get_species_with_given_atoms()
#print(get_species_with_given_atom( 'H', ['H2O', 'OH', 'O2', 'H2', 'O', 'H', 'N2']))
#print(get_species_with_given_atom( 'O', ['H2O', 'OH', 'O2', 'H2', 'O', 'H', 'N2']))
#print(get_species_with_given_atom( 'N', ['H2O', 'OH', 'O2', 'H2', 'O', 'H', 'N2']))
#print(get_species_with_given_atom( 'He', ['ABCHe', 'ABCHe2', 'He2ABC', 'HeABC', 'AHeBC', 'AHe2BC', 'ABC']))

## checks for get_atoms_for_species()
#print(get_atoms_for_species('H2SO4'))
#print(get_atoms_for_species('CH3COOH'))

## checks for reaction_to_dicts(reaction) function
#print(reaction_to_dicts(reactions[0]))
#print(reaction_to_dicts(reactions[1]))
#print(reaction_to_dicts(reactions[2]))
#print(reaction_to_dicts(reactions[3]))

## checks for BURCAT_parser()
#burcat_parsed = BURCAT_parser(species_list)
#for sp in burcat_parsed:
#    print(sp)
#    for T_dat in burcat_parsed[sp]:
#        print(T_dat)
#    print("\n")
