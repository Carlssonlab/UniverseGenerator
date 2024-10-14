/*
    Author: Andreas Luttens
    Contact: andreas.luttens@gmail.com
    Date: June 22, 2023
    Description: Write instructions to build superstructures
*/

// Basic utilities
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>
#include <tuple>
#include <set>
#include <utility>
#include <iterator>

// RDKit functions
#include <GraphMol/GraphMol.h>
#include <GraphMol/Canon.h> // Check for symmetry later
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <GraphMol/MolStandardize/Normalize.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>

// A map between bondorders and bondtypes
std::unordered_map<unsigned int, RDKit::Bond::BondType> bond_types = {{1u, RDKit::Bond::SINGLE},
                                                                      {2u, RDKit::Bond::DOUBLE},
                                                                      {3u, RDKit::Bond::TRIPLE}};

// A map between bondtypes and bondorders
std::unordered_map<RDKit::Bond::BondType, unsigned int> rev_bond_types = {{RDKit::Bond::SINGLE, 1u},
                                                                          {RDKit::Bond::DOUBLE, 2u},
                                                                          {RDKit::Bond::TRIPLE, 3u}};

// A map between couplings and bondtypes
std::unordered_map<std::string, RDKit::Bond::BondType> bond_types_str = {{"single", RDKit::Bond::SINGLE},
                                                                         {"double", RDKit::Bond::DOUBLE},
                                                                         {"triple", RDKit::Bond::TRIPLE},
                                                                         {"single_fusion", RDKit::Bond::SINGLE}};

// A map between activation tag integers and their type of coupling
std::unordered_map<unsigned int,std::string> activation_types = {{71u, "single"},
                                                                 {72u, "double"},
                                                                 {73u, "triple"},
                                                                 {74u, "spiro"},
                                                                 {75u, "single_fusion"},
                                                                 {76u, "double_fusion"}};

//
void generateCombinationsHelper(const std::vector<int>& nums, int k, int start,\
                                std::vector<int>& currentCombination,\
                                std::vector<std::vector<int>>& combinations)
{
    if (static_cast<int>(currentCombination.size()) == k)
    {
        combinations.push_back(currentCombination);
        return;
    }

    for (auto it = nums.begin(); it != nums.end(); ++it)
    {
        int num = *it;
        if (num >= start)
        {
            currentCombination.push_back(num);
            generateCombinationsHelper(nums, k, num + 1, currentCombination, combinations);
            currentCombination.pop_back();
        }
    }
}

//
std::vector<std::vector<int>> generateCombinations(const std::vector<int>& nums, int k)
{
    std::vector<std::vector<int>> combinations;
    std::vector<int> currentCombination;
    generateCombinationsHelper(nums, k, INT_MIN, currentCombination, combinations);
    return combinations;
}

// Function to retrieve atom adjacencies
std::unordered_map<unsigned int, std::unordered_set<unsigned int>> get_atom_adjacencies(const RDKit::ROMol& mol)
{
    // Initialize an empty map to store the atom adjacencies
    std::unordered_map<unsigned int, std::unordered_set<unsigned int>> atom_adjacencies;

    // Iterate over all atoms in the molecule
    const unsigned int num_atoms = mol.getNumAtoms();
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx)
    {
        // Retrieve the atom with this index
        const RDKit::Atom *atom = mol.getAtomWithIdx(atom_idx);

        // Initialize an empty set to store indices of adjacent atoms of this atom
        atom_adjacencies.emplace(atom_idx,std::unordered_set<unsigned int>{});

        // Iterate over all neighboring atoms
        for(const auto &neighbor_atom: mol.atomNeighbors(atom))
        {
            // Get the index of this current atom
            const unsigned int neighbor_atom_idx = neighbor_atom->getIdx();

            // Insert the atom's index in the array of indices
            atom_adjacencies[atom_idx].emplace(neighbor_atom_idx);
        }
    }

    return atom_adjacencies;
}

// Function to get information on bonds given an RDKit molecule
std::unordered_map<unsigned int, std::pair<unsigned int,unsigned int>> get_bond_information(const RDKit::ROMol& mol)
{
    // Initialize an empty dictionary to store the bond information
    std::unordered_map<unsigned int, std::pair<unsigned int,unsigned int>> bond_information;

    // Iterate over all bonds in the molecule
    unsigned int num_bonds = mol.getNumBonds();
    for (unsigned int bond_idx = 0 ; bond_idx < num_bonds ; ++bond_idx)
    {
        // Retrieve the bond with this index
        const RDKit::Bond *bond = mol.getBondWithIdx(bond_idx);

        // Get the indices of the atoms forming this bond
        const unsigned int bond_bgn_idx = bond->getBeginAtomIdx();
        const unsigned int bond_end_idx = bond->getEndAtomIdx();

        // Create pair of atom indices
        std::pair<unsigned int,unsigned int> atom_pair = std::make_pair(bond_bgn_idx, bond_end_idx);

        // Store key value (bond_idx, atom_idx pair) pairs
        bond_information.emplace(bond_idx, atom_pair);
    }

    return bond_information;
}

// Function to retrieve bond adjacencies given an RDKit molecule
std::unordered_map<unsigned int, std::unordered_set<unsigned int>> get_bond_adjacencies(const RDKit::ROMol& mol)
{
    // initialize an empty dictionary to store the bond adjacencies
    std::unordered_map<unsigned int, std::unordered_set<unsigned int>> bond_adjacencies;

    // Iterate over all bonds in the molecule
    unsigned int num_bonds = mol.getNumBonds();
    for (unsigned int bond_idx = 0 ; bond_idx < num_bonds ; ++bond_idx)
    {
        // Retrieve the bond with this index
        const RDKit::Bond *bond = mol.getBondWithIdx(bond_idx);

        // Get the atoms forming the bond
        const RDKit::Atom *bond_bgn = bond->getBeginAtom();
        const RDKit::Atom *bond_end = bond->getEndAtom();

        // Initialize an array for the bond index to store indices of adjacent bonds
        bond_adjacencies.emplace(bond_idx,std::unordered_set<unsigned int>{});

        // Iterate over all bonds through bond begin atom
        for(const auto &neighbor_bond: mol.atomBonds(bond_bgn))
        {
            // Get the index of this neighbor bond
            unsigned int neighbor_bond_idx = neighbor_bond->getIdx();

            // insert the bond's index in the array of indices
            bond_adjacencies[bond_idx].emplace(neighbor_bond_idx);
        }

        // Iterate over all bonds through bond end atom
        for(const auto &neighbor_bond: mol.atomBonds(bond_end))
        {
            // Get the index of this neighbor bond
            unsigned int neighbor_bond_idx = neighbor_bond->getIdx();

            // insert the bond's index in the array of indices
            bond_adjacencies[bond_idx].emplace(neighbor_bond_idx);
        }
    }

    return bond_adjacencies;
}

// Function to get CIP ranks of each atom in a given RDKit molecule
std::unordered_map<unsigned int,unsigned int> get_atom_symmetries(const RDKit::ROMol& mol)
{
    // Initialize a map between the atom index and its CIP rank
    std::unordered_map<unsigned int,unsigned int> symmetry_map = {};

    // Iterate over all atoms in the molecule
    unsigned int num_atoms = mol.getNumAtoms();
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx) {
        
        // Get the atom with this index
        const RDKit::Atom *atom = mol.getAtomWithIdx(atom_idx);
        
        // Find the atom's CIP rank
        unsigned int rank;
        atom->getProp("_CIPRank", rank);

        // Store the atom and its CIP rank in the map
        symmetry_map.emplace(atom->getIdx(), rank);
    }

    return symmetry_map;
}

// Function to get the symmetry factor of each atom given their CIP ranks
std::unordered_map<unsigned int,unsigned int> get_symmetry_factors(std::unordered_map<unsigned int,unsigned int>& atom_symmetries)
{
    // Initialize a map between between the atom and its symmetry factor (sigma)
    std::unordered_map<unsigned int,unsigned int> symmetry_factors = {};

    // A vector to store all CIP ranks seen in a molecule
    std::vector<unsigned int> symmetry_classes = {};

    // Iterate over all atoms and their symmetries
    for (const auto& atom_sym_pair : atom_symmetries) {
        
        // Store the CIP rank
        symmetry_classes.emplace_back(atom_sym_pair.second);
    }

    // Iterate over all atoms and their symmetries
    for (const auto& atom_sym_pair : atom_symmetries) {
        
        // Store the atom and how many times its CIP rank was observed
        symmetry_factors.emplace(atom_sym_pair.first, count(symmetry_classes.begin(), symmetry_classes.end(), atom_sym_pair.second));
    }
    return symmetry_factors;
}

// Function to get all CIP ranks in a molecule and the atoms that have those CIP ranks
std::unordered_map<unsigned int,std::vector<unsigned int>> get_symmetry_classes(std::unordered_map<unsigned int,unsigned int>& atom_symmetries)
{
    // Initialize a map between CIP ranks and atom indices that have this CIP rank
    std::unordered_map<unsigned int,std::vector<unsigned int>> symmetry_classes = {};

    // Iterate over all atoms and their symmetries
    for (const auto& atom_sym_pair : atom_symmetries) {
        
        // If the CIP rank is not seen
        if (symmetry_classes.count(atom_sym_pair.second) == 0)
        {
            // Store it and the corresponding atom that has this rank
            symmetry_classes.emplace(atom_sym_pair.second, std::vector<unsigned int>{atom_sym_pair.first});
        }
        else
        {
            // Store the new atom that also has this rank
            symmetry_classes[atom_sym_pair.second].emplace_back(atom_sym_pair.first);
        }
    }
    return symmetry_classes;
}

// Generates all partitions of an integer n into k parts.
std::vector<std::vector<int>> generate_partitions(int& n, int& k)
{    
    // n: The integer to partition.
    // k: The number of parts in the partition.

    // Initialize vector for storing valid partitions
    std::vector<std::vector<int>> partitions = {};
    
    // Create valid range
    std::vector<int> range = {};
    for (int i = 0; i <= n + k - 1; i++)
    {
        range.emplace_back(i);
    }
    
    // Generate all possible combinations of k-1 separators
    std::vector<std::vector<int>> separator_positions = generateCombinations(range, k - 1);

    // Iterate over separator positions
    for (const std::vector<int>& separator : separator_positions)
    {    
        // Initialize the partition for this separator
        std::vector<int> partition = {};

        // Initialize vector for distance calculations
        std::vector<int> head = {-1};
        std::vector<int> tail = {n + k - 1};

        // Counter for the separator size
        // int separator_length = 1;

        // insert the separator sequence in the head
        for (const int& i : separator)
        {
            head.emplace_back(i);
            //++separator_length; //TODO I'm not increasing the separator length, but still works?
        }

        // insert the reverse separator sequence in the tail
        for (auto it = separator.rbegin(); it != separator.rend(); ++it)
        {
            int element = *it;
            tail.insert(tail.begin(), element);
        }

        // Compute distances between adjacent separators, this is the partition
        for (int i = 0; i < k; i++)
        {
            partition.emplace_back(tail[i] - head[i] - 1);
        }

        bool valid_partition = true;
        int sum_partition = 0;
        for (const int& i : partition)
        {
            if (i <= 0)
            {
                valid_partition = false;
                break;
            }
            sum_partition+=i;
        }
        if (valid_partition && sum_partition == n)
        {
            partitions.emplace_back(partition);
        }
    }

    return partitions;
}

// Helper function to do stepwise sorting
bool stepwise_sorter(const std::pair<int, int>& pair1, const std::pair<int, int>& pair2) 
{   

    if (pair1.first != pair2.first)
    // Sort by the first item
    {
        return pair1.first < pair2.first;
    } 
    else
    // If the first items are equal, sort by the second item
    {
        return pair1.second < pair2.second;
    }
}

// Check if the parsed atom is bound to an activator tag
bool is_atom_bound_to_activator(const RDKit::ROMol &mol, const RDKit::Atom &atom)
{
    // Element numbers to exclude as regular heavy atoms
    std::unordered_set<unsigned int> exclude_elem = {71u,72u,73u,74u};

    // Iterate over all neighboring atoms
    for(const auto &neighbor_atom: mol.atomNeighbors(&atom))
    {
        //
        const unsigned int neighbor_atom_elem = neighbor_atom->getAtomicNum();

        //
        if (exclude_elem.find(neighbor_atom_elem) != exclude_elem.end())
        {
            return true;
        }
        
    }
    return false;
}

// TODO
std::pair<std::vector<unsigned int>,std::vector<unsigned int>> is_activated_bond(const RDKit::ROMol &mol, const RDKit::Bond &bond)
{

    std::vector<unsigned int> tag_numbers = {75u,76u};

    //
    std::vector<unsigned int> bgn_activators = {};
    std::vector<unsigned int> end_activators = {};

    //
    std::vector<unsigned int> bgn_isotopes = {};
    std::vector<unsigned int> end_isotopes = {};

    // Get the atoms forming the bond
    const RDKit::Atom *bond_bgn = bond.getBeginAtom();
    const RDKit::Atom *bond_end = bond.getEndAtom();

    // Iterate over all neighboring atoms
    for(const auto &neighbor_atom: mol.atomNeighbors(bond_bgn))
    {
        //
        const unsigned int neighbor_atom_elem = neighbor_atom->getAtomicNum();
        const unsigned int neighbor_isotope = neighbor_atom->getIsotope();

        //
        if (count(tag_numbers.begin(), tag_numbers.end(), neighbor_atom_elem) > 0)
        {
            //std::cout << "Tag!" << std::endl;
            bgn_activators.emplace_back(neighbor_atom_elem);
            bgn_isotopes.emplace_back(neighbor_isotope);
        }   
    }

    // Iterate over all neighboring atoms
    for(const auto &neighbor_atom: mol.atomNeighbors(bond_end))
    {
        //
        const unsigned int neighbor_atom_elem = neighbor_atom->getAtomicNum();
        const unsigned int neighbor_isotope = neighbor_atom->getIsotope();

        //
        if (count(tag_numbers.begin(), tag_numbers.end(), neighbor_atom_elem) > 0)
        {
            //std::cout << "Tag!" << std::endl;
            //
            end_activators.emplace_back(neighbor_atom_elem);
            end_isotopes.emplace_back(neighbor_isotope);
        }   
    }

    // If either of the atoms has no tags, the bond is not activated
    if (bgn_activators.size() == 0 || end_activators.size() == 0)
    {
        // The bond is not activated
        return std::make_pair(std::vector<unsigned int>(),std::vector<unsigned int>());
    }
    
    // Sort the activator tag vectors
    std::sort(bgn_activators.begin(), bgn_activators.end());
    std::sort(end_activators.begin(), end_activators.end());

    // Sort the isotope vectors
    std::sort(bgn_isotopes.begin(), bgn_isotopes.end());
    std::sort(end_isotopes.begin(), end_isotopes.end());
    
    /*
    for (const unsigned int& num : bgn_isotopes)
    {
        std::cout << num << "\t";
    }
    std::cout << std::endl;

    for (const unsigned int& num : end_isotopes)
    {
        std::cout << num << "\t";
    }
    std::cout << std::endl;
    */

    // What does this do?
    bgn_activators.erase(std::unique(bgn_activators.begin(), bgn_activators.end()), bgn_activators.end());
    end_activators.erase(std::unique(end_activators.begin(), end_activators.end()), end_activators.end());
    
    bgn_isotopes.erase(std::unique(bgn_isotopes.begin(), bgn_isotopes.end()), bgn_isotopes.end());
    end_isotopes.erase(std::unique(end_isotopes.begin(), end_isotopes.end()), end_isotopes.end());
    
    // Template syntax : std::includes(vec1.begin(), vec1.end(), vec2.begin(), vec2.end())
    //if (std::equal(bgn_activators.begin(), bgn_activators.end(), end_activators.begin()) && \
    //    std::equal(bgn_isotopes.begin(), bgn_isotopes.end(), end_isotopes.begin()))
    
    std::vector<unsigned int> small_activators;
    std::vector<unsigned int> large_activators;
    std::vector<unsigned int> small_isotopes;
    std::vector<unsigned int> large_isotopes;
    if (bgn_activators.size() > end_activators.size())
    {
        //std::cout << "is larger" << std::endl;
        small_activators = end_activators;
        large_activators = bgn_activators;
    }
    else
    {
        //std::cout << "is smaller" << std::endl;
        small_activators = bgn_activators;
        large_activators = end_activators;
    }
    if(bgn_isotopes.size() > end_isotopes.size())
    {
        small_isotopes = end_isotopes;
        large_isotopes = bgn_isotopes;
    }
    else
    {
        small_isotopes = bgn_isotopes;
        large_isotopes = end_isotopes;
    }
    bool is_subset_activators = true;
    bool is_subset_isotopes = true;
    
    //std::cout << "s_a: " << small_activators.size() << std::endl;
    //std::cout << "l_a: " << large_activators.size() << std::endl;
    //std::cout << "s_i: " << small_isotopes.size() << std::endl;
    //std::cout << "l_i: " << large_isotopes.size() << std::endl;

    for (const unsigned int& act_num : small_activators)
    {
        if (count(large_activators.begin(), large_activators.end(), act_num) == 0)
        {
            is_subset_activators = false;
            //std::cout << "not subset activators" << std::endl;
            break;
        }
    }
    for (const unsigned int& iso_num : small_isotopes)
    {
        if (count(large_isotopes.begin(), large_isotopes.end(), iso_num) == 0)
        {
            is_subset_isotopes = false;
            //std::cout << "not subset isotopes" << std::endl;
            break;
        }
    }
    
    //
    if (is_subset_isotopes && is_subset_activators)
    {
        //std::cout << bond.getIdx() << " Is activated !!!" << std::endl;
        return std::make_pair(bgn_activators,end_activators);
    }    
    else
    {
        //std::cout << bond.getIdx() << " Not activated !!!" << std::endl;
        // The bond is not activated
        return std::make_pair(std::vector<unsigned int>(),std::vector<unsigned int>());
    }
}

// Check if the parsed bond is bound to activator tags
bool is_bond_bound_to_activators(const RDKit::ROMol &mol, const RDKit::Bond &bond, \
                                 const std::unordered_map<unsigned int, std::unordered_set<unsigned int>> &bond_adjacencies,\
                                 const std::unordered_map<unsigned int, std::unordered_set<unsigned int>> &atom_adjacencies)
{
    // Element numbers to exclude as regular heavy atoms
    std::unordered_set<unsigned int> exclude_elem = {75u,76u};

    // Check if this bond is activated!
    std::pair<std::vector<unsigned int>,std::vector<unsigned int>> bond_info = is_activated_bond(mol, bond);
    if (bond_info.first.size() == 0 && bond_info.second.size() == 0)
    {
        //std::cout << "bond is not activated :-(" << std::endl;
        return false;
    }
    else
    {
        ;//std::cout << "bond is activated :-)" << std::endl;
    }
    
    // Get the atoms forming the bond
    const RDKit::Atom *bond_bgn = bond.getBeginAtom();
    const RDKit::Atom *bond_end = bond.getEndAtom();

    // Get the indices of the atoms that form the bond
    const unsigned int bond_bgn_idx = bond_bgn->getIdx();
    const unsigned int bond_end_idx = bond_end->getIdx();

    // Iterate over all neighbor bond indices
    for (const unsigned int& neighbor_bond_idx: bond_adjacencies.at(bond.getIdx()))
    {
        //std::cout << "Neighbor bond index is:\t" << neighbor_bond_idx << std::endl;

        // Get the bond with this index
        const RDKit::Bond *neighbor_bond = mol.getBondWithIdx(neighbor_bond_idx);

        // Check if the bond is activated
        std::pair<std::vector<unsigned int>,std::vector<unsigned int>> neighbor_bond_info = is_activated_bond(mol, *neighbor_bond);
        if (neighbor_bond_info.first.size() == 0 && neighbor_bond_info.second.size() == 0)
        {
            //std::cout << "Neighbor was not activated" << std::endl;
            continue;
        }

        //std::cout << "I SURVIVED" << std::endl;

        // At this bond, our bond and a neighboring bond are activated
        const RDKit::Atom *neighbor_bond_bgn = neighbor_bond->getBeginAtom();
        const RDKit::Atom *neighbor_bond_end = neighbor_bond->getEndAtom();
        
        //
        const unsigned int neighbor_bond_bgn_idx = neighbor_bond_bgn->getIdx();
        const unsigned int neighbor_bond_end_idx = neighbor_bond_end->getIdx();
        
        // Store all atom indices
        std::vector<unsigned int> atom_indices = {bond_bgn_idx,bond_end_idx,neighbor_bond_bgn_idx,neighbor_bond_end_idx};

        unsigned int common_atom_idx;
        // Iterate over all atom indices in the bond and its neighbor
        for (const unsigned int& j : atom_indices)
        {
            // Find the index of the atom the neighboring bonds have in common
            if (count(atom_indices.begin(), atom_indices.end(),j) == 2)
            {
                common_atom_idx = j;
            }
        }

        // Iterate over all neighboring atoms of the common atom
        unsigned int sum = 0;
        for (const unsigned int& neighbor_atom_idx: atom_adjacencies.at(common_atom_idx))
        {
            //
            const RDKit::Atom *neighbor_atom = mol.getAtomWithIdx(neighbor_atom_idx);
            const unsigned int neighbor_atom_elem = neighbor_atom->getAtomicNum();

            //
            if (exclude_elem.find(neighbor_atom_elem) != exclude_elem.end())
            {
                ++sum;
            }
        }

        // At max two rings to be fused onto an atom
        if (sum < 2)
        {
            return false;
        }
    }

    return true;
}

// Check what type of activation tag is bound to parsed atom
std::vector<std::pair<unsigned int, unsigned int>> determine_atom_activation(const RDKit::ROMol &mol, const RDKit::Atom &atom,\
                                                                             std::unordered_set<std::string> &seen_activations)
{
    // Vector of activation pairs: activation tag atom idx and activation tag element
    std::vector<std::pair<unsigned int, unsigned int>> operators = {};

    // Iterate over all neighboring atoms of this atom
    for(const auto &neighbor_atom: mol.atomNeighbors(&atom))
    {
        // Get the atomic number of the neighbor
        const unsigned int neighbor_atom_elem = neighbor_atom->getAtomicNum();

        // Is the neighbor an activation tag
        if (activation_types.count(neighbor_atom_elem) > 0)
        {
            // Get the atomic index of the neighbor
            const unsigned int neighbor_idx = neighbor_atom->getIdx();

            // Get the isotope
            const unsigned int neighbor_isotope = neighbor_atom->getIsotope();

            // Convert information about activation tag to a string
            std::string activation = std::to_string(neighbor_isotope) + "_" + std::to_string(neighbor_atom_elem);

            // Check if the activation information string has been observed
            if (seen_activations.count(activation) == 0)
            {
                // Store a pair with the obtained info
                operators.emplace_back(std::make_pair(neighbor_idx, neighbor_atom_elem));

                // Store the activation information string
                seen_activations.emplace(activation);
            }
        }
    }
    return operators;
}

// Check what type of activation tag is bound to parsed bond
std::vector<std::pair<unsigned int, unsigned int>> determine_bond_activation(const RDKit::ROMol &mol, const RDKit::Bond &bond,\
                                                                             std::unordered_set<std::string> &seen_activations)
{
    // Vector of activation pairs: activation tag bond idx and activation tag element
    std::vector<std::pair<unsigned int, unsigned int>> operators = {};

    // Get the atoms forming the bond
    const RDKit::Atom *bond_bgn = bond.getBeginAtom();

    // Iterate over all neighboring atoms
    for(const auto &neighbor_atom: mol.atomNeighbors(bond_bgn))
    {
        // Get the atomic number of the neighbor
        const unsigned int neighbor_atom_elem = neighbor_atom->getAtomicNum();

        // Is the neighbor an activation tag
        if (activation_types.count(neighbor_atom_elem) > 0)
        {
            // Get the atomic index of the neighbor
            const unsigned int neighbor_idx = neighbor_atom->getIdx();

            // Get the isotope
            const unsigned int neighbor_isotope = neighbor_atom->getIsotope();

            // Convert information about activation tag to a string
            std::string activation = std::to_string(neighbor_isotope) + "_" + std::to_string(neighbor_atom_elem);

            // Check if the activation information string has been observed
            if (seen_activations.count(activation) == 0)
            {
                // Store a pair with the obtained info
                operators.emplace_back(std::make_pair(neighbor_idx, neighbor_atom_elem));

                // Store the activation information string
                seen_activations.emplace(activation);
            }
        }   
    }
    return operators;
}

// Standardize the given fingerprint based on atom symmetries and symmetry classes
std::string standardize_fingerprint(std::vector<std::tuple<unsigned int,unsigned int, unsigned int>>& fingerprint,\
                                    std::unordered_map<unsigned int,unsigned int>& atom_symmetries,\
                                    std::unordered_map<unsigned int,std::vector<unsigned int>>& symmetry_classes)
{
    // This is the output string
    std::string standardize_fingerprint;

    // Vector corresponding to standardized fingerprint after sort
    std::vector<std::tuple<int,int,int>> standardized_vector = {};

    // Iterate over all activation tuples
    for (const auto& information_tuple : fingerprint) 
    {     
        // Check if there is an atom with equivalent symmetry
        std::vector<unsigned int> eval_set = symmetry_classes[atom_symmetries[std::get<0>(information_tuple)]];

        // Get the minimal element
        int min_elem = *(--eval_set.rend());

        //
        standardized_vector.emplace_back(std::make_tuple(min_elem, std::get<1>(information_tuple), std::get<2>(information_tuple)));
    }

    // Sort the standardized fingerprint
    std::sort(standardized_vector.begin(), standardized_vector.end()); //, stepwise_sorter);

    // Format the fingerprint string
    for (const auto& information_tuple : standardized_vector)
    {
        //
        std::string token = std::to_string(std::get<0>(information_tuple)) + "@" + \
                            std::to_string(std::get<1>(information_tuple)) + "#" + \
                            std::to_string(std::get<2>(information_tuple)) + "-";
        standardize_fingerprint.append(token);
    }

    // Remove the last junction character
    standardize_fingerprint.erase(standardize_fingerprint.size() - 1);

    return standardize_fingerprint;
}

// Count number of heavy atoms in molecule
// but exclude activator tag atoms
unsigned int count_scaffold_heavy_atoms(const RDKit::ROMol& mol)
{

    // Element numbers to exclude as regular heavy atoms
    std::unordered_set<int> exclude_elem = {1,71,72,73,74,75,76};

    // Count number of heavy atoms
    unsigned int num_heavy = 0;

    // Iterate over all atoms in the molecule
    const unsigned int num_atoms = mol.getNumAtoms();
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx) 
    {
        // Retrieve the molecule's atom
        const RDKit::Atom *atom = mol.getAtomWithIdx(atom_idx);

        // Find out what element the atom is
        int elem_num = static_cast<int>(atom->getAtomicNum());

        // Check if the atom is a regular heavy atom
        if (exclude_elem.find(elem_num) == exclude_elem.end())
        {
            // Increase the counter
            ++num_heavy;
        }
        
    }
    return num_heavy;
}

// Driver function
int main(int argc, char *argv[])
{
    
    // File names
    std::string in_file_path;
    std::string out_file_path;
    
    // Total number of heavy atoms in molecule
    unsigned int num_heavy;

    // Boolean to control verbose output
    bool verbose = false;

    // String to inform user
    std::string usage_string = " -i <input_file_path>\n\
                            -o <output_file_path>\n\
                            -n <num_heavy>\n\
                            -v <verbose>";

    // Argument parser
    int opt;
    while ((opt = getopt(argc, argv, "i:o:n:v:h")) != -1)
    {
        switch (opt)
        {
            case 'i':
                in_file_path = optarg;
                break;
            case 'o':
                out_file_path = optarg;
                break;
            case 'n':
                num_heavy = static_cast<unsigned int>(std::stoi(optarg));
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                std::cout << "Usage: " << argv[0] << usage_string << std::endl;
                return 0;
            default:
                std::cerr << "Usage: " << argv[0] << usage_string << std::endl;
                return 1;
        }
    }

    // Verify path
    if (in_file_path.empty())
    {
        std::cout << "Usage: " << argv[0] << usage_string << std::endl;
        return 1;
    }

    // Verify path
    if (out_file_path.empty())
    {
        std::cout << "Usage: " << argv[0] << usage_string << std::endl;
        return 1;
    }
    
    // Verify you can open input file
    std::ifstream input_stream(in_file_path);
    if (!input_stream)
    {
        std::cout << "Error opening input file: " << in_file_path << std::endl;
        return 1;
    }

    // Verify you can open output file
    std::ofstream output_stream(out_file_path);
    if (!output_stream)
    {
        std::cout << "Error opening output file: " << out_file_path << std::endl;
        return 1;
    }

    // String for reading lines
    std::string line;

    // The minimal number of atoms a substituent needs to have for
    // the intended activation
    // TODO read this in and make modular
    std::unordered_map<std::string,int> minimum_heavy_atoms;
    minimum_heavy_atoms["single"] = 1;
    minimum_heavy_atoms["double"] = 1;
    minimum_heavy_atoms["triple"] = 1;
    minimum_heavy_atoms["single_fusion"] = 3 - 2; // two atoms consumed by fusing bonds
    minimum_heavy_atoms["double_fusion"] = 5 - 2; // two atoms consumed by fusing bonds
    minimum_heavy_atoms["spiro"] = 3 - 1; // one atom consumed by merging rings in a point

    // Read the input file line by line
    while (std::getline(input_stream, line))
    {
        // Turn the line (SMILES) into a molecule
        std::shared_ptr<RDKit::ROMol> scaffold(RDKit::SmilesToMol(line));

        // This is the point where you set the hydrogens to explicit mode
        std::shared_ptr<RDKit::ROMol> protonated_scaffold(RDKit::MolOps::addHs(*scaffold));

        // Count the number of heavy atoms in the starting scaffold, but do not consider 
        // activator tags
        unsigned int num_scaffold_heavy = count_scaffold_heavy_atoms(*protonated_scaffold);

        // How many new heavy atoms will we introduce in total
        int num_introducing_heavy = num_heavy - num_scaffold_heavy;

        // Some verbose if requested
        if (verbose)
        {
            std::cout << "Introducing:\t" << num_introducing_heavy << " Heavy:\t" << num_heavy << " scaffold:\t" << num_scaffold_heavy << std::endl;
        }

        // Perform hidden symmetry analysis by assigning labels to atoms
        // cleanIt=True,force=True,flagPossibleStereoCenters=True
        RDKit::MolOps::assignStereochemistry(*protonated_scaffold,true,true,true);

        // Perform symmetry analysis
        std::unordered_map<unsigned int,unsigned int> atom_symmetries = get_atom_symmetries(*protonated_scaffold);

        // Retrieve the symmetry factors of the atoms
        std::unordered_map<unsigned int,unsigned int> symmetry_factors = get_symmetry_factors(atom_symmetries);

        // Retrieve the symmetry classes and all atoms belonging to these classes
        std::unordered_map<unsigned int,std::vector<unsigned int>> symmetry_classes = get_symmetry_classes(atom_symmetries);

        // Retrieve atom and bond adjacency maps
        std::unordered_map<unsigned int, std::unordered_set<unsigned int>> atom_adjacencies = get_atom_adjacencies(*protonated_scaffold);
        std::unordered_map<unsigned int, std::unordered_set<unsigned int>> bond_adjacencies = get_bond_adjacencies(*protonated_scaffold);

        // A set to store activation patterns
        std::unordered_set<std::string> seen_activations;

        // Find all atoms that are connected to an activator, determine their type of activation
        std::map<int,std::vector<std::pair<unsigned int, unsigned int>>> activation_atom_library = {};

        // Get the number of atoms in this molecule
        const unsigned int num_atoms = protonated_scaffold->getNumAtoms();

        // Iterate over all atoms in the molecule
        for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx) 
        {
            // Retrieve the atom with this index
            const RDKit::Atom *activate_atom = protonated_scaffold->getAtomWithIdx(atom_idx);

            // Check if this atom carries an activation tag
            if (is_atom_bound_to_activator(*protonated_scaffold, *activate_atom))
            {
                // Map the type of activation with the atom's index
                activation_atom_library.emplace(atom_idx, determine_atom_activation(*protonated_scaffold, *activate_atom, seen_activations));
            }           
        }

        // Find all bonds that are connected to activators, determine the type of activation
        // When to activated bonds are adjacent and have one bond in between, that bond is not
        // treated as activated, as it would leave two activators tags on the end not used
        std::map<int,std::vector<std::pair<unsigned int, unsigned int>>> activation_bond_library = {};

        // Get the number of bonds in this molecule
        const unsigned int num_bonds = protonated_scaffold->getNumBonds();

        // Iterate over all bonds in the molecule
        for (unsigned int bond_idx = 0 ; bond_idx < num_bonds ; ++bond_idx) 
        {
            // Retrieve the bond with this index
            const RDKit::Bond *activate_bond = protonated_scaffold->getBondWithIdx(bond_idx);

            // Check if this bond carries activation tags
            if (is_bond_bound_to_activators(*protonated_scaffold, *activate_bond, bond_adjacencies, atom_adjacencies))
            {
                // Map the type of activation with the bond's index
                activation_bond_library.emplace(bond_idx, determine_bond_activation(*protonated_scaffold, *activate_bond, seen_activations));
            }    
        }

        // Some verbose if requested
        if (verbose)
        {
            std::cout << "Num activated atoms:\t" << activation_atom_library.size() << std::endl;
            std::cout << "Num activated bonds:\t" << activation_bond_library.size() << std::endl;
        }
        
        // Num_connections and connection_indices should be the same length
        // Get all indices of atoms that have an activator tag
        std::vector<int> connection_atom_indices = {};
        std::vector<int> connection_bond_indices = {};

        // Num_connections and connection_indices should be the same length
        // Get all indices of atoms that have an activator tag
        std::vector<int> connection_atom_neighbors = {};
        std::vector<int> connection_bond_neighbors = {};

        // Num_connections and connection_indices should be the same length
        // Get all indices of atoms that have an activator tag
        std::vector<int> connection_atom_operations = {};
        std::vector<int> connection_bond_operations = {};

        // Iterate over all atom activations
        for (const auto& activation_pair : activation_atom_library)
        {       
            // Iterate over all activations for this atom
            int num_atom_activations = static_cast<int>(activation_atom_library[activation_pair.first].size());
            for (int i = 0; i < num_atom_activations; i++) 
            {    
                // Add the index of this atom
                connection_atom_indices.emplace_back(activation_pair.first);
                // Activation tag atom index
                connection_atom_neighbors.emplace_back(activation_pair.second[i].first);
                // Activation tag atom element (says what type of activation)
                connection_atom_operations.emplace_back(activation_pair.second[i].second);
            }
        }

        // Iterate over all bond activations
        for (const auto& activation_pair : activation_bond_library)
        {       
            // Iterate over all activations for this bond
            int num_bond_activations = static_cast<int>(activation_bond_library[activation_pair.first].size());
            for (int i = 0; i < num_bond_activations; i++)
            {    
                // Add the index of this bond
                connection_bond_indices.emplace_back(activation_pair.first);
                // Activation tag atom indices
                connection_bond_neighbors.emplace_back(activation_pair.second[i].first);
                // Activation tag atom element (says what type of activation)
                connection_bond_operations.emplace_back(activation_pair.second[i].second);
            }
        }

        // Debug purposes
        //for (const auto& activation_pair : activation_atom_library)
        //{
        //    int num_atom_activations = static_cast<int>(activation_atom_library[activation_pair.first].size());
        //    for (int i = 0; i < num_atom_activations; i++) 
        //    {
        //        ;//std::cout << activation_pair.second[i].second << std::endl;
        //    } 
        //}

        // A vector to store all atom and bond (bgn) indices of sites that are activated
        std::vector<int> connection_indices = {};
        // Iterate over the activated atom sites
        for (const int& index : connection_atom_indices)
        {
            connection_indices.emplace_back(index);
        }
        // Iterate over the activated bond (bgn) sites
        for (const int& index : connection_bond_indices)
        {
            connection_indices.emplace_back(index);
        }

        // A vector to store all atom and bond (bgn) activation tag indices
        std::vector<int> connection_neighbors = {};
        // Iterate over the atom activation tag indices
        for (const int& neighbor : connection_atom_neighbors)
        {
            connection_neighbors.emplace_back(neighbor);
        }
        // Iterate over the bond (bgn) activation tag indices
        for (const int& neighbor : connection_bond_neighbors)
        {
            connection_neighbors.emplace_back(neighbor);
        }

        // A vector to store all atom and bond (bgn) activation tag elements
        std::vector<int> connection_operations = {};
        // Iterate over the atom activation tag elements
        for (const int& operation : connection_atom_operations)
        {
            connection_operations.emplace_back(operation);
        }
        // Iterate over the bond (bgn) activation tag elements
        for (const int& operation : connection_bond_operations)
        {
            connection_operations.emplace_back(operation);
        }

        // Numbers of atom activations and bond activations
        int num_atom_activations = static_cast<int>(connection_atom_indices.size());
        int num_bond_activations = static_cast<int>(connection_bond_indices.size());

        int num_atom_op_activations = static_cast<int>(connection_atom_operations.size());
        int num_bond_op_activations = static_cast<int>(connection_bond_operations.size());

        // Total number of activations
        int num_connections = num_atom_activations + num_bond_activations;
        int num_operations = num_atom_op_activations + num_bond_op_activations;

        // Some verbose if requested
        if (verbose)
        {
            std::cout << "Num indices atoms:\t" << num_atom_activations << std::endl;
            std::cout << "Num indices bonds:\t" << num_bond_activations << std::endl;
            std::cout << "Number of connection sites:\t" << num_connections << std::endl;
            std::cout << "Num operations atoms:\t" << num_atom_op_activations << std::endl;
            std::cout << "Num operations bonds:\t" << num_bond_op_activations << std::endl;
            std::cout << "Number of connection sites:\t" << num_operations << std::endl;
        }

        // Avoid doing the same partitions again, if the activated atoms are symmetrically identical
        std::unordered_set<std::string> seen_activation_fingerprints = {};

        // Iterate over all possible combinations to attach n new heavy atoms to k connection points
        // For instance, distributing 5 atoms over 3 connection points can be done in the following way:
        // [1,2,2]
        // [2,1,2]
        // [2,2,1]
        // [1,1,3]
        // [1,3,1]
        // [3,1,1]

        // Each connection needs to get at least 1 new heavy atom, otherwise no substitution happened
        std::vector<std::vector<int>> partitions = generate_partitions(num_introducing_heavy, num_connections);

        // Iterate over the atom partitions
        for (const std::vector<int>& atom_distribution_over_connections : partitions)
        {    
            // Construct activation fingerprint from number of added atoms and activated position indices
            std::vector<std::tuple<unsigned int,unsigned int,unsigned int>> fingerprint = {};
            int a = 0;
            for (const int& num_atoms_to_add : atom_distribution_over_connections)
            {
                // Tuple contains site index + number of atoms to add + what type of activation
                fingerprint.emplace_back(std::make_tuple(connection_indices[a], num_atoms_to_add, connection_operations[a]));
                ++a;
            }

            // See if this fingerprint has been evaluated before in terms of symmetrical atoms
            //std::string standard_fingerprint = standardize_fingerprint(fingerprint, atom_symmetries, symmetry_classes);

            // If this operation has not been observed before
            //if (seen_activation_fingerprints.count(standard_fingerprint) == 0)
            //{
                // Store the activation fingerprint
                //seen_activation_fingerprints.emplace(standard_fingerprint);
            //}
            // This configuration has done before, skip loop and save time
            //else
            //{
                //continue; // implement later
            //}
            
            /*
            // For debugging purposes now
            int count = -1;
            int num_p = atom_distribution_over_connections.size();
            std::cout << "Partition: [";
            for (const int& p : atom_distribution_over_connections)
            {
                ++count;
                std::cout << p;
                if (count < num_p)
                {
                    std::cout << ",";
                }                
            }
            std::cout << "]" << std::endl;
            */

            // Boolean to control if we can actually modify the molecule
            bool can_forge = true;

            // String variable that holds the type of activation info
            std::string activation_type;

            // Create a work molecule that we can modify by attaching substituents
            std::shared_ptr<RDKit::RWMol> scaffold2substitute(new RDKit::RWMol(*protonated_scaffold));
            
            //
            std::unordered_set<std::string> seen_tags;

            // Iterate over all connections in this partition
            // c is the index of the activated site in the operation array, 
            // position 0, 1, 2, but not the index of atoms or bonds in the molecule!
            int c = 0; 

            // A string containing the instructions to make superstructures
            std::string job_string = line; // Starts of with the activated scaffold SMILES

            // Iterate over all operations that add atom to the scaffold
            for (const int& num_atoms_to_add : atom_distribution_over_connections)
            {    
                // Verify that this connection point gets the minimum number
                // of heavy atoms needed for the activation. 
                int connection_index;
                int operation_index;
                int neighbor_index;

                // TODO
                std::vector<std::pair<unsigned int, unsigned int>> activation_pattern;

                // Check if the activation site relates to atoms or bonds
                // We're coupling to a single atom
                if (c < num_atom_activations)
                {
                    activation_type = "ATOM";
                    // For the current activation index, find the corresponding atom index
                    connection_index = connection_indices[c];
                    neighbor_index = connection_neighbors[c];
                    operation_index = connection_operations[c];
                }
                // We're fusing a ring
                else
                {
                    activation_type = "BOND";
                    // For the current activation index, find the corresponding bond index
                    connection_index = connection_indices[c];
                    neighbor_index = connection_neighbors[c];
                    operation_index = connection_operations[c];
                }

                /*
                // Iterate over all activations in this pattern
                for (const std::pair<unsigned int,unsigned int>& activation : activation_pattern)
                {
                    // The minimum number of heavy atoms an incoming substituent must have
                    int num_req_heavy_atoms = minimum_heavy_atoms[activation_types[activation.second]];

                    //std::cout << "num:\t" << num_req_heavy_atoms << "," << num_atoms_to_add << std::endl;

                    // If the activation requires more atoms than were intended, break
                    if (num_req_heavy_atoms > num_atoms_to_add)
                    {
                        //std::cout << "not enough atoms" << std::endl;
                        // Make sure we break the parent loop too
                        can_forge = false;
                        break;
                    }
                }
                */

                // The minimum number of heavy atoms an incoming substituent must have
                int num_req_heavy_atoms = minimum_heavy_atoms[activation_types[operation_index]];

                // If the activation requires more atoms than we intend adding
                if (num_req_heavy_atoms > num_atoms_to_add)
                {
                    // Skip activation pattern entirely
                    can_forge = false;
                    break;
                }

                // At least one atom in this activation partition was not valid
                // break this partition loop
                //if (!can_forge)
                //{
                //    break;
                //}
                
                //std::cout << "act size is:\t" << activation_pattern.size() << std::endl;

                // Don't think i need to loop over a pattern here, just assign the atom
                
                /*
                int k = 0;
                // Iterate over all activations in this pattern
                for (const std::pair<unsigned int,unsigned int>& activation : activation_pattern)
                {
                    // Find the atom
                    const RDKit::Atom *activation_tag = scaffold2substitute->getAtomWithIdx(activation.first);

                    // Retrieve the isotope 
                    const unsigned int isotope_number = activation_tag->getIsotope();

                    // Format the job string
                    std::string sub_string = " ***" + std::to_string(k) + " [" + std::to_string(isotope_number) + \
                    activation_tag->getSymbol() + "] " + std::to_string(num_atoms_to_add) + \
                    "." + activation_types[activation.second] + ".ism";

                    std::cout << sub_string << std::endl;

                    if (seen_tags.count(sub_string) == 0)
                    {
                        job_string.append(sub_string);
                        seen_tags.emplace(sub_string);
                    }

                    // Increment the counter
                    ++k;
                }
                */
               
                // Get the atom with this index
                const RDKit::Atom *activation_tag = scaffold2substitute->getAtomWithIdx(neighbor_index);

                // Retrieve the isotope of the activation tag
                const unsigned int isotope_number = activation_tag->getIsotope();

                // Format the job string
                std::string sub_string = " [" + std::to_string(isotope_number) + \
                activation_tag->getSymbol() + "] " + std::to_string(num_atoms_to_add) + \
                "." + activation_types[operation_index] + ".ism";

                // If this operation is not part of the job string yet
                if (seen_tags.count(sub_string) == 0)
                {
                    // Make it part of the job string
                    job_string.append(sub_string);
                    // Store the operation in the seen tags
                    seen_tags.emplace(sub_string);
                }

                // Increase the index counter
                ++c;
            }

            // If it is a valid job string
            if (can_forge)
            {
                // Write the job stirng to the output file
                output_stream << job_string << std::endl;
            }
        }
    }

    // Close the file streams
    input_stream.close();
    output_stream.close();

    // Signal success
    return 0;
}