/*
    Author: Andreas Luttens
    Contact: andreas.luttens@gmail.com
    Date: June 22, 2023
    Description: Introduce double and triple bonds in given saturated hydrocarbons
*/

// Basic utilities
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>
#include <utility>
#include <functional> 
#include <unistd.h> 

// RDKit functions 
#include <GraphMol/GraphMol.h>
#include <GraphMol/Canon.h> // Check for symmetry later
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

// A map between bondorders and corresponding bondtypes
std::unordered_map<unsigned int, RDKit::Bond::BondType> bond_types = {{1u, RDKit::Bond::SINGLE},
                                                                      {2u, RDKit::Bond::DOUBLE},
                                                                      {3u, RDKit::Bond::TRIPLE}};

//
void generateCombinationsHelper(const std::unordered_set<unsigned int>& nums, int k, int start, 
                                std::vector<int>& currentCombination, 
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
std::vector<std::vector<int>> generateCombinations(const std::unordered_set<unsigned int>& nums, int k)
{
    std::vector<std::vector<int>> combinations;
    std::vector<int> currentCombination;
    generateCombinationsHelper(nums, k, INT_MIN, currentCombination, combinations);
    return combinations;
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

// TODO This function says nothing about which bridgehead
bool can_unsat_bridgehead(std::unordered_set<unsigned int>& bridgeheads, 
                          std::vector<std::unordered_set<unsigned int>>& cycles, 
                          const unsigned int& num_bredt_violations)
{
    //
    unsigned int limit;
    if (num_bredt_violations == 0)
    {
        limit = 8;
    }
    else
    {
        limit = 10;
    }
    
    //
    for (const unsigned int& bridgehead_idx : bridgeheads) 
    {    
        // 
        unsigned int max_ring_size = 0;

        for (const std::unordered_set<unsigned int>& cycle : cycles) 
        {

            // is FIND faster, just return first, not all hits
            if (cycle.find(bridgehead_idx) != cycle.end())
            {
                // Find the ringsize
                unsigned int ring_size = static_cast<unsigned int>(cycle.size());
                
                if (ring_size > max_ring_size)
                {
                    max_ring_size = ring_size;
                }   
            }
        }

        //
        if (max_ring_size < limit)
        {
            return false;
        } 
    }

    //
    return true;
}


// Function does a quick check, all atoms in bond must be carbon
// and the bondorder needs to matched the given bondorder
bool is_viable_bond(const RDKit::Bond& bond, const unsigned int& new_bond_order){

    // Check the bondorder, needs to be a sigma bond
    // We don't turn a double bond into triple bonds
    // we go from single to triple in one step
    // we also exclude aromatic bonds this way
    const unsigned int bondtype = bond.getBondType();
    if (bondtype != 1)
    {
        return false;
    }
    
    // Get the atoms forming the bond
    const RDKit::Atom *bond_bgn = bond.getBeginAtom();
    const RDKit::Atom *bond_end = bond.getEndAtom();
    
    // Get the elements of these atoms
    const unsigned int bond_bgn_elem = bond_bgn->getAtomicNum();
    const unsigned int bond_end_elem = bond_bgn->getAtomicNum();

    // If either atom is not a carbon, false
    if (bond_bgn_elem != 6 || bond_end_elem != 6)
    {
        return false;
    }

    // If we make a triple bond, we need two hydrogens per atom
    // if we make a double bond, we need one hydrogen per atom
    if (bond_bgn->getTotalNumHs() < (new_bond_order - 1) || bond_end->getTotalNumHs() < (new_bond_order - 1))
    {
        return false;
    }

    return true;
}

// Function that checks if unsaturation can be introduced in four-membered ring, bond has to be exo, not endo
bool can_unsat_four_ring(const RDKit::ROMol& mol, const int& begin_idx, const int& end_idx, std::unordered_set<unsigned int>& non_sp3_indices)
{
    // Both atoms cannot be in the same R4, would be endo unsaturation
    if (mol.getRingInfo()->areAtomsInSameRingOfSize(begin_idx, end_idx, 4))
    {
        return false;
    }

    // Iterate over rings in the molecule
    for(const auto &ring: mol.getRingInfo()->atomRings())
    {    
        // If it is not a four-mem ring, skip
        if (ring.size() != 4)
        {
            continue;
        }
        
        // Iterate over the atoms in the current four-mem ring
        // check if one has sp2 hybridization
        for(const int &atom_idx: ring)
        {   
            // If the current atom is in the bond (begin_idx, end_idx)
            // skip
            if (atom_idx == begin_idx || atom_idx == end_idx)
            {
                continue;
            }

            // If this atom and begin_idx are in the same four-mem ring
            if (mol.getRingInfo()->areAtomsInSameRingOfSize(begin_idx, atom_idx, 4))
            {                
                // Get the hybridization state of this atom
                if (non_sp3_indices.count(atom_idx) > 0)

                // It must be sp3 hybridized
                {
                    return false;
                }    
            }

            // If this atom and end_idx are in the same four-mem ring
            if (mol.getRingInfo()->areAtomsInSameRingOfSize(end_idx, atom_idx, 4))
            {
                // Get the hybridization state
                if (non_sp3_indices.count(atom_idx) > 0)
                {
                    return false;
                }    
            }
        }
    }

    return true;
}

//TODO
std::pair<std::unordered_set<unsigned int>,std::vector<std::unordered_set<unsigned int>>> find_bridgehead_atoms(RDKit::ROMol& mol)
{
    // Set containing all bridgehead atom indices
    std::unordered_set<unsigned int> bridgehead_indices = {};

    // SSSR indices
    const std::vector<RDKit::RingInfo::INT_VECT> sssr = mol.getRingInfo()->atomRings();

    // pattern to match sp3 hybridized ring atoms with 3 or 4 neighbors
    const RDKit::RWMol *pattern = RDKit::SmartsToMol("[x3,x4:1]"); // can move this outside of function

    // Init a vector for matches
    std::vector<RDKit::MatchVectType> matches;    

    // Perform substructure match
    RDKit::SubstructMatch(mol, *pattern, matches);
    
    //
    std::vector<unsigned int> standard_matches = {};

    // Get the number of rings
    const unsigned int num_rings = static_cast<unsigned int>(sssr.size());

    //
    std::set<std::vector<unsigned int>> all_intersections;
    
    // Iterate over all cycles
    for (unsigned int a_ring_idx = 0 ; a_ring_idx < num_rings ; ++a_ring_idx)
    {    
        for (unsigned int b_ring_idx = 0 ; b_ring_idx < num_rings ; ++b_ring_idx)
        {    
            if (a_ring_idx > b_ring_idx)
            {     
                std::vector<int> intersection;
                
                RDKit::Intersect(sssr[a_ring_idx], sssr[b_ring_idx], intersection);

                std::vector<unsigned int> add_intersection = {};

                // Add all elements
                for (const int& inter_idx : intersection){
                    add_intersection.emplace_back(static_cast<unsigned int>(inter_idx));
                }

                all_intersections.emplace(add_intersection);
            }
        }
    }

    // Iterate over all ring intersections
    for (const auto& intersection : all_intersections)
    {    
        
        // Condition for bridgehead
        if (intersection.size() > 2)
        {
            
            for (const auto& atom_idx : intersection)
            {    
                const RDKit::Atom *atom = mol.getAtomWithIdx(atom_idx);

                // Correct type
                std::unordered_set<unsigned int> addition = {};

                // Iterate over the neighbor atoms of the current atom
                for(const auto &neighbor_atom: mol.atomNeighbors(atom))
                {    
                    // Correct type
                    const unsigned int neighbor_atom_idx = neighbor_atom->getIdx();

                    if (count(intersection.begin(), intersection.end(), neighbor_atom_idx) == 0 && 
                        count(standard_matches.begin(), standard_matches.end(), atom_idx) > 0)
                    {
                        addition.emplace(atom_idx);
                    }
                }

                // Update the bridgehead atom indices
                bridgehead_indices.insert(addition.begin(), addition.end());
            }
        }
    }

    // Change SSSR into cycles here
    std::vector<std::unordered_set<unsigned int>> cycles = {};

    for (unsigned int ring_idx = 0 ; ring_idx < num_rings ; ++ring_idx)
    {
        std::unordered_set<unsigned int> atom_set = {};
        for(const auto &atom_idx: sssr[ring_idx])
        {
            atom_set.emplace(static_cast<unsigned int>(atom_idx));
        }
        cycles.emplace_back(atom_set);
    }

    std::pair<std::unordered_set<unsigned int>,std::vector<std::unordered_set<unsigned int>>> information = std::make_pair(bridgehead_indices, cycles);

    return information;
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

// A function to identify bonds that are equivalent in terms of their composing atoms' CIP rank
std::unordered_map<unsigned int,std::vector<unsigned int>> get_equivalent_bonds(const RDKit::ROMol& mol, 
                                                                                std::unordered_map<unsigned int,unsigned int>& atom_symmetries){

    // 
    std::map<std::pair<unsigned int,unsigned int>,std::vector<unsigned int>> bond_symmetries = {};

    // Iterate over all bonds in the molecule
    unsigned int num_bonds = mol.getNumBonds();
    for (unsigned int bond_idx = 0 ; bond_idx < num_bonds ; ++bond_idx) {

        // Get the bond with the current bond index
        const RDKit::Bond *bond = mol.getBondWithIdx(bond_idx);

        // Get the atom indices of the atoms forming this bond
        const unsigned int bond_bgn_idx = bond->getBeginAtom()->getIdx();
        const unsigned int bond_end_idx = bond->getEndAtom()->getIdx();

        // Get the atom symmetries
        const unsigned int bond_bgn_sym = atom_symmetries[bond_bgn_idx];
        const unsigned int bond_end_sym = atom_symmetries[bond_end_idx];

        //TODO can I make this a const
        // Make a pair of the atom symmetries, which will become the bond symmetry
        std::pair<unsigned int,unsigned int> bond_symmetry = std::make_pair(bond_bgn_sym,bond_end_sym);

        // Make the reverse pair as well
        std::pair<unsigned int,unsigned int> rev_bond_symmetry = std::make_pair(bond_end_sym,bond_bgn_sym);

        // If this bond symmetry has not been seen before
        if (bond_symmetries.count(bond_symmetry) == 0)
        {
            // Add the symmetry, and an empty set to store pairs later
            bond_symmetries.emplace(bond_symmetry, std::vector<unsigned int>{});
            bond_symmetries.emplace(rev_bond_symmetry, std::vector<unsigned int>{});

        }

        // Add the bond symmetry and corresponding bond index
        bond_symmetries[bond_symmetry].emplace_back(bond_idx);

        // We can add the reverse symmetry as well
        bond_symmetries[rev_bond_symmetry].emplace_back(bond_idx);
                
    }
    
    // Equivalent bonds should be sorted!, so std::vector
    std::unordered_map<unsigned int, std::vector<unsigned int>> equivalent_bonds;
    // Iterate over all bond symmetries
    for (const auto& sym_pair : bond_symmetries)
    {
        //
        std::vector<unsigned int> sym_bonds = sym_pair.second;
        std::sort(sym_bonds.begin(),sym_bonds.end());
        // Iterate over all bonds with this bond symmetry
        for (const auto& sym_bond: sym_bonds){

            equivalent_bonds.emplace(sym_bond, sym_bonds);
        }
    }

    //
    return equivalent_bonds;
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

// Standardize the given operation based on equivalent bonds
std::vector<int> standardize_operation(const std::vector<int>& operation, std::unordered_map<unsigned int,std::vector<unsigned int>>& equivalent_bonds)
{
    std::vector<int> standard_operation;

    // Iterate over all bond indices that will be unsaturated
    for (const int& bond_idx : operation) 
    { 
        // Get the minimal equivalent bond index, since the vector is sorted, it is the first element
        // if there were no equivalent bonds, it will return itself as sole (and first) element
        int min_equivalent_idx = equivalent_bonds[bond_idx].front();

        //
        standard_operation.emplace_back(min_equivalent_idx);
    }

    // sort the standardized fingerprint
    sort(standard_operation.begin(), standard_operation.end()); //, stepwise_sorter);

    return standard_operation;
}

//
std::unordered_set<unsigned int> get_nonsp3_indices(RDKit::ROMol& mol)
{
    //
    std::unordered_set<unsigned int> indices = {};

    // Iterate over all atoms in the molecule
    unsigned int num_atoms = mol.getNumAtoms();
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx) {

        // Get the atom with the current atom index
        const RDKit::Atom *atom = mol.getAtomWithIdx(atom_idx);

        if(atom->getHybridization() != RDKit::Atom::HybridizationType::SP3)
        {
            indices.emplace(atom_idx);
        }
    }
    return indices;
}

// Function to generate all unsaturations for a given input molecule and bondorder
std::unordered_set<std::string> generate_unsaturations(RDKit::ROMol& mol, unsigned int new_bond_order)
{    
    // A set of unique SMILES string for this hydrocarbon
    std::unordered_set<std::string> unique_smiles;

    // Add the original smiles
    unique_smiles.emplace(RDKit::MolToSmiles(mol));

    // Perform ring analysis
    RDKit::MolOps::symmetrizeSSSR(mol);

    // Store ring info, can I do a const here? This molecule won't change
    RDKit::RingInfo* ring_info = mol.getRingInfo();

    // Perform hidden symmetry analysis by assigning labels to atoms
    // cleanIt=True,force=True,flagPossibleStereoCenters=True
    RDKit::MolOps::assignStereochemistry(mol,true,true,true);

    // Perform atom symmetry analysis
    std::unordered_map<unsigned int,unsigned int> atom_symmetries = get_atom_symmetries(mol);

    // Get all equivalent bonds in terms of atom symmetry
    std::unordered_map<unsigned int,std::vector<unsigned int>> equivalent_bonds = get_equivalent_bonds(mol, atom_symmetries);

    // Initialize a set to store all viable C-C bonds
    std::unordered_set<unsigned int> viable_bond_indices = {};

    // Retrieve each bond's begin and end atom indices
    std::unordered_map<unsigned int,std::pair<unsigned int,unsigned int>> bond_information = get_bond_information(mol);

    // Retrieve each bond's adjacent baonds
    std::unordered_map<unsigned int,std::unordered_set<unsigned int>> bond_adjacencies = get_bond_adjacencies(mol);

    // Find the bridgehead atoms and list all rings in the molecule
    std::pair<std::unordered_set<unsigned int>,std::vector<std::unordered_set<unsigned int>>> information = find_bridgehead_atoms(mol); // turn this into a void!
    std::unordered_set<unsigned int> bridgeheads = information.first;
    std::vector<std::unordered_set<unsigned int>> cycles = information.second;

    // Find the carbon atoms that are not sp3 hybridized
    std::unordered_set<unsigned int> non_sp3_indices = get_nonsp3_indices(mol);

    // Iterate over all bonds in the molecule
    unsigned int num_bonds = mol.getNumBonds();
    for (unsigned int bond_idx = 0 ; bond_idx < num_bonds ; ++bond_idx)
    {
        // Get the bond with the current bond index
        const RDKit::Bond *bond = mol.getBondWithIdx(bond_idx);

        // Check whether this bond is viable for unsaturation
        if (!is_viable_bond(*bond, new_bond_order))
        {
            // Skip this bond if not viable
            continue;
        }
        
        // Get the atoms forming the bond
        const RDKit::Atom *bond_bgn = bond->getBeginAtom();
        const RDKit::Atom *bond_end = bond->getEndAtom();

        // Get the atom indices of these atoms
        const unsigned int bond_bgn_idx = bond_bgn->getIdx();
        const unsigned int bond_end_idx = bond_end->getIdx();

        // Bool to keep track that no adjacent double bonds are introduced
        bool is_allene_bond = false;

        // Iterate over all bonds adjacent to the current bond
        for (const unsigned int& neighbor_bond_idx : bond_adjacencies[bond_idx])
        {    
            // Get the neighbor bond
            const RDKit::Bond *neighbor_bond = mol.getBondWithIdx(neighbor_bond_idx);

            // If the neighboring bonds is not a sigma bond
            if (neighbor_bond->getBondType() != 1)
            {
                // Allene formation
                is_allene_bond = true;
                break;
            }  
        }

        // If a single pair of bonds is adjacent
        if (is_allene_bond)
        {
            continue;
        }
        
        // If the bond is in a ring
        if (ring_info->numBondRings(bond_idx))
        {
            
            // If the new bondorder is three
            if (new_bond_order == 3)
            {
                // The smallest ringsize must be 9 or higher
                if (ring_info->minBondRingSize(bond_idx) < 9)
                {
                    continue;
                }      
            }
            
            // Bridgeheads can only carry a double bond in some cases
            // a triple bond is never possible

            // Check if either atom in the bond is a bridgehead
            if (bridgeheads.find(bond_bgn_idx) != bridgeheads.end() || bridgeheads.find(bond_end_idx) != bridgeheads.end())
            {
                // Check if bond violates Bredt's rule 
                // Bridgeheads and cycles are constant for this molecule, so
                // I should be able to move this function outside of the loop?
                if (!can_unsat_bridgehead(bridgeheads, cycles, 0))
                {
                    continue;
                }
            }
        }

        // Check if either atom in the bond is part of a three-mem ring
        if (ring_info->isAtomInRingOfSize(bond_bgn_idx,3) || ring_info->isAtomInRingOfSize(bond_end_idx,3))
        {
            continue;
        }

        // Check if either atom in the bond is part of a four-mem ring
        if (ring_info->isAtomInRingOfSize(bond_bgn_idx,4) || ring_info->isAtomInRingOfSize(bond_end_idx,4))
        {
            // Can only have one unsaturation in this four-mem ring, exo double bond
            if (!can_unsat_four_ring(mol, bond_bgn_idx, bond_end_idx, non_sp3_indices))
            {
                continue;
            }
        }

        // Add bond index to viable bond indices
        viable_bond_indices.emplace(bond_idx);
    }

    // Designer tip: try to remove as many unviable bonds as you can before the combination generator
    // A lower number of viable bonds will result in a significant reduction in possible combinations
    // and thus a faster code iterating through these loops!

    // Get the number of viable bonds
    const unsigned int num_viable_bonds = static_cast<unsigned int>(viable_bond_indices.size());

    // Unsaturate 1 bond, to maximum number of unsaturations (all putative viable)
    for (unsigned int num_unsat = 1; num_unsat <= num_viable_bonds; ++num_unsat)
    {
        // Store all viable unsaturation operations
        std::unordered_set<int> viable_operations = {};

        // Avoid doing the same operations, if the bonds or atoms are symmetrically identical
        std::set<std::vector<int>> seen_unsaturation_patterns = {};

        // Generate all combinations to unsaturate input molecule N times
        const std::vector<std::vector<int>> combinations = generateCombinations(viable_bond_indices, num_unsat); // change into unsigned ints

        // Iterate over each combination of unsaturations
        for (const std::vector<int>& unsaturation_pattern : combinations)
        {    
            // See if this unsaturation pattern has been evaluated before in terms of symmetrical atoms
            std::vector<int> standard_unsaturation_pattern = standardize_operation(unsaturation_pattern, equivalent_bonds);

            // Store the unsaturation pattern
            if (seen_unsaturation_patterns.count(standard_unsaturation_pattern) == 0)
            {
                seen_unsaturation_patterns.emplace(standard_unsaturation_pattern);
            }

            // An equivalent operation has been done before, skip loop and save time
            else
            {
                ; //continue; Fix this so the loop is avoided when a similar fingerprint has been seen
            }

            // Assume the operation is viable
            bool is_viable_operation = true;

            // Do a bond-adjacency check, no allenes be formed!
            bool contains_adjacent_unsaturations = false;

            // Number of bond unsaturations in this pattern
            const int num_unsaturations = static_cast<int>(unsaturation_pattern.size());

            // Analyse the upper triangular submatrix of unsaturations (i,j)
            // Iterate over the indices of bonds (i) that will be unsaturated 
            for (int i = 0; i < num_unsaturations; ++i)
            {    
                // Get the adjacent bond indices of the bond with index i 
                const std::unordered_set<unsigned int> adjacent_bonds = bond_adjacencies[unsaturation_pattern[i]];

                // Iterate over the indices of atoms (j) that will be unsaturated
                for (int j = i + 1; j < num_unsaturations; ++j)
                {
                    // Check if the adjacent bond will be unsaturated as well
                    if (adjacent_bonds.count(unsaturation_pattern[j]) > 0)
                    {
                        contains_adjacent_unsaturations = true;
                        break;
                    }
                }

                // If the unsaturation pattern contains two adjacent bonds
                if (contains_adjacent_unsaturations)
                {
                    // Stop the adjacency check, break loops
                    break;
                }
            }

            // If the pattern contains two adjacent bonds
            if (contains_adjacent_unsaturations)
            {
                // Skip to the next pattern
                continue;
            }

            // If the unsaturation operation is viable
            if (is_viable_operation)
            {
                // Create a work molecule that we can unsaturate by changing bondtypes
                std::shared_ptr<RDKit::RWMol> unsaturate_mol(new RDKit::RWMol(mol));

                // Create a set to store atom indices that have become non-sp3 hybridized
                std::unordered_set<unsigned int> unsat_non_sp3_indices = non_sp3_indices;

                // Counter to restrict triple bonds in rings, there are zero triple bonds
                // at the start, they get introduced at the last stage
                unsigned int num_trip_ring_bonds = 0; // Should be a map with the ring's index

                // Counter to restrict the number of Bredt's rule violations
                unsigned int num_bredt_violations = 0;

                // Iterate over all unsaturation operations
                for (const int& unsat_bond_idx : unsaturation_pattern)
                {    
                    // Assume we can unsaturate the current bond
                    bool can_unsaturate = true;

                    // std:pair adjustment TODO
                    const int unsat_bgn_idx = bond_information[unsat_bond_idx].first;
                    const int unsat_end_idx = bond_information[unsat_bond_idx].second;

                    // Get the bond to unsaturate
                    RDKit::Bond *unsat_bond = unsaturate_mol->getBondWithIdx(unsat_bond_idx);

                    // Check if we can introduce triple bond in ring
                    // Minimal ringsize for triple bond has to be 9
                    // When multiple triple bonds are installed, the ring needs to be bigger
                    if (new_bond_order == 3 && ring_info->numBondRings(unsat_bond_idx))
                    {
                        // Check the minimal ring size of this bond 
                        const unsigned int min_ring_size = ring_info->minBondRingSize(unsat_bond_idx);

                        // If ringsize < 9 ; can have 0 triple bonds
                        if (min_ring_size < 9)
                        {
                            can_unsaturate = false;
                            break;
                        } 

                        // If 9 <= ringsize < 11 
                        else if (9 <= min_ring_size && min_ring_size < 11)
                        {
                            // Can have 1 triple bond
                            if (num_trip_ring_bonds >= 1)
                            {
                                can_unsaturate = false;
                                break;
                            }

                            // Increment number of triple bonds in this ring
                            else
                            {
                                ++num_trip_ring_bonds;
                            }
                        }

                        // If 11 <= ringsize
                        else
                        {
                            // Can have 2 triple bonds
                            if (num_trip_ring_bonds >= 2)
                            {
                                can_unsaturate = false;
                                break;
                            }

                            // Increment number of triple bonds in this ring
                            else
                            {
                                ++num_trip_ring_bonds;
                            }
                        }
                    }

                    // If an Sp2 was introduced in this operation
                    // Check if we can still add a new Sp2
                    if (ring_info->isAtomInRingOfSize(unsat_bgn_idx,4) || ring_info->isAtomInRingOfSize(unsat_end_idx,4))
                    {
                        // Can only have one unsaturation in this four-mem ring, exo double bond
                        if (!can_unsat_four_ring(*unsaturate_mol, unsat_bgn_idx, unsat_end_idx, unsat_non_sp3_indices)) // the problem is here, you iterate over rings from *mol, not unsaturate_mol, mol has no sp2's at
                        {
                            can_unsaturate = false;
                            break;
                        }
                    }

                    // Bredt's Rule violations, no double bonds with bridgehead atoms
                    // with some exceptions depending on ringsizes
                    // If either atom in the bond is a bridgehead atom
                    if (bridgeheads.count(unsat_bgn_idx) > 0 || bridgeheads.count(unsat_end_idx) > 0)
                    {
                        // For the current number of Bredt's rule violations, can we unsaturate this bond
                        if (!can_unsat_bridgehead(bridgeheads, cycles, num_bredt_violations))
                        {
                            can_unsaturate = false;
                            break;
                        }
                        // If we can, we increase the number of violations counter
                        else
                        {
                            ++num_bredt_violations;
                        }
                    }
                
                    // If we can unsaturate this bond
                    if (can_unsaturate)
                    {
                        // Change the bond order of this bond
                        unsat_bond->setBondType(bond_types[new_bond_order]);
                        
                        // Add bond's bgn and end atoms as non-sp3 atoms
                        unsat_non_sp3_indices.emplace(unsat_bgn_idx);
                        unsat_non_sp3_indices.emplace(unsat_end_idx);
                    }
                }

                // Helps fix implicit/explicit hydrogens
                RDKit::MolOps::sanitizeMol(*unsaturate_mol);

                // Convert the molecule back to a regular RDKit mol obj
                RDKit::ROMol sanitized_mol(*unsaturate_mol);

                // Convert the unsaturated molecule to its corresponding SMILES string
                const std::string unsaturate_smiles = RDKit::MolToSmiles(sanitized_mol);

                // Store the SMILES string 
                unique_smiles.emplace(unsaturate_smiles);
            }
        }
    }

    // Return the unique SMILES strings
    return unique_smiles;
}

// Driver function
int main(int argc, char *argv[])
{
    // File names
    std::string in_file_path;
    std::string out_file_path;

    // Argument parser
    int opt;
    while ((opt = getopt(argc, argv, "i:o:h")) != -1)
    {
        switch (opt)
        {
            case 'i':
                in_file_path = optarg;
                break;
            case 'o':
                out_file_path = optarg;
                break;
            case 'h':
                std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path>" << std::endl;
                return 0;
            default:
                std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path>" << std::endl;
                return 1;
        }
    }

    // Verify path
    if (in_file_path.empty())
    {
        std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path>" << std::endl;
        return 1;
    }

    // Verify path
    if (out_file_path.empty())
    {
        std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path>" << std::endl;
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

    // Ensure uniqueness of the output via string hashing
    std::unordered_set<size_t> seen_hashes;

    // Create a hashing object
    std::hash<std::string> stringHasher;

    // Read the input file line by line
    while (std::getline(input_stream, line))
    {
        // Turn the line (SMILES) into a molecule
        std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(line));
        
        // Generate all double-bond unsaturations for the current input molecule 
        std::unordered_set<std::string> double_bond_smiles = generate_unsaturations(*mol, 2);

        // For each modified molecule, generate triple-bond unsaturations
        for (const std::string& double_smiles : double_bond_smiles)
        {
            // Turn the unsaturated SMILES into a molecule
            std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(double_smiles));

            // Generate all triple-bond unsaturations for the current input molecule
            std::unordered_set<std::string> triple_bond_smiles = generate_unsaturations(*mol, 3);

            // For each modified molecule, generate triple-bond unsaturations
            for (const std::string& triple_smiles : triple_bond_smiles)
            {
                // Hash the SMILES string
                size_t hashed_smiles = stringHasher(triple_smiles);

                // If you have not seen this hash
                if (seen_hashes.find(hashed_smiles) == seen_hashes.end()) 
                {
                    // Write the unsaturated SMILES string to output file
                    output_stream << triple_smiles + "\n";

                    // Store the hashed SMILES string
                    seen_hashes.insert(hashed_smiles);
                }
            }
        }
    }

    // Close the file streams
    input_stream.close();
    output_stream.close();

    // Signal success 
    return 0;
}