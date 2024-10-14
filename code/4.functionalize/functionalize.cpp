/*
    Author: Andreas Luttens
    Contact: andreas.luttens@gmail.com
    Date: June 22, 2023
    Description: Introduce heteroatoms (N,O) in hydrocarbons
*/

// Basic utilities
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>
#include <functional> 
#include <utility>
#include <unistd.h> 

// RDKit functions
#include <GraphMol/GraphMol.h>
#include <GraphMol/Canon.h> // Check for symmetry later
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <GraphMol/MolStandardize/Normalize.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

// A map between atomic symbols and atomic numbers
std::unordered_map<std::string, unsigned int> atomic_numbers_map = {{"N", 7}, {"O", 8}};

// A map between atomic symbols and number of hydrogens needed to change a C to this element
std::unordered_map<std::string, unsigned int> min_hydrogen_map = {{"N", 1}, {"O", 2}};

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

// A function to identify bonds that are equivalent in terms of their composing atoms' CIP rank
std::unordered_map<unsigned int,std::vector<unsigned int>> get_equivalent_bonds(const RDKit::ROMol& mol, std::unordered_map<unsigned int,unsigned int>& atom_symmetries){

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

// A function to identify atoms that are equivalent in terms of their CIP rank
std::unordered_map<unsigned int,std::vector<unsigned int>> get_equivalent_atoms(const RDKit::ROMol& mol, std::unordered_map<unsigned int,unsigned int>& atom_symmetries)
{
    //
    std::map<unsigned int,std::vector<unsigned int>> symmetries = {};

    // Iterate over all atoms in the molecule
    unsigned int num_atoms = mol.getNumAtoms();
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx)
    {
        // Get the atom symmetry
        const unsigned int atom_sym = atom_symmetries[atom_idx];

        // If this atom symmetry has not been seen before
        if (symmetries.count(atom_sym) == 0)
        {
            // Add the symmetry, and an empty set to store pairs later
            symmetries.emplace(atom_sym, std::vector<unsigned int>{});
        }

        // Add the atom symmetry and corresponding atom index
        symmetries[atom_sym].emplace_back(atom_idx);                
    }
    
    // Equivalent atoms should be sorted!
    std::unordered_map<unsigned int, std::vector<unsigned int>> equivalent_atoms;
    // Iterate over all bond symmetries
    for (const auto& sym_pair : symmetries)
    {
        //
        std::vector<unsigned int> sym_atoms = sym_pair.second;
        std::sort(sym_atoms.begin(),sym_atoms.end());

        // Iterate over all atoms with this bond symmetry
        for (const auto& sym_bond: sym_atoms)
        {
            equivalent_atoms.emplace(sym_bond, sym_atoms);
        }
    }

    //
    return equivalent_atoms;
}

// Function retrieve information on all carbons (bound to heteroatoms, or bound to sp-hyb atoms)
void get_bound_carbons_information(const RDKit::ROMol& mol, std::unordered_set<int>& hetero_bound_carbons, std::unordered_set<int>& sphyb_bound_carbons)
{
    // Get the number of atoms in this molecule
    const unsigned int num_atoms = mol.getNumAtoms();
    
    // Iterate over all atoms in the molecule
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx) 
    {
        // Retrieve the atom with this index
        const RDKit::Atom *atom = mol.getAtomWithIdx(atom_idx);

        // Verify the atom is a carbon
        const unsigned int elem_num = atom->getAtomicNum();

        // If it isn't a carbon
        if (elem_num != 6)
        {
            continue;
        }
        
        // Iterate over all neighbor atoms of this carbon
        for(const auto &neighbor_atom: mol.atomNeighbors(atom))
        {
            // Retrieve the atomic number of this neighbor atom
            const unsigned int neighbor_elem_num = neighbor_atom->getAtomicNum();

            // Check if the neighbor atom is an oxygen or nitrogen
            if (neighbor_elem_num == 7 || neighbor_elem_num == 8)
            {
                hetero_bound_carbons.emplace(atom_idx);
            }

            // Retrieve the hybridization type of the neighbor atom
            const RDKit::Atom::HybridizationType hybridization = neighbor_atom->getHybridization();

            // Check if the neighbor atom is sp hybridized
            if (hybridization == RDKit::Atom::HybridizationType::SP)
            {
                sphyb_bound_carbons.emplace(atom_idx);
            }
        }
    }
}

//
std::unordered_map<unsigned int, unsigned int> get_bound_hydrogens(const RDKit::ROMol& mol)
{
    std::unordered_map<unsigned int,unsigned int> hydrogen_map = {};    

    const unsigned int num_atoms = mol.getNumAtoms();
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx) 
    {
        const RDKit::Atom *atom = mol.getAtomWithIdx(atom_idx);

        // Get the hybridization type
        const unsigned int h_count = atom->getTotalNumHs();

        hydrogen_map.emplace(atom_idx, h_count);
    }

    return hydrogen_map;
}

//
std::pair<std::unordered_set<unsigned int>,std::vector<std::unordered_set<unsigned int>>> find_bridgehead_atoms(RDKit::ROMol& mol)
{
    // Set containing all bridgehead atom indices
    std::unordered_set<unsigned int> bridgehead_indices = {};

    // SSSR indices
    const std::vector<RDKit::RingInfo::INT_VECT> sssr = mol.getRingInfo()->atomRings();

    // pattern to match sp3 hybridized ring atoms with 3 or 4 neighbors
    const RDKit::RWMol *pattern = RDKit::SmartsToMol("[x3,x4:1]"); // can move this outside of function, its const
    
    std::vector<RDKit::MatchVectType> matches;    

    // Perform substructure match
    RDKit::SubstructMatch(mol, *pattern, matches);
    
    //
    std::vector<unsigned int> standard_matches = {};

    // Loop over the matches
    for(size_t i = 0 ; i < matches.size() ; ++i) 
    {
        standard_matches.emplace_back(matches[i][0].second);
    }

    // Get the number of rings
    const unsigned int num_rings = static_cast<unsigned int>(sssr.size());

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

// Standardize the given operation based on equivalent bonds
std::vector<int> standardize_operation(const std::vector<int>& operation, std::unordered_map<unsigned int,std::vector<unsigned int>>& equivalent_atoms){

    std::vector<int> standard_operation;

    // Iterate over all atom indices that will be mutated
    for (const int& atom_idx : operation) { 

        // Get the minimal equivalent atom index, since the vector is sorted, it is the first element
        // if there were no equivalent atoms, it will return itself as sole (and first) element
        int min_equivalent_idx = equivalent_atoms[atom_idx].front();

        //
        standard_operation.emplace_back(min_equivalent_idx);
    }

    // sort the standardized fingerprint
    sort(standard_operation.begin(), standard_operation.end()); //, stepwise_sorter);

    return standard_operation;
}

// Function to check if heteroatom insertion is allowed based on Bredt's rule violation
bool can_bredt_amide(const RDKit::ROMol& mol, const int& atom_idx){

    // Get the minimal ring size
    const unsigned int min_ring_size = mol.getRingInfo()->minAtomRingSize(atom_idx);

    // Fetch the bridgehead nitrogen
    const RDKit::Atom *bridgehead = mol.getAtomWithIdx(atom_idx);

    // Iterate over all neighbor atoms of the bridgehead nitrogen
    for (const auto &neighbor_atom: mol.atomNeighbors(bridgehead)) {

        // Get the hybridization type
        const RDKit::Atom::HybridizationType hybridization = neighbor_atom->getHybridization();

        // Check if the neighbor atom is sp2 hybridized and not aromatic
        if (hybridization == RDKit::Atom::HybridizationType::SP2 && !neighbor_atom->getIsAromatic())
        {
            // Have to check the ring size
            if (min_ring_size < 9)
            {
                return false;
            }
        }
    }

    return true;
}

// Function check if the atom in a molecule is a carbon and if the minimum number of bound hydrogens matches the parsed number
bool is_viable_atom(const RDKit::ROMol& mol, const unsigned int& atom_idx, const unsigned int& min_hydrogens,
                    const std::unordered_set<int>& hetero_bound_carbons,
                    const std::unordered_set<int>& sphyb_bound_carbons, 
                    std::unordered_map<unsigned int,std::unordered_set<unsigned int>>& atom_adjacencies, 
                    const std::unordered_set<unsigned int>& bridgeheads, const bool& verbose)
{

    // Fetch the atom with the parsed index
    const RDKit::Atom *atom = mol.getAtomWithIdx(atom_idx);

    // Verify the atom is a carbon
    const unsigned int elem_num = atom->getAtomicNum();
    if (elem_num != 6)
    {
        return false;
    }
    
    // Verify the would-be valence is possible with the number of bound hydrogens
    if (atom->getTotalNumHs() < min_hydrogens)
    {
        return false;
    }

    // The carbon cannot be bound to a heteroatom or a sp-hybridized atom (hetero-hetero etc)
    if (hetero_bound_carbons.find(atom_idx) != hetero_bound_carbons.end() || sphyb_bound_carbons.find(atom_idx) != sphyb_bound_carbons.end())
    { 
        return false;
    }

    // Here we check if a N-C-O is about to be formed // Collapse all of this in one function
    // We're introducing oxygens
    if (min_hydrogens == 2u)
    {
        // Iterate over all neighbors
        for (const unsigned int& neighbor_idx : atom_adjacencies[atom_idx])
        {
            // Check if the neighbor is bound to a hetero atom
            if (hetero_bound_carbons.find(neighbor_idx) != hetero_bound_carbons.end())
            {
                // Check the hybridization state of this neighbor atom
                const RDKit::Atom::HybridizationType hybridization = mol.getAtomWithIdx(neighbor_idx)->getHybridization();

                // This will be an invalid pattern
                if (hybridization == RDKit::Atom::HybridizationType::SP3)
                {
                    return false;
                }
            }
        }
    }

    // Bredt's rule for twisted amides
    if (min_hydrogens == 1 && bridgeheads.find(atom_idx) != bridgeheads.end())
    {
        // Bridgehead amides: If a “non-zero” bridgehead nitrogen is bound to a nonaromatic sp2 atom, 
        // the ring sizes of the smallest set of smallest rings will be checked. 
        // At least one ring of the nitrogen must be ≥9
        if (!can_bredt_amide(mol, atom_idx))
        {
            // Skip this atom
            return false;
        }   
    }
    
    return true;
}

// Function to read file with chemical patterns (SMARTS) to filter out unwanted molecules
void read_filter_file(const std::string& filter_file, std::vector<RDKit::RWMol *>& filters)
{
    // Verify you can open the provided file
    std::ifstream filter_stream(filter_file);
    if (!filter_stream) 
    {
        // Couldn't open, return error message and exit
        std::cerr << "Error opening file: " << filter_file << std::endl;
        return;
    }

    // String to line from file
    std::string line;

    // Read in the filter file line by line
    while (std::getline(filter_stream, line)) 
    {
        // Read line into string stream
        std::istringstream iss(line);

        // Two strings, one for the SMARTS pattern, and one for its description
        std::string rule, smarts;

        // Assuming columns are separated by tabs ('\t')
        if (iss >> rule >> smarts) 
        {
            // Turn the SMARTS pattern into a RDKit Mol
            RDKit::RWMol *pattern = RDKit::SmartsToMol(smarts.c_str());
             
            // Store pattern into vector of RDKit mols
            filters.emplace_back(pattern);
        }

        // Could not parse this SMARTS, return error message
        else 
        {
            std::cout << "Error parsing SMARTS pattern for: " << rule << std::endl;
        }
    }

    // Close the stream
    filter_stream.close();
}

// Function to generate all mutations for a given input molecule and heteroatom
std::unordered_set<std::string> generate_mutations(RDKit::ROMol& mol, const std::string element, bool verbose) 
{    
    // Get the atomic number for the element we want to insert
    const unsigned int atom_num = atomic_numbers_map[element];

    // Get the minimal number of hydrogens we need for a carbon to
    // mutate into the element we want to insert
    unsigned int min_hydrogen = min_hydrogen_map[element];

    // A set of unique SMILES strings for this hydrocarbon
    std::unordered_set<std::string> unique_smiles;

    // Add the input smiles
    unique_smiles.emplace(RDKit::MolToSmiles(mol));

    // Perform ring analysis
    RDKit::MolOps::findSSSR(mol);

    // Store ring info, can I do a const here? This molecule won't change
    RDKit::RingInfo* ring_info = mol.getRingInfo();

    // Perform hidden symmetry analysis by assigning labels to atoms
    // cleanIt=True,force=True,flagPossibleStereoCenters=True
    RDKit::MolOps::assignStereochemistry(mol,true,true,true);

    // Perform symmetry analysis
    std::unordered_map<unsigned int,unsigned int> atom_symmetries = get_atom_symmetries(mol);
    
    // Get all equivalent atoms
    std::unordered_map<unsigned int,std::vector<unsigned int>> equivalent_atoms = get_equivalent_atoms(mol, atom_symmetries);

    // Find all adjacent heavy atoms
    std::unordered_map<unsigned int,std::unordered_set<unsigned int>> atom_adjacencies = get_atom_adjacencies(mol);
    
    // A set to store indices of carbons that are already bound to heteroatoms
    std::unordered_set<int> hetero_bound_carbons = {};

    // A set to store indices of carbons that are bound to sp hybrized atoms
    std::unordered_set<int> sphyb_bound_carbons = {};

    // Find all carbons that are bound to a hetero-atom, or an sp-hybridized carbon
    get_bound_carbons_information(mol, hetero_bound_carbons, sphyb_bound_carbons);

    // Create a bound hydrogen count map
    std::unordered_map<unsigned int, unsigned int> hydrogen_map = get_bound_hydrogens(mol);

    // Find the bridgehead atoms and list all rings in the molecule
    std::pair<std::unordered_set<unsigned int>,std::vector<std::unordered_set<unsigned int>>> information = find_bridgehead_atoms(mol);
    std::unordered_set<unsigned int> bridgeheads = information.first;
    std::vector<std::unordered_set<unsigned int>> cycles = information.second;

    // A store to store all atom indices of all (putatively valid) carbon atoms
    std::unordered_set<unsigned int> carbon_indices = {};

    // Get the number of atoms in this molecule
    const unsigned int num_atoms = mol.getNumAtoms();

    // Iterate over all atoms in the molecule
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx)
    {
       
        // Has to be a mutatable atom
        if (!is_viable_atom(mol, atom_idx, min_hydrogen, hetero_bound_carbons,  
                            sphyb_bound_carbons, atom_adjacencies, bridgeheads, verbose))
        {            
            continue;
        }

        // Store this atom's index as viable
        carbon_indices.emplace(atom_idx);
    }

    // Get the number of carbon atoms to mutate
    const unsigned int num_mutations = static_cast<unsigned int>(carbon_indices.size());

    // Mutate 1 to N carbon atoms in this molecule
    for (unsigned int n = 1; n <= num_mutations; ++n)
    {
        // Avoid doing the same operations, if the atoms are symmetrically identical
        std::set<std::vector<int>> seen_mutation_patterns = {};

        // Generate all combinations to mutate input molecule in N different spots
        const std::vector<std::vector<int>> combinations = generateCombinations(carbon_indices, n); 

        // Iterate over each combination of mutations
        for (const std::vector<int>& mutation_pattern : combinations)
        {    
            // See if this mutation pattern has been evaluated before in terms of symmetrical atoms
            std::vector<int> standard_mutation_pattern = standardize_operation(mutation_pattern, equivalent_atoms);

            // Store the mutation pattern
            if (seen_mutation_patterns.count(standard_mutation_pattern) == 0)
            {
                seen_mutation_patterns.emplace(standard_mutation_pattern);
            }

            // an equivalent operation has been done before, skip loop and save time
            // do this before any evaluation of valid operations
            else
            {   
                ; //continue
            }

            // Do a atom-adjacency check, no hetero_hetero to be formed!
            bool contains_adjacent_mutations = false;

            // Do a common neighbor check, no aminals etc to be formed!
            bool invalid_sharing_neighbor = false;

            // Number of atom mutations in this pattern
            const int num_mutations = static_cast<int>(mutation_pattern.size());

            // Analyse the upper triangular submatrix of mutations (i,j)
            // Iterate over the indices of atom (i) that will be mutated 
            for (int i = 0; i < num_mutations; ++i) 
            {    
                // Get the adjacent atom indices of the atom with index i 
                const std::unordered_set<unsigned int> adjacent_atoms = atom_adjacencies[mutation_pattern[i]];

                // Iterate over the indices of atoms (j) that will be mutated
                for (int j = i + 1; j < num_mutations; ++j) 
                {
                    // Check if the mutations are adjacent
                    if (adjacent_atoms.count(mutation_pattern[j]) > 0)
                    {
                        contains_adjacent_mutations = true;
                        break;
                    }

                    // Check if they share a common neighbor
                    for (const unsigned int& k : atom_adjacencies[mutation_pattern[j]])
                    {
                        // Share a common neighbor
                        if (adjacent_atoms.count(k) > 0)
                        {
                            // Get the hybridization state of the common neighbor
                            RDKit::Atom::HybridizationType hybridization = mol.getAtomWithIdx(k)->getHybridization();

                            // If the common neighbor is sp3 hybridized
                            if (hybridization == RDKit::Atom::HybridizationType::SP3)
                            {
                                // Check if they are in the same ring
                                if (ring_info->areAtomsInSameRing(mutation_pattern[i], mutation_pattern[j]))
                                {
                                    // Can't have N-@Csp3-@N
                                    if (element == "N")
                                    {
                                        invalid_sharing_neighbor = true;
                                        break;
                                    }

                                    // Can have O-@Csp3-@O
                                    else
                                    {
                                        // If the ringsize is larger than 4
                                        if (ring_info->areAtomsInSameRingOfSize(mutation_pattern[i], mutation_pattern[j],4u))
                                        {
                                            invalid_sharing_neighbor = true;
                                            break;
                                        }
                                    }
                                }

                                // They're not in the same ring
                                else
                                {
                                    invalid_sharing_neighbor = true;
                                    break;
                                }              
                            }
                        }
                    
                    }

                    // If they're sharing an invalid neighbor, break loop
                    if (invalid_sharing_neighbor)
                    {
                        break;
                    }
                }

                // if the pattern contains two adjacent atoms
                // OR
                // there is an invalid XCX pattern
                if (contains_adjacent_mutations || invalid_sharing_neighbor)
                {
                    // Skip this mutation pattern
                    continue;
                }
            }

            // if the pattern contains two adjacent atoms
            if (contains_adjacent_mutations || invalid_sharing_neighbor)
            {
                // Skip to the next pattern   
                continue;
            }

            // Create a work molecule that we can mutate by changing atomic numbers
            std::unique_ptr<RDKit::RWMol> mutate_mol(new RDKit::RWMol(mol));

            // Iterate over all atoms to mutate
            for (const int& mut_atom_idx : mutation_pattern)
            {
                // Get the atom that will be mutated
                RDKit::Atom *mut_atom = mutate_mol->getAtomWithIdx(mut_atom_idx);

                // Change the atomic number of this carbon into the mutant
                mut_atom->setAtomicNum(atom_num);
            }

            // Helps fix implicit/explicit hydrogens
            RDKit::MolOps::sanitizeMol(*mutate_mol); 

            // Initialize an enumerator of tautomers for this functionalized molecule
            RDKit::MolStandardize::TautomerEnumerator enumerator;
            
            // Standardize the molecule into a single tautomer
            RDKit::ROMol *standard_tautomer_mol = enumerator.canonicalize(RDKit::ROMol(*mutate_mol)); 

            // Convert the functionalized molecule to its corresponding SMILES string
            const std::string mutate_smiles = RDKit::MolToSmiles(*standard_tautomer_mol);

            // Store the SMILES string
            unique_smiles.emplace(mutate_smiles);
        }
    }

    // Return the unique SMILES strings
    return unique_smiles;
}

// Check if the given molecule contains any of the unwanted chemical patterns
bool functional_group_filter(const RDKit::ROMol& mol, std::vector<RDKit::RWMol *>& filters)
{
    // Iterate over all patterns in the filter
    for (const auto& pattern : filters)
    {
        // To assign pattern matches
        RDKit::MatchVectType hits;

        // If the mutate mol hits a pattern
        if(RDKit::SubstructMatch(mol, *pattern, hits)) 
        {
            return false;
        }
    }
    return true;
}

// Driver function
int main(int argc, char *argv[])
{
    
    // File names
    std::string in_file_path;
    std::string out_file_path;
    std::string filter_file_path;

    // Boolean to control verbose output
    bool verbose = false;

    // Argument parser
    int opt;
    while ((opt = getopt(argc, argv, "i:o:f:v:h")) != -1)
    {
        switch (opt)
        {
            case 'i':
                in_file_path = optarg;
                break;
            case 'o':
                out_file_path = optarg;
                break;
            case 'f':
                filter_file_path = optarg;
                break;
            case 'v':
                verbose = static_cast<bool>(std::stoi(optarg));
                break;
            case 'h':
                std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -f <filter_file_path>" << std::endl;
                return 0;
            default:
                std::cerr << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -f <filter_file_path>" << std::endl;
                return 1;
        }
    }

    // Verify path
    if (in_file_path.empty())
    {
        std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -f <filter_file_path>" << std::endl;
        return 1;
    }

    // Verify path
    if (out_file_path.empty())
    {
        std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -f <filter_file_path>" << std::endl;
        return 1;
    }

    // Verify path
    if (filter_file_path.empty())
    {
        std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -f <filter_file_path>" << std::endl;
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

    // Create a SMARTS pattern filter array
    // from cheapest to most expensive calculation
    std::vector<RDKit::RWMol*> filters;

    // Read in the filter file
    read_filter_file(filter_file_path, filters);

    // String for reading lines, there are some bugs with regular openeye file opening and reading
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
        
        // Combinatorially introduce nitrogen atoms for the current input molecule 
        std::unordered_set<std::string> nitrogen_mutant_smiles = generate_mutations(*mol, "N", verbose);

        // Iterate over all canonical tautomers of nitrogenated mutants
        for (const std::string& nitrogen_smiles : nitrogen_mutant_smiles)
        { 
            // Turn the mutated SMILES into a molecule
            std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(nitrogen_smiles));
            
            // Combinatorially introduce oxygen atoms into the molecule
            std::unordered_set<std::string> oxygen_mutant_smiles = generate_mutations(*mol, "O", verbose);

            // Iterate over all canonical tautomers of oxygenated mutants
            for (const std::string& oxygen_smiles : oxygen_mutant_smiles)
            { 
                // Hash the SMILES string
                size_t hashed_smiles = stringHasher(oxygen_smiles);

                // If you have not seen this hash
                if (seen_hashes.find(hashed_smiles) == seen_hashes.end()) 
                {
                    // Turn the mutated SMILES into a molecule
                    std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(oxygen_smiles));

                    // If the molecule survives all chemical pattern filters
                    if (functional_group_filter(*mol, filters))
                    {
                        // Write the functionalized SMILES string to output file
                        output_stream << oxygen_smiles + "\n";

                        // Store the hashed SMILES string
                        seen_hashes.insert(hashed_smiles);
                    }
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