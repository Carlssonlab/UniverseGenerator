/*
    Author: Andreas Luttens
    Contact: andreas.luttens@gmail.com
    Date: June 22, 2023
    Description: Generate all combinations you can activate a scaffold for superstructure generation
*/

// Basic utilities
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <functional> 
#include <vector>
#include <utility>
#include <unistd.h> 

// RDKit functions
#include <GraphMol/GraphMol.h>
#include <GraphMol/Canon.h> // Check for symmetry later
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

// A map between bondorders and bondtypes
std::unordered_map<unsigned int, RDKit::Bond::BondType> bond_types = {{1u, RDKit::Bond::SINGLE},
                                                                      {2u, RDKit::Bond::DOUBLE},
                                                                      {3u, RDKit::Bond::TRIPLE}};

// A map between bondtypes and bondorders
std::unordered_map<RDKit::Bond::BondType, unsigned int> rev_bond_types = {{RDKit::Bond::SINGLE, 1u},
                                                                          {RDKit::Bond::DOUBLE, 2u},
                                                                          {RDKit::Bond::TRIPLE, 3u}};

//
void generateCombinationsHelper(const std::unordered_set<unsigned int>& nums, int k, int start, std::vector<int>& currentCombination, std::vector<std::vector<int>>& combinations) {
    if (static_cast<int>(currentCombination.size()) == k) {
        combinations.push_back(currentCombination);
        return;
    }

    for (auto it = nums.begin(); it != nums.end(); ++it) {
        int num = *it;
        if (num >= start) {
            currentCombination.push_back(num);
            generateCombinationsHelper(nums, k, num + 1, currentCombination, combinations);
            currentCombination.pop_back();
        }
    }
}

//
std::vector<std::vector<int>> generateCombinations(const std::unordered_set<unsigned int>& nums, int k) {
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

// Function to map how many hydrogens are bound to each heavy atom
std::unordered_map<unsigned int,unsigned int> get_hydrogen_map(RDKit::ROMol& mol)
{
    // A map between heavy atoms and their number of bound hydrogens
    std::unordered_map<unsigned int,unsigned int> hydrogen_map = {};

    // Iterate over all heavy atoms
    for (const auto& atom : mol.atoms())
    {
        // Counter
        unsigned int num_hydrogens = 0;

        // Iterate over the neighbors of this atom
        for (const auto& neigbor : mol.atomNeighbors(atom))
        {
            // If it is a hydrogen
            if (neigbor->getAtomicNum() == 1u)
            {
                // Increment the counter
                ++num_hydrogens;
            }
        }

        // Map the number of hydrogens of this heavy atom
        hydrogen_map.emplace(atom->getIdx(), num_hydrogens);
    }
    return hydrogen_map;
}

// Count the number of activations of the same type
unsigned int count_num_same_activations(const RDKit::ROMol& mol, const unsigned int atomic_number, const unsigned int degeneracy)
{
    // Counter
    unsigned int num_same_activations = 0;
    
    // Get the number of atoms in this molecule
    const unsigned int num_atoms = mol.getNumAtoms();

    // Iterate over all molecules in the molecule
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx)
    {
        // If the atom has the queried atomic number
        if (mol.getAtomWithIdx(atom_idx)->getAtomicNum() == atomic_number) 
        {
            // Increment the counter
            ++num_same_activations;
        }
    }

    // If the activation type is bond (fusion), half the count of activations 
    // as there are two tags per bond activation
    // this number is the degeneracy
    
    return num_same_activations/degeneracy;
}

// Introduce activation tags for (single, double, or triple) bonds on all viable hydrogens
std::unordered_set<std::string> generate_coupling_activations(const RDKit::ROMol& mol, unsigned int new_bond_order, 
                                                              std::unordered_map<unsigned int,std::unordered_set<unsigned int>>& atom_adjacencies, 
                                                              std::unordered_map<unsigned int,unsigned int>& hydrogen_map, 
                                                              const bool verbose) 
{    
    // A set of unique SMILES strings
    std::unordered_set<std::string> unique_smiles = {};
    
    // A set containing atom indices of valid hydrogens
    std::unordered_set<unsigned int> valid_hydrogens = {};
   
    // Count the number of times the current activation has been introduced
    // and adjust the isotope number of the activation tag accordingly
    const unsigned int num_same_activations = count_num_same_activations(mol, 70u + new_bond_order, 1u);

    // Get the number of atoms in this molecule
    const unsigned int num_atoms = mol.getNumAtoms();

    // Iterate over all atoms in the molecule
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx)
    {
        // Has to be a hydrogen
        if (mol.getAtomWithIdx(atom_idx)->getAtomicNum() != 1u) // would it be faster it were heavy atoms?
        {
            continue;
        }

        // Adjacent heavy atom needs to have enough hydrogens
        if (hydrogen_map[*atom_adjacencies[atom_idx].begin()] < new_bond_order)
        {
            continue;
        }

        // Store this as a valid hydrogen
        valid_hydrogens.emplace(atom_idx);
    }

    // Get the number of hydrogens to activate
    const unsigned int num_activations = static_cast<unsigned int>(valid_hydrogens.size());

    // Activate 1 to N hydrogens
    for (unsigned int n = 1; n <= num_activations; ++n)
    {
        // Generate all combinations to activate input molecule N times
        const std::vector<std::vector<int>> combinations = generateCombinations(valid_hydrogens, n); 

        // Iterate over each combination of activations
        for (const std::vector<int>& activation_pattern : combinations) // 
        {

            // Symmetry breaker for speedup
            //TODO, is this needed?

            // Assume we can activate
            bool can_activate = true;

            // Iterate over all atoms in the molecule
            for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx)
            {
                // Has to be a heavy atom
                if (mol.getAtomWithIdx(atom_idx)->getAtomicNum() == 1u)
                {
                    continue;
                }

                // Keep track of how many activations this heavy atom gets
                unsigned int num_activations = 0;

                // Iterate over all neighboring atoms of the heavy atom
                for (const unsigned int& adj_idx : atom_adjacencies[atom_idx])
                {
                    num_activations += count(activation_pattern.begin(), activation_pattern.end(), adj_idx);
                }

                // Heavy atom needs to have enough hydrogens to allow n activations
                if (hydrogen_map[atom_idx] < num_activations * new_bond_order)
                {
                    can_activate = false;
                    break;
                }
            }

            // If this molecule cannot be activated in this way, proceed to next configuration
            if (!can_activate)
            {
                continue;
            }

            // Keep track of the heavy atoms that have been seen before
            std::unordered_map<unsigned int,unsigned int> seen_heavy = {};

            // Create a work mol object we can transform with activations
            std::shared_ptr<RDKit::RWMol> activate_mol(new RDKit::RWMol(mol));

            // Isotope counter
            unsigned int isotope_counter = 0;

            // Iterate over all activated positions
            for (const int& activate_idx : activation_pattern) // activate_idx should not be unique, possibly two activations on 1 atom
            {
                // Get the hydrogen atom inside the transform molecule
                RDKit::Atom *activate_hydrogen = activate_mol->getAtomWithIdx(activate_idx);

                // Get the index of the heavy atom the activated hydrogen is bound to
                const unsigned int heavy_idx = *atom_adjacencies[activate_idx].begin();

                // If the heavy atom has not been seen
                if (seen_heavy.find(heavy_idx) == seen_heavy.end())
                {
                    // Count it once
                    seen_heavy.emplace(heavy_idx, 1);
                }
                // It has been seen
                else
                {
                    // Increment the counter
                    ++seen_heavy[heavy_idx];
                }

                // Get the heavy atom that will carry the activation tag
                RDKit::Atom *activate_heavy = activate_mol->getAtomWithIdx(heavy_idx);

                // Get the bond between the hydrogen and heavy atoms
                RDKit::Bond *activate_bond = activate_mol->getBondBetweenAtoms(activate_idx, heavy_idx);

                // Change the hydrogen into the activated tag
                activate_hydrogen->setAtomicNum(70u + new_bond_order);

                // Change the isotope depending on how many same activations have been observed before
                activate_hydrogen->setIsotope(num_same_activations + isotope_counter);

                // Calculate the amount of new hydrogens to install
                int num_new_hydrogens = std::max(0, static_cast<int>(hydrogen_map[heavy_idx]) - static_cast<int>(seen_heavy[heavy_idx]*new_bond_order));

                // Change the number of explicit hydrogens on the activated heavy atom
                activate_heavy->setNumExplicitHs(num_new_hydrogens);
                
                //activate_heavy->setNumExplicitHs(0);

                activate_heavy->setNoImplicit(false);

                // Change the bond order between activated tag and heavy atom
                activate_bond->setBondType(bond_types[new_bond_order]);

                // Increase the isotope counter
                ++isotope_counter;
            }
        
            // Remove explicit hydrogens again
            RDKit::MolOps::removeHs(*activate_mol);

            // Helps fix hydrogens
            RDKit::MolOps::sanitizeMol(*activate_mol);

            // Convert the molecule back to a regular RDKit mol obj
            RDKit::ROMol sanitized_mol(*activate_mol);

            // Convert the activated molecule to its corresponding SMILES string
            const std::string activate_smiles = RDKit::MolToSmiles(sanitized_mol);

            // Add the SMILES to the set of activated molecules
            unique_smiles.emplace(activate_smiles);
        }
    }

    return unique_smiles;
}

// Introduce activation tags for fusion (single, no double right now) bonds on all viable hydrogens
std::unordered_set<std::string> generate_fusion_activations(const RDKit::ROMol& mol, unsigned int new_bond_order, 
                                                           std::unordered_map<unsigned int,std::unordered_set<unsigned int>>& bond_adjacencies, 
                                                           std::unordered_map<unsigned int, std::pair<unsigned int,unsigned int>>& bond_information,
                                                           std::unordered_map<unsigned int,unsigned int>& hydrogen_map, 
                                                           const bool verbose) 
{    
    // A set of unique SMILES strings
    std::unordered_set<std::string> unique_smiles = {};

    // A set containing atom indices of valid hydrogens
    std::unordered_set<unsigned int> valid_bond_indices = {};

    // Count the number of times the current activation has been introduced
    // and adjust the isotope number of the activation tag accordingly
    const unsigned int num_same_activations = count_num_same_activations(mol, 74u + new_bond_order, 2u);

    // Get the number of bonds in this molecule
    const unsigned int num_bonds = mol.getNumBonds();

    // Iterate over all bonds in the molecule
    for (unsigned int bond_idx = 0 ; bond_idx < num_bonds ; ++bond_idx)
    {  
        // Get the current bond corresponding to the index
        const RDKit::Bond *bond = mol.getBondWithIdx(bond_idx);

        // Get the atoms composing the bond
        const RDKit::Atom *bgn = bond->getBeginAtom();
        const RDKit::Atom *end = bond->getEndAtom();
        
        // Get the atoms' indices
        const unsigned int bgn_idx = bgn->getIdx();
        const unsigned int end_idx = end->getIdx();

        // The atoms need to have enough hydrogens
        if (hydrogen_map[bgn_idx] < 1 || hydrogen_map[end_idx] < 1)
        {
            continue;
        }

        // Has same bond order
        if (rev_bond_types[bond->getBondType()] != new_bond_order)
        {
            continue;
        }

        // Store this as a valid bond
        valid_bond_indices.emplace(bond_idx);
    }

    // Get the number of bonds that we can activate
    const unsigned int num_activations = static_cast<unsigned int>(valid_bond_indices.size());

    // Activate 1 to N bonds
    for (unsigned int n = 1; n <= num_activations; ++n)
    {
        // Generate all combinations to activate the input molecule N times
        const std::vector<std::vector<int>> combinations = generateCombinations(valid_bond_indices, n); 

        // Iterate over each combination of activations
        for (const std::vector<int>& activation_pattern : combinations) // 
        {
            
            // Assume we can activate
            bool can_activate = true;

            // Iterate over the indices of bond (i) that will be activated 
            for (int i = 0; i < n; ++i) 
            {    
                // Get the adjacent bond indices of the bond with index i 
                const std::unordered_set<unsigned int> adjacent_bonds = bond_adjacencies[activation_pattern[i]];

                // Iterate over the indices of bonds (j) that will be activated
                for (int j = i + 1; j < n; ++j) 
                {
                    // Check if the activations are adjacent
                    if (adjacent_bonds.count(activation_pattern[j]) > 0)
                    {
                        
                        // No adjacent double bonds
                        if (new_bond_order > 1)
                        {
                            can_activate = false;
                            break;
                        }
                        
                        // To store the atom index of eventual adjacent activated bonds
                        int shared_atom_idx = -1;

                        // Voodoo code I wrote after 4 cans of redbull
                        // For the bond (i), find the atoms composing it
                        std::pair<unsigned int,unsigned int> atom_pair = bond_information[i];
                        
                        // Atom indices composing the bond
                        const unsigned int r = atom_pair.first;
                        const unsigned int s = atom_pair.second;
                        
                        // For the adjacent bond (j), figure out the index of the atom they share
                        if (bond_information[j].first == r || bond_information[j].second == r)
                        {
                            shared_atom_idx = r;
                            break;
                        }
                        if (bond_information[j].first == s || bond_information[j].second == s)
                        {
                            shared_atom_idx = s;
                            break;
                        }

                        // If activated bonds do not share an atom
                        // Adjacent bonds that don't share, makes no sense
                        if (shared_atom_idx == -1)
                        {
                            // Go to next case
                            break;
                        }

                        // Activated bonds share an atom, AND
                        // there are not sufficient hydrogens to activate
                        // all bonds under this configuration
                        if (hydrogen_map[shared_atom_idx] < 2u) // Need two hydrogens!
                        {
                            // This activation pattern is not viable
                            can_activate = false;
                            break;
                        }
                    }
                }

                // If pattern is not viable
                if (!can_activate)
                {
                    // Move on to next pattern
                    break;
                }    
            }

            // If pattern is not viable
            if (!can_activate)
            {
                // Move on to next molecule
                continue;
            }            

            // Keep track of activated bonds
            std::unordered_map<unsigned int, unsigned int> seen_activate = {};

            // Create a work mol object we can transform with activations
            std::shared_ptr<RDKit::RWMol> activate_mol(new RDKit::RWMol(mol));

            // Isotope counter: we need this to be able to distinguish between
            // ring fusion patterns if two bonds have a bond in between that 
            // has two activation tags, but isn't supposed to be activated
            unsigned int isotope_counter = 0;

            // Iterate over all indices of bonds that will be activated
            for (const int& activate_idx : activation_pattern) 
            {
                // Get the activated bond
                const RDKit::Bond *activate_bond = activate_mol->getBondWithIdx(activate_idx);

                // And the atoms composing that bond
                RDKit::Atom *activate_bgn = activate_bond->getBeginAtom();
                RDKit::Atom *activate_end = activate_bond->getEndAtom();

                // Get those atoms' indices
                const unsigned int activate_bgn_idx = activate_bgn->getIdx();
                const unsigned int activate_end_idx = activate_end->getIdx();

                // To fix the hydrogen count later on
                if (seen_activate.find(activate_bgn_idx) == seen_activate.end())
                {
                    seen_activate.emplace(activate_bgn_idx, 1);
                }
                else
                {
                    ++seen_activate[activate_bgn_idx];
                }

                // To fix the hydrogen count later on
                if (seen_activate.find(activate_end_idx) == seen_activate.end())
                {
                    seen_activate.emplace(activate_end_idx, 1);
                }
                else
                {
                    ++seen_activate[activate_end_idx];
                }

                // Find the hydrogen bound to the bond's begin atom
                for (auto& neighbor : activate_mol->atomNeighbors(activate_bgn))
                {
                    // If the neighbor atom is a hydrogen
                    if (neighbor->getAtomicNum() == 1u)
                    {
                        // Turn it into an activation tag
                        neighbor->setAtomicNum(74u + new_bond_order);
                        
                        // Change the isotope depending on how many same activations have been observed before
                        neighbor->setIsotope(num_same_activations + isotope_counter);
                        break;
                    }
                }
                
                // Find the hydrogen bound the the bond's end atom
                for (auto& neighbor : activate_mol->atomNeighbors(activate_end))
                {
                    // If the neighbor atom is a hydrogen
                    if (neighbor->getAtomicNum() == 1u)
                    {
                        // Turn it into an activation tag
                        neighbor->setAtomicNum(74u + new_bond_order);
                        neighbor->setIsotope(num_same_activations + isotope_counter);
                        break;
                    }
                }

                // Start curating explicit hydrogens
                int num_bgn_hydrogens = std::max(0, static_cast<int>(hydrogen_map[activate_bgn_idx]) - \
                                                    static_cast<int>(seen_activate[activate_bgn_idx]));
                int num_end_hydrogens = std::max(0, static_cast<int>(hydrogen_map[activate_end_idx]) - \
                                                    static_cast<int>(seen_activate[activate_end_idx]));

                // Change the number of explicit hydrogens on the activated heavy atom
                activate_bgn->setNumExplicitHs(num_bgn_hydrogens);
                activate_bgn->setNoImplicit(false);
                activate_end->setNumExplicitHs(num_end_hydrogens);
                activate_end->setNoImplicit(false);

                // Increase the isotope counter
                ++isotope_counter;
            }

            // Remove explicit hydrogens again
            RDKit::MolOps::removeHs(*activate_mol);

            // Helps fix hydrogens
            RDKit::MolOps::sanitizeMol(*activate_mol);

            // Convert the molecule back to a regular RDKit mol obj
            RDKit::ROMol sanitized_mol(*activate_mol);

            // Convert the activated molecule to its corresponding SMILES string
            const std::string activate_smiles = RDKit::MolToSmiles(sanitized_mol);

            // Add the SMILES to the set of activated molecules
            unique_smiles.emplace(activate_smiles);
        }
    }
    
    //
    return unique_smiles;
}

// Introduce activation tags for spiro rings on all viable hydrogens
std::unordered_set<std::string> generate_spiro_activations(const RDKit::ROMol& mol, 
                                                           std::unordered_map<unsigned int,std::unordered_set<unsigned int>>& atom_adjacencies, 
                                                           std::unordered_map<unsigned int,unsigned int>& hydrogen_map, 
                                                           const bool verbose) 
{    
    // A set of unique SMILES strings
    std::unordered_set<std::string> unique_smiles = {};

    // A set containing atom indices of valid hydrogens
    std::unordered_set<unsigned int> valid_spiro_indices = {};

    // Count the number of times the current activation has been introduced
    // and adjust the isotope number of the activation tag accordingly
    const unsigned int num_same_activations = count_num_same_activations(mol, 74u, 1u);

    // Iterate over all heavy atoms in the molecule
    const unsigned int num_atoms = mol.getNumAtoms();
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx)
    {
        // Has to be a heavy atom
        if (mol.getAtomWithIdx(atom_idx)->getAtomicNum() == 1u)
        {
            continue;
        }

        // Heavy atom needs to have enough hydrogens
        if (hydrogen_map[atom_idx] < 2)
        {
            continue;
        }

        valid_spiro_indices.emplace(atom_idx);
    }

    // Get the number heavy atoms that we can activate
    const unsigned int num_activations = static_cast<unsigned int>(valid_spiro_indices.size());

    // Activate 1 to N atoms
    for (unsigned int n = 1; n <= num_activations; ++n)
    {
        
        // Generate all combinations to mutate input molecule N times
        const std::vector<std::vector<int>> combinations = generateCombinations(valid_spiro_indices, n); 

        // Iterate over each combination of activations
        for (const std::vector<int>& activation_pattern : combinations) // 
        {

            // Symmetry breaker for speedup
            //TODO

            // Create a work mol object we can transform with activations
            std::shared_ptr<RDKit::RWMol> activate_mol(new RDKit::RWMol(mol));

            // Isotope counter
            unsigned int isotope_counter = 0;
            
            // Iterate over all activations in this pattern
            for (const int& activate_idx : activation_pattern) 
            {
                // Get the heavy atom inside the transform molecule
                RDKit::Atom *activate_heavy = activate_mol->getAtomWithIdx(activate_idx);

                // Iterate over the hydrogens and activate them, don't have to change the bond_order
                int num_tagged = 0;

                // Iterate over all neighbors of the current atom 
                for (auto& activate_neighbor : activate_mol->atomNeighbors(activate_heavy))
                {
                    
                    // If the neighbor is a hydrogen
                    if (activate_neighbor->getAtomicNum() == 1u)
                    {
                        // Change the hydrogen into the activated tag
                        activate_neighbor->setAtomicNum(74u);
                        
                        // Change the isotope depending on how many same activations have been observed before
                        activate_neighbor->setIsotope(num_same_activations + isotope_counter);
                        ++num_tagged;

                        // If two hydrogens are activated for spiro
                        if (num_tagged > 1)
                        {
                            break;
                        }
                    }
                }

                // Increase the isotope counter
                ++isotope_counter;
            }
        
            // Remove explicit hydrogens again
            RDKit::MolOps::removeHs(*activate_mol);

            // Helps fix hydrogens
            RDKit::MolOps::sanitizeMol(*activate_mol);

            // Convert the molecule back to a regular RDKit mol obj
            RDKit::ROMol sanitized_mol(*activate_mol);

            // Convert the activated molecule to its corresponding SMILES string
            const std::string activate_smiles = RDKit::MolToSmiles(sanitized_mol);

            // Add the SMILES to the set of activated molecules
            unique_smiles.emplace(activate_smiles);
        }
    }
    return unique_smiles;
}

// Wrapper to generate all activations in a single function
void generate_activations(const RDKit::ROMol& mol, std::unordered_set<std::string>& all_smiles,
                          std::unordered_map<unsigned int,std::unordered_set<unsigned int>>& atom_adjacencies, 
                          std::unordered_map<unsigned int,std::unordered_set<unsigned int>>& bond_adjacencies,
                          std::unordered_map<unsigned int, std::pair<unsigned int,unsigned int>>& bond_information,
                          std::unordered_map<unsigned int,unsigned int>& hydrogen_map, bool verbose)
{
    // A set to store SMILES of activated scaffolds
    std::unordered_set<std::string> activated_smiles;

    // Generate all single coupling activations
    activated_smiles = generate_coupling_activations(mol, 1u, atom_adjacencies, hydrogen_map, verbose);

    all_smiles.insert(activated_smiles.begin(), activated_smiles.end());

    // Generate all double coupling activations
    activated_smiles = generate_coupling_activations(mol, 2u, atom_adjacencies, hydrogen_map, verbose);

    all_smiles.insert(activated_smiles.begin(), activated_smiles.end());

    // Generate all triple coupling activations
    activated_smiles = generate_coupling_activations(mol, 3u, atom_adjacencies, hydrogen_map, verbose);

    all_smiles.insert(activated_smiles.begin(), activated_smiles.end());
    
    // Generate all spiro cyclization activations
    activated_smiles = generate_spiro_activations(mol, atom_adjacencies, hydrogen_map, verbose);
    
    all_smiles.insert(activated_smiles.begin(), activated_smiles.end());

    // Generate all fusion activations
    activated_smiles = generate_fusion_activations(mol, 1u, bond_adjacencies, bond_information, hydrogen_map, verbose);

    all_smiles.insert(activated_smiles.begin(), activated_smiles.end());

    // Skip double fusions
    //activated_smiles = generate_fusion_activations(mol, 2u, bond_adjacencies, bond_information, hydrogen_map, verbose);

    //all_smiles.insert(activated_smiles.begin(), activated_smiles.end());
}

// Clear all isotope information from a molecule
std::string clear_isotope_info(const std::string &smiles)
{
    // Turn the line (SMILES) into a molecule
    std::shared_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(smiles));

    // Iterate over all atoms in the molecule
    const unsigned int num_atoms = mol->getNumAtoms();
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx)
    {
        //
        RDKit::Atom *atom = mol->getAtomWithIdx(atom_idx);
        
        // Remove the isotope information
        atom->setIsotope(0);
    }

    // Helps fix hydrogens
    RDKit::MolOps::sanitizeMol(*mol); 

    // Convert the molecule back to a regular RDKit mol obj
    RDKit::ROMol sanitized_mol(*mol);

    // Convert the curated molecule to its corresponding SMILES string
    const std::string curated_smiles = RDKit::MolToSmiles(sanitized_mol);

    return curated_smiles;
}

// Check whether the SMILES are valid in terms of activation patterns
bool validate_smiles(const std::string &smiles)
{
    // Strings representing activation tags
    const std::string spiro_str = "W";
    const std::string fusion_str = "Re";

    // Activation counter
    unsigned int num_tags = 0;

    // Used for shifting char lengths in SMILES string
    size_t spiro_pos = 0;
    size_t fusion_pos = 0;

    // Find all instances of "W" string
    while ((spiro_pos = smiles.find(spiro_str, spiro_pos)) != std::string::npos) 
    {
        // Increment counter
        ++num_tags;

        spiro_pos += spiro_str.length();
    }

    // Spiro activation needs two tags
    if (num_tags % 2 != 0)
    {
        return false;
    }
    
    // Find all instances of "Re" string
    while ((fusion_pos = smiles.find(fusion_str, fusion_pos)) != std::string::npos) 
    {
        // Increment counter
        ++num_tags;

        fusion_pos += fusion_str.length();
    }

    // Fusion activation needs two tags
    if (num_tags % 2 != 0)
    {
        return false;
    }

    // Is valid
    return true;
}

// Driver function
int main(int argc, char *argv[])
{
    
    // File names
    std::string in_file_path;
    std::string out_file_path;

    // Boolean to control verbose output
    bool verbose = false;

    // Argument parser
    int opt;
    while ((opt = getopt(argc, argv, "i:o:v:h")) != -1)
    {
        switch (opt)
        {
            case 'i':
                in_file_path = optarg;
                break;
            case 'o':
                out_file_path = optarg;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path>" << std::endl;
                return 0;
            default:
                std::cerr << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path>" << std::endl;
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
        std::cerr << "Error opening input file: " << in_file_path << std::endl;
        return 1;
    }

    // Verify you can open output file
    std::ofstream output_stream(out_file_path);
    if (!output_stream)
    {
        std::cerr << "Error opening output file: " << out_file_path << std::endl;
        return 1;
    }

    // String for reading lines
    std::string line;

    // Read the input file line by line
    while (std::getline(input_stream, line))
    {
        // Store SMILES of all activated scaffolds
        std::unordered_set<std::string> all_smiles;

        // Turn the line (SMILES) into a molecule
        std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(line));

        // This is the point where you set the hydrogens to explicit mode
        std::shared_ptr<RDKit::ROMol> protonated_mol(RDKit::MolOps::addHs(*mol));

        // Map containing hydrogen counts for each heavy atom
        std::unordered_map<unsigned int,unsigned int> hydrogen_map = get_hydrogen_map(*protonated_mol);

        // Map between atom indices and indices of adjacent atoms
        std::unordered_map<unsigned int,std::unordered_set<unsigned int>> atom_adjacencies = get_atom_adjacencies(*protonated_mol);

        // Map between bond indices and indices of adjacent bonds
        std::unordered_map<unsigned int,std::unordered_set<unsigned int>> bond_adjacencies = get_bond_adjacencies(*protonated_mol);

        // Map betweene bond indices and indices of atoms composing each bond
        std::unordered_map<unsigned int,std::pair<unsigned int,unsigned int>> bond_information = get_bond_information(*protonated_mol);

        // First you populate a set (all_smiles) with SMILES of scaffolds that have a single type of activation
        generate_activations(*protonated_mol, all_smiles, atom_adjacencies, bond_adjacencies, bond_information, hydrogen_map, verbose);

        // Copy elements from unordered_set to a vector
        std::vector<std::string> sortedElements(all_smiles.begin(), all_smiles.end());

        // Sort the vector
        std::sort(sortedElements.begin(), sortedElements.end());

        // For each activated molecule, reactivate in all configurations
        for (const std::string& smiles : sortedElements)
        { 
           
            // Turn the line (SMILES) into a molecule
            std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));

            // This is the point where you set the hydrogens to explicit mode
            std::shared_ptr<RDKit::ROMol> protonated_mol(RDKit::MolOps::addHs(*mol));

            // Map containing hydrogen counts for each heavy atom
            hydrogen_map = get_hydrogen_map(*protonated_mol);

            // Map between atom indices and indices of adjacent atoms
            atom_adjacencies = get_atom_adjacencies(*protonated_mol);

            // Map between bond indices and indices of adjacent bonds
            bond_adjacencies = get_bond_adjacencies(*protonated_mol);
            
            // Map betweene bond indices and indices of atoms composing each bond
            bond_information = get_bond_information(*protonated_mol);

            // Put a scaffold with multiple activations in the same set
            generate_activations(*protonated_mol, all_smiles, atom_adjacencies, bond_adjacencies, bond_information, hydrogen_map, verbose);
        }

        // A set to keep track of valid SMILES (with isotope info)
        std::unordered_set<std::string> clear_isotope_smiles;

        // Iterate over all activated molecules
        for (const std::string& smiles : all_smiles)
        { 
            std::string curated_smiles = clear_isotope_info(smiles);

            if (clear_isotope_smiles.find(curated_smiles) == clear_isotope_smiles.end())
            {
                clear_isotope_smiles.emplace(curated_smiles);

                if (validate_smiles(curated_smiles))
                {
                    output_stream << smiles + "\n"; 
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