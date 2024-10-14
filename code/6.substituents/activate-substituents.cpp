/*
    Author: Andreas Luttens
    Contact: andreas.luttens@gmail.com
    Date: June 22, 2023
    Description: Generate all activated analogues of input molecules
*/

// Basic utilities
#include <iostream>
#include <sstream>
#include <fstream>
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
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

// A map with bondtypes and bondorders
std::unordered_map<unsigned int, RDKit::Bond::BondType> bond_types = {{1u, RDKit::Bond::SINGLE},
                                                                      {2u, RDKit::Bond::DOUBLE},
                                                                      {3u, RDKit::Bond::TRIPLE}};

                                                                      // A map with bondtypes and bondorders
std::unordered_map<RDKit::Bond::BondType, unsigned int> rev_bond_types = {{RDKit::Bond::SINGLE, 1u},
                                                                          {RDKit::Bond::DOUBLE, 2u},
                                                                          {RDKit::Bond::TRIPLE, 3u}};

// Function to retrieve atom adjacencies
std::unordered_map<unsigned int, std::unordered_set<unsigned int>> get_atom_adjacencies(const RDKit::ROMol& mol){

    // initialize an empty dictionary to store the atom adjacencies
    std::unordered_map<unsigned int, std::unordered_set<unsigned int>> atom_adjacencies;

    // Iterate over all atoms in the molecule
    const unsigned int num_atoms = mol.getNumAtoms();
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx) {
        const RDKit::Atom *atom = mol.getAtomWithIdx(atom_idx);

        // initialize an array for the atoms index to store indices of adjacent atoms
        atom_adjacencies.emplace(atom_idx,std::unordered_set<unsigned int>{});

        // Iterate over all neighboring atoms
        for(const auto &neighbor_atom: mol.atomNeighbors(atom)) {
        
            const unsigned int neighbor_atom_idx = neighbor_atom->getIdx();

            // insert the atom's index in the array of indices
            atom_adjacencies[atom_idx].emplace(neighbor_atom_idx);
        }
    }
    return atom_adjacencies;
}

// Function to get information on bonds
std::unordered_map<unsigned int, std::pair<unsigned int,unsigned int>> get_bond_information(const RDKit::ROMol& mol)
{
    // Initialize an empty dictionary to store the bond information
    std::unordered_map<unsigned int, std::pair<unsigned int,unsigned int>> bond_information;

    // Iterate over all bonds in the molecule
    unsigned int num_bonds = mol.getNumBonds();
    for (unsigned int bond_idx = 0 ; bond_idx < num_bonds ; ++bond_idx)
    {
        const RDKit::Bond *bond = mol.getBondWithIdx(bond_idx);

        // Get the indices of the atoms forming this bond
        const unsigned int bond_bgn_idx = bond->getBeginAtomIdx();
        const unsigned int bond_end_idx = bond->getEndAtomIdx();

        // Create information pair
        std::pair<unsigned int,unsigned int> atom_pair = std::make_pair(bond_bgn_idx, bond_end_idx);

        // inserting key-value pairs
        bond_information.emplace(bond_idx, atom_pair);
    }

    return bond_information;
}

// Function to retrieve bond adjacencies
std::unordered_map<unsigned int, std::unordered_set<unsigned int>> get_bond_adjacencies(const RDKit::ROMol& mol)
{
    // initialize an empty dictionary to store the bond adjacencies
    std::unordered_map<unsigned int, std::unordered_set<unsigned int>> bond_adjacencies;

    // Iterate over all bonds in the molecule
    unsigned int num_bonds = mol.getNumBonds();
    for (unsigned int bond_idx = 0 ; bond_idx < num_bonds ; ++bond_idx)
    {
        const RDKit::Bond *bond = mol.getBondWithIdx(bond_idx);

        // Get the atoms forming the bond
        const RDKit::Atom *bond_bgn = bond->getBeginAtom();
        const RDKit::Atom *bond_end = bond->getEndAtom();

        // initialize an array for the bond index to store indices of adjacent bonds
        bond_adjacencies.emplace(bond_idx,std::unordered_set<unsigned int>{});

        // Iterate over all bonds through bond begin atom
        for(const auto &neighbor_bond: mol.atomBonds(bond_bgn))
        {
            unsigned int neighbor_bond_idx = neighbor_bond->getIdx();

            // insert the bond's index in the array of indices
            bond_adjacencies[bond_idx].emplace(neighbor_bond_idx);
        }

        // Iterate over all bonds through bond end atom
        for(const auto &neighbor_bond: mol.atomBonds(bond_end))
        {
            unsigned int neighbor_bond_idx = neighbor_bond->getIdx();

            // insert the bond's index in the array of indices
            bond_adjacencies[bond_idx].emplace(neighbor_bond_idx);
        }
    }

    return bond_adjacencies;
}

std::unordered_map<unsigned int,unsigned int> get_atom_symmetries(const RDKit::ROMol& mol)
{
    //
    std::unordered_map<unsigned int,unsigned int> symmetry_map = {};

    // Iterate over all atoms in the molecule
    const unsigned int num_atoms = mol.getNumAtoms();
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx)
    {
        //
        const RDKit::Atom *atom = mol.getAtomWithIdx(atom_idx);
        
        //
        unsigned int rank;
        atom->getProp("_CIPRank", rank);
        symmetry_map.emplace(atom->getIdx(), rank);
    }

    return symmetry_map;
}

//
std::unordered_map<unsigned int,std::vector<unsigned int>> get_equivalent_bonds(const RDKit::ROMol& mol, std::unordered_map<unsigned int,unsigned int>& atom_symmetries)
{
    std::map<std::pair<unsigned int,unsigned int>,std::vector<unsigned int>> bond_symmetries = {};

    // Iterate over all bonds in the molecule
    unsigned int num_bonds = mol.getNumBonds();
    for (unsigned int bond_idx = 0 ; bond_idx < num_bonds ; ++bond_idx)
    {
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
        for (const auto& sym_bond: sym_bonds)
        {

            equivalent_bonds.emplace(sym_bond, sym_bonds);
        }
    }

    //
    return equivalent_bonds;
}

//
std::unordered_map<unsigned int,std::vector<unsigned int>> get_equivalent_atoms(const RDKit::ROMol& mol, std::unordered_map<unsigned int,unsigned int>& atom_symmetries){

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
    
    // Equivalent atoms should be sorted!, so std::vector
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

std::unordered_set<unsigned int> get_unique_indices(std::unordered_map<unsigned int, std::vector<unsigned int>>& information)
{

    std::unordered_set<unsigned int> unique_indices = {};

    for (const auto& pair : information)
    {
        if (pair.second.size() == 1)
        {
            unique_indices.emplace(pair.first);
        }
    }
    return unique_indices;
}


std::unordered_set<unsigned int> get_unique_bonds(const RDKit::ROMol& mol, std::unordered_map<unsigned int,unsigned int>& atom_symmetries)
{

    std::unordered_set<unsigned int> unique_bonds = {};

    std::set<std::pair<unsigned int,unsigned int>> bond_symmetries = {};

    // Iterate over all bonds in the molecule
    unsigned int num_bonds = mol.getNumBonds();
    for (unsigned int bond_idx = 0 ; bond_idx < num_bonds ; ++bond_idx)
    {

        // Get the bond with the current bond index
        const RDKit::Bond *bond = mol.getBondWithIdx(bond_idx);

        // Get the atom indices of the atoms forming this bond
        const unsigned int bond_bgn_idx = bond->getBeginAtom()->getIdx();
        const unsigned int bond_end_idx = bond->getEndAtom()->getIdx();

        // Get the atom symmetries
        const unsigned int bond_bgn_sym = atom_symmetries[bond_bgn_idx];
        const unsigned int bond_end_sym = atom_symmetries[bond_end_idx];

        // Make a pair of the atom symmetries, which will become the bond symmetry
        std::pair<unsigned int,unsigned int> bond_symmetry = std::make_pair(bond_bgn_sym,bond_end_sym);

        // Make the reverse pair as well
        std::pair<unsigned int,unsigned int> rev_bond_symmetry = std::make_pair(bond_end_sym,bond_bgn_sym);

        // If this bond symmetry has not been seen before
        if (bond_symmetries.find(bond_symmetry) == bond_symmetries.end())
        {
            // Add the symmetry, and an empty set to store pairs later
            bond_symmetries.emplace(bond_symmetry);
            bond_symmetries.emplace(rev_bond_symmetry);

            // Add the bond index as a unique bond, bonds with same symmetries
            // will not be added in future loops, making this one unique
            unique_bonds.emplace(bond_idx);
        }
    }

    return unique_bonds;
}

std::unordered_set<unsigned int> get_unique_atoms(const RDKit::ROMol& mol, std::unordered_map<unsigned int,unsigned int>& atom_symmetries)
{
    std::unordered_set<unsigned int> unique_atoms = {};

    std::unordered_set<unsigned int> seen_symmetries = {};

    // Iterate over all atoms in the molecule
    unsigned int num_atoms = mol.getNumAtoms();
    for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx)
    {
        // Get the atom symmetries
        const unsigned int atom_sym = atom_symmetries[atom_idx];

        // If this bond symmetry has not been seen before
        if (seen_symmetries.find(atom_sym) == seen_symmetries.end())
        {
            // Add the symmetry, and an empty set to store pairs later
            seen_symmetries.emplace(atom_sym);

            // Add the atom index as a unique atom, atoms with same symmetries
            // will not be added in future loops, making this one unique
            unique_atoms.emplace(atom_idx);
        }
    }

    return unique_atoms;
}

std::unordered_map<unsigned int,unsigned int> get_hydrogen_map(RDKit::ROMol& mol)
{
    std::unordered_map<unsigned int,unsigned int> hydrogen_map = {};

    for (const auto& atom : mol.atoms())
    {
        unsigned int num_hydrogens = 0;
        for (const auto& neigbor : mol.atomNeighbors(atom))
        {
            if (neigbor->getAtomicNum() == 1u)
            {
                ++num_hydrogens;
            }
        }
        hydrogen_map.emplace(atom->getIdx(), num_hydrogens);
    }
    return hydrogen_map;
}

std::unordered_set<std::string> generate_coupling_activations(const RDKit::ROMol& mol, unsigned int new_bond_order, 
                                                              const std::unordered_set<unsigned int>& unique_bonds, 
                                                              std::unordered_map<unsigned int,unsigned int>& hydrogen_map, 
                                                              const bool verbose) 
{    
    // a set of unique SMILES strings
    std::unordered_set<std::string> unique_smiles = {};
    
    // Iterate over all bonds in the molecule
    for (const unsigned int& bond_idx : unique_bonds) 
    {
        // The bond must be unique
        if (unique_bonds.find(bond_idx) == unique_bonds.end())
        {
            continue;
        }

        const RDKit::Bond *bond = mol.getBondWithIdx(bond_idx);
        const RDKit::Atom *bgn = bond->getBeginAtom();
        const RDKit::Atom *end = bond->getEndAtom();

        const unsigned int bgn_idx = bgn->getIdx();
        const unsigned int end_idx = end->getIdx();

        unsigned int hydrogen_idx = 0;
        unsigned int heavy_idx = 0;

        if (bgn->getAtomicNum() == 1u)
        {
            hydrogen_idx = bgn_idx;
            heavy_idx = end_idx;
        }
        
        else if (end->getAtomicNum() == 1u)
        {
            hydrogen_idx = end_idx;
            heavy_idx = bgn_idx;
        }

        // The bond contains no hydrogen atom
        else
        {
            continue;
        }

        // Check if the heavy atom has enough hydrogens
        if (hydrogen_map[heavy_idx] < new_bond_order)
        {
            continue;
        }

        // Create a work mol object we can transform with activations
        std::shared_ptr<RDKit::RWMol> activate_mol(new RDKit::RWMol(mol));

        // Get the hydrogen atom inside the transform molecule
        RDKit::Atom *activate_hydrogen = activate_mol->getAtomWithIdx(hydrogen_idx);

        // Get the heavy atom that will carry the activation tag
        RDKit::Atom *activate_heavy = activate_mol->getAtomWithIdx(heavy_idx);

        // Get the bond between the hydrogen and heavy atoms
        RDKit::Bond *activate_bond = activate_mol->getBondBetweenAtoms(hydrogen_idx, heavy_idx);

        // Change the hydrogen into the activated tag
        activate_hydrogen->setAtomicNum(70u + new_bond_order);

        // Change the number of explicit hydrogens on the activated heavy atom
        activate_heavy->setNumExplicitHs(hydrogen_map[heavy_idx] - new_bond_order);

        activate_heavy->setNoImplicit(false);

        // Change the bond order between activated tag and heavy atom
        activate_bond->setBondType(bond_types[new_bond_order]);

        // Remove explicit hydrogens again
        RDKit::MolOps::removeHs(*activate_mol);

        // Help fix hydrogens
        RDKit::MolOps::sanitizeMol(*activate_mol); // does this fix the hydrogens?

        // Can I skip this?
        RDKit::ROMol sanitized_mol(*activate_mol);

        const std::string activate_smiles = RDKit::MolToSmiles(sanitized_mol);

        unique_smiles.emplace(activate_smiles);
    }

    return unique_smiles;
}

std::unordered_set<std::string> generate_fusion_activations(const RDKit::ROMol& mol, unsigned int new_bond_order, 
                                                            const std::unordered_set<unsigned int>& unique_bonds, 
                                                            std::unordered_map<unsigned int,unsigned int>& hydrogen_map, 
                                                            const bool verbose) 
{    
    // a set of unique SMILES strings
    std::unordered_set<std::string> unique_smiles = {};

    // Iterate over all bonds in the molecule
    for (const unsigned int& bond_idx : unique_bonds) 
    {
        // The bond must be unique to avoid repetition
        if (unique_bonds.find(bond_idx) == unique_bonds.end())
        {
            continue;
        }

        const RDKit::Bond *bond = mol.getBondWithIdx(bond_idx);
        const RDKit::Atom *bgn = bond->getBeginAtom();
        const RDKit::Atom *end = bond->getEndAtom();

        // if either one of the atoms forming the bond is a hydrogen, skip
        if (bgn->getAtomicNum() == 1u || end->getAtomicNum() == 1u)
        {
            continue;
        }

        // the bond_order has to match the query
        if (rev_bond_types[bond->getBondType()] != new_bond_order)
        {
            continue;
        }
        
        //
        const unsigned int bgn_idx = bgn->getIdx();
        const unsigned int end_idx = end->getIdx();

        const unsigned int num_bgn_hydrogens = hydrogen_map[bgn_idx];
        const unsigned int num_end_hydrogens = hydrogen_map[end_idx];

        // Both atoms in the bond need enough hydrogens to activate
        if (num_bgn_hydrogens < 1 || num_end_hydrogens < 1)
        {
            continue;
        }
        
        // Create a work mol object we can transform with activations
        std::shared_ptr<RDKit::RWMol> activate_mol(new RDKit::RWMol(mol));

        // Get the bgn atom inside the transform molecule
        RDKit::Atom *activate_bgn = activate_mol->getAtomWithIdx(bgn_idx);
        
        // Get the bgn atom inside the transform molecule
        RDKit::Atom *activate_end = activate_mol->getAtomWithIdx(end_idx);

        for (const auto& neighbor : activate_mol->atomNeighbors(activate_bgn))
        {
            if (neighbor->getAtomicNum() == 1u)
            {
                neighbor->setAtomicNum(74u + new_bond_order);
                break;
            }
        }
        
        for (const auto& neighbor : activate_mol->atomNeighbors(activate_end))
        {
            if (neighbor->getAtomicNum() == 1u)
            {
                neighbor->setAtomicNum(74u + new_bond_order);
                break;
            }
        }

        // Change the number of explicit hydrogens on the activated heavy atom
        activate_bgn->setNumExplicitHs(hydrogen_map[bgn_idx] - 1);
        activate_bgn->setNoImplicit(false);
        
        activate_end->setNumExplicitHs(hydrogen_map[end_idx] - 1);
        activate_end->setNoImplicit(false);
        
        // Remove explicit hydrogens again
        RDKit::MolOps::removeHs(*activate_mol);

        // Help fix hydrogens
        RDKit::MolOps::sanitizeMol(*activate_mol); // does this fix the hydrogens?

        // Can I skip this?
        RDKit::ROMol sanitized_mol(*activate_mol);

        const std::string activate_smiles = RDKit::MolToSmiles(sanitized_mol);

        // Insert the activated smiles
        unique_smiles.emplace(activate_smiles);
    }

    return unique_smiles;
}

// Function to build spiro activations of substituent molecule
std::unordered_set<std::string> generate_spiro_activations(const RDKit::ROMol& mol, 
                                                           const std::unordered_set<unsigned int>& unique_bonds, 
                                                           std::unordered_map<unsigned int,unsigned int>& hydrogen_map, 
                                                           const bool verbose) 
{    
    // a set of unique SMILES strings
    std::unordered_set<std::string> unique_smiles = {};

    // Iterate over all bonds in the molecule
    for (const unsigned int& bond_idx : unique_bonds) 
    {
    
        // The bond must be unique to avoid repetition
        if (unique_bonds.find(bond_idx) == unique_bonds.end())
        {
            continue;
        }

        // Get the current bond, and the atoms forming it
        const RDKit::Bond *bond = mol.getBondWithIdx(bond_idx);
        const RDKit::Atom *bgn = bond->getBeginAtom();
        const RDKit::Atom *end = bond->getEndAtom();

        // Get the atoms' indices
        const unsigned int bgn_idx = bgn->getIdx();
        const unsigned int end_idx = end->getIdx();

        // Figure out if the bond contains a hydrogen
        unsigned int heavy_idx = 0;
        if (bgn->getAtomicNum() == 1u)
        {
            heavy_idx = end_idx;
        }
        else if (end->getAtomicNum() == 1u)
        {
            heavy_idx = bgn_idx;
        }

        // The bond contains no hydrogen atom
        else
        {
            continue;
        }
        
        // We need two hydrogens on the heavy atom for spiro
        if (hydrogen_map[heavy_idx] < 2u)
        {
            continue;
        }
        
        // Create a work mol object we can transform with activations
        std::shared_ptr<RDKit::RWMol> activate_mol(new RDKit::RWMol(mol));
        
        // Get the heavy atom
        RDKit::Atom *activate_heavy = activate_mol->getAtomWithIdx(heavy_idx);

        unsigned int num_spiro_atoms = 0;

        // Iterate over the neighbors of the heavy atom
        for (const auto& neighbor : activate_mol->atomNeighbors(activate_heavy))
        {   
            // If it's a hydrogen
            if (neighbor->getAtomicNum() == 1u)
            {
                // Turn it into a tag
                neighbor->setAtomicNum(74u);
                ++num_spiro_atoms;

                // Increment counter, we need two tags
                if (num_spiro_atoms > 1)
                {
                    break;
                }
            }
        }

        // Change the number of explicit hydrogens on the activated heavy atom
        activate_heavy->setNumExplicitHs(hydrogen_map[heavy_idx] - 2);
        activate_heavy->setNoImplicit(false);
        
        // Remove explicit hydrogens again
        RDKit::MolOps::removeHs(*activate_mol);

        // Help fix hydrogens
        RDKit::MolOps::sanitizeMol(*activate_mol); // does this fix the hydrogens?

        // Can I skip this?
        RDKit::ROMol sanitized_mol(*activate_mol);

        const std::string activate_smiles = RDKit::MolToSmiles(sanitized_mol);

        // Insert these smiles, we cannot have unactivated smiles this time
        unique_smiles.emplace(activate_smiles);
    }

    return unique_smiles;
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

    // Verify you can open input file
    if (in_file_path.empty())
    {
        std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path>" << std::endl;
        return 1;
    }

    // Verify you can open output file
    if (out_file_path.empty())
    {
        std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path>" << std::endl;
        return 1;
    }

    // Open input file stream
    std::ifstream input_stream(in_file_path);
    if (!input_stream)
    {
        std::cerr << "Error opening input file: " << in_file_path << std::endl;
        return 1;
    }

    // Open output file stream
    std::ofstream output_stream(out_file_path);
    if (!output_stream)
    {
        std::cerr << "Error opening output file: " << out_file_path << std::endl;
        return 1;
    }

    // String for reading lines
    std::string line;

    // Counter to track number of completed compounds
    unsigned int n_activated = 0;

    // Ensure uniqueness of the output via string hashing
    std::unordered_set<size_t> seen_hashes;

    // Create a hashing object
    std::hash<std::string> stringHasher;

    // Read the input file line by line
    while (std::getline(input_stream, line))
    {
        // Turn the line (SMILES) into a molecule
        std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(line));

        // This is the point where you set the hydrogens to explicit mode
        std::shared_ptr<RDKit::ROMol> protonated_mol(RDKit::MolOps::addHs(*mol));

        // Perform hidden symmetry analysis by assigning labels to atoms
        // cleanIt=True,force=True,flagPossibleStereoCenters=True
        RDKit::MolOps::assignStereochemistry(*protonated_mol,true,true,true);

        // Map containing hydrogen counts for each heavy atom
        std::unordered_map<unsigned int,unsigned int> hydrogen_map = get_hydrogen_map(*protonated_mol);

        // Perform symmetry analysis
        std::unordered_map<unsigned int,unsigned int> atom_symmetries = get_atom_symmetries(*protonated_mol);
    
        // All unique bonds
        const std::unordered_set<unsigned int> unique_bonds = get_unique_bonds(*protonated_mol, atom_symmetries);

        // Store all SMILES
        std::unordered_set<std::string> all_smiles;

        std::unordered_set<std::string> activated_smiles;

        // Find all sigma bond activations for this molecule 
        activated_smiles = generate_coupling_activations(*protonated_mol, 1u, unique_bonds, hydrogen_map, verbose);

        all_smiles.insert(activated_smiles.begin(), activated_smiles.end());

        // Find all double bond activations for this molecule 
        activated_smiles = generate_coupling_activations(*protonated_mol, 2u, unique_bonds, hydrogen_map, verbose);

        all_smiles.insert(activated_smiles.begin(), activated_smiles.end());

        // Find all triple bond activations for this molecule 
        activated_smiles = generate_coupling_activations(*protonated_mol, 3u, unique_bonds, hydrogen_map, verbose);

        all_smiles.insert(activated_smiles.begin(), activated_smiles.end());

        // generate all single bond fusion activations for the current input molecule 
        activated_smiles = generate_fusion_activations(*protonated_mol, 1u, unique_bonds, hydrogen_map, verbose);

        all_smiles.insert(activated_smiles.begin(), activated_smiles.end());

        // generate all double bond fusion activations for the current input molecule 
        //activated_smiles = generate_fusion_activations(*protonated_mol, 2u, unique_bonds, hydrogen_map, verbose); skip for now

        //all_smiles.insert(activated_smiles.begin(), activated_smiles.end());

        // generate all spiro cyclization activations for the current input molecule 
        activated_smiles = generate_spiro_activations(*protonated_mol, unique_bonds, hydrogen_map, verbose);

        all_smiles.insert(activated_smiles.begin(), activated_smiles.end());

        // Iterate over all canonical tautomers of decorated molecules
        for (const std::string& smiles : all_smiles)
        { 
            // Hash the SMILES string
            size_t hashed_smiles = stringHasher(smiles);

            // If you have not seen this hash
            if (seen_hashes.find(hashed_smiles) == seen_hashes.end()) 
            {
                output_stream << smiles + "\n";
            }
        }

        // Verbose, counter
        ++n_activated;
        if (verbose && n_activated%1000 == 0)
        {
            std::cout << "Activated " << n_activated << " compounds!" << std::endl;
        }
    }

    // Close the file streams
    input_stream.close();
    output_stream.close();

    // Signal success 
    return 0;
}