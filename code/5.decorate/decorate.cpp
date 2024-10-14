/*
    Author: Andreas Luttens
    Date: June 22, 2023
    Description: Introduce decorations in functionalized molecules
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
#include <filesystem>

// RDKit functions
#include <GraphMol/GraphMol.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <GraphMol/MolStandardize/Normalize.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

namespace fs = std::filesystem;

// Vector of strings representing chemical symbols we can have in decorated molecules
const std::vector<std::string> heavy_atoms = {"c", "C", "n", "N", "o", "O", "F", "s", "S", "Cl", "Br", "I"};

// Function that eliminates molecules that have more than two heavy halogens (Br, I)
void halogen_limiter(std::set<std::string>& molecules)
{
    // Strings representing bromine and iodine chemical symbols
    const std::string bromine_str = "Br";
    const std::string iodine_str = "I";

    // Iterate over all SMILES in the set of molecules
    for (auto smiles = molecules.begin(); smiles != molecules.end();)
    {
        // A halogen counter
        unsigned int num_halogens = 0;

        // Used for shifting char lengths in SMILES string
        size_t bromine_pos = 0;
        size_t iodine_pos = 0;

        // Find all instances of "I" string
        while ((iodine_pos = smiles->find(iodine_str, iodine_pos)) != std::string::npos) 
        {
            // Increment counter
            ++num_halogens;

            // If you found more than two, stop and eliminate
            if (num_halogens > 2)
            {
                break;
            }
            iodine_pos += iodine_str.length();
        }

        // Find all instances of "Br" string
        while ((bromine_pos = smiles->find(bromine_str, bromine_pos)) != std::string::npos) 
        {
            // Increment counter
            ++num_halogens;

            // If you found more than two, stop and eliminate
            if (num_halogens > 2)
            {
                break;
            }
            bromine_pos += bromine_str.length();
        }

        // Elimination
        if (num_halogens > 2)
        {
            smiles = molecules.erase(smiles);
        }

        // Go to the next SMILES string
        else
        {
            ++smiles;
        }
    }
}

// Function to eliminate molecules that have a nonaromatic-bond between heteroatoms
void remove_hetero_sigma(std::set<std::string>& molecules)
{
    // Atomic numbers of nitrogen and oxygen
    const std::unordered_set<int> heteros = {7,8};

    // Iterate only over newly generated molecules from the previous iteration
    for (auto smiles = molecules.begin(); smiles != molecules.end();) 
    {

        // Turn the line (SMILES) into a molecule
        RDKit::ROMol *mol(RDKit::SmilesToMol(*smiles));

        // Bool to store information
        bool is_hetero_sigma = false;

        // Find the number of bonds in this molecule
        const unsigned int num_bonds = mol->getNumBonds();

        // Iterate over all bonds in the molecule
        for (unsigned int bond_idx = 0 ; bond_idx < num_bonds ; ++bond_idx) 
        {
            // Retrieve the bond with this index
            const RDKit::Bond *bond = mol->getBondWithIdx(bond_idx);

            // If the bond is not aromatic TODO, what about oximes that are made later on?
            if (!bond->getIsAromatic()) 
            {
                // Check the atomic numbers of the bonded atoms
                const unsigned int bond_bgn_elem = bond->getBeginAtom()->getAtomicNum();
                const unsigned int bond_end_elem = bond->getEndAtom()->getAtomicNum();

                // Check if both atoms are hetero atoms
                if (heteros.find(bond_bgn_elem) != heteros.end() && heteros.find(bond_end_elem) != heteros.end())
                {
                    // Hetero-hetero
                    is_hetero_sigma = true;
                    break;
                }
            }   
        }

        // Eliminate if it has a hetero-hetero
        if (is_hetero_sigma)
        {
            smiles = molecules.erase(smiles);
        }

        // Go to the next SMILES string
        else{
            ++smiles;
        }
    }
}

// Function that returns canonicalized tautomers of given molecules
std::set<std::string> standardize_tautomers(std::set<std::string>& molecules)
{
    // Set to store standardized SMILES
    std::set<std::string> standardized_molecules;

    // Iterate over all given SMILES strings
    for (const auto& smiles : molecules)
    {
        // Turn it into a molecule
        std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
        
        // Initialize an enumerator of tautomers for this decorated molecule
        RDKit::MolStandardize::TautomerEnumerator enumerator;
        
        // Canonicalize the SMILES
        RDKit::ROMol *standard_mol = enumerator.canonicalize(*mol);
        
        // Store the canonicalized SMILES
        standardized_molecules.emplace(RDKit::MolToSmiles(*standard_mol));
    }

    return standardized_molecules;
}

// Read a file of reaction SMIRKS and corresponding names, store these pairs into a vector
void read_reaction_file(const std::string& reaction_file, std::vector<std::pair<std::string,RDKit::ChemicalReaction*>>& reactions)
{
    // Open the file containing the reaction SMIRKS patterns
    std::ifstream reaction_stream(reaction_file);
    if (!reaction_stream) 
    {
        std::cerr << "Error opening file: " << reaction_file << std::endl;
        return;
    }
    
    // String for the reaction names and SMIRKS patterns
    std::string reaction_name;
    std::string smirks_string;

    // Assuming columns are separated by tabs ('\t')
    while (reaction_stream >> reaction_name >> smirks_string) 
    {   
        // Create chemical reaction from reaction SMIRKS string
        RDKit::ChemicalReaction *rxn = RDKit::RxnSmartsToChemicalReaction(smirks_string.c_str());
        rxn->initReactantMatchers();
        
        // Store chemical reaction and its name
        reactions.emplace_back(std::make_pair(reaction_name, rxn));
    }

    // Close the stream
    reaction_stream.close();
}

// Function to insert all decoration onto a given input molecule and chemical reactions
std::set<std::string> generate_decorations(RDKit::ROMol& mol, 
                                           std::vector<std::pair<std::string,RDKit::ChemicalReaction*>>& reaction_library, 
                                           bool verbose)
{
    // A set of unique SMILES strings
    std::set<std::string> total_smiles = {RDKit::MolToSmiles(mol)};

    // Apply each reaction on a growing set of product molecule combinations
    for (const auto& reaction : reaction_library)
    {
        // Iterate over all molecules
        for (const std::string& smiles_to_convert : total_smiles) 
        {                
            // Create a work molecule that we can decorate by applying reactions
            RDKit::ROMol *mol_to_convert(RDKit::SmilesToMol(smiles_to_convert));

            // Store all SMILES
            std::set<std::string> pre_reaction_smiles = {};

            std::set<std::string> post_reaction_smiles = {};

            // Create empty environment for reaction to take place
            RDKit::MOL_SPTR_VECT flasks;

            // Add the compound as a reagent
            flasks.push_back(RDKit::ROMOL_SPTR(mol_to_convert));

            // Keep reacting the current molecule on all reaction sites, and make 
            // all possible products for the current reaction
            while (true)
            {
                // Retrieve the converted molecules
                std::vector<RDKit::MOL_SPTR_VECT> converted_molecules;

                // Run the reaction at 1 reaction site at a time
                converted_molecules = reaction.second->runReactants(flasks,1);

                // Empty the flasks
                flasks.clear();

                // Iterate over all formed products, we can kill one loop since we know the size, (1)
                for (const auto& product_vector : converted_molecules) // single loop if runReactant
                {
                    //
                    for (const auto& product : product_vector) 
                    {
                        // Obtain product SMILES
                        std::string converted_smiles = RDKit::MolToSmiles(*product);

                        // Store SMILES
                        post_reaction_smiles.emplace(converted_smiles);

                        // Add product back into flasks for further reactions
                        flasks.push_back(RDKit::ROMOL_SPTR(product));                        
                    }   
                }

                // If no new molecules were generated, break out of the while loop
                if (pre_reaction_smiles.size() == post_reaction_smiles.size())
                {
                    break;
                }

                //
                pre_reaction_smiles.insert(post_reaction_smiles.begin(), post_reaction_smiles.end());
            }

            // If the reaction was an iodination or bromination
            // We need to enfore a maximum two iodines or bromines in one molecule
            // These names are hardcoded right now but can be made modular ...
            if (reaction.first == "cOH->cBr" || reaction.first == "cOH->cI")
            {
                halogen_limiter(post_reaction_smiles);
            }

            // Clean up aromatic nitrogen introduction products
            // only for aza reactions, that way oxime formation is not affected
            // These names are hardcoded right now but can be made modular ...
            if (reaction.first == "aza1" || reaction.first == "aza2" || reaction.first == "aza3")
            {
                remove_hetero_sigma(post_reaction_smiles);
            }

            // Update the final smiles set
            total_smiles.insert(post_reaction_smiles.begin(), post_reaction_smiles.end());
        }
    }

    return total_smiles;
}

// Count heavy atoms without SMILES strings into mol objects
unsigned int count_heavy_atoms(const std::string& smiles)
{
    // A counter
    unsigned int num_heavy_atoms = 0;

    // Iterate over the heavy atom symbol strings
    for (const std::string& elem_str : heavy_atoms)
    {
        // Used to shift a specific number of characters based on symbol
        size_t elem_pos = 0;

        // If the string is found in the SMILES
        while ((elem_pos = smiles.find(elem_str, elem_pos)) != std::string::npos) 
        {
            // Increase the heavy atom count
            ++num_heavy_atoms;
            elem_pos += elem_str.length();
        }
    }
    
    return num_heavy_atoms;
}

// Driver function
int main(int argc, char *argv[])
{
    
    // File names
    std::string in_file_path;
    std::string out_file_path;
    std::string plus_out_file_path = "";
    std::string reaction_file_path;

    // Integer to compare heavy atom count of generated molecules
    int file_num_heavy_atoms = 0;

    // Boolean to control verbose output
    bool verbose = false;

    // Argument parser
    int opt;
    while ((opt = getopt(argc, argv, "i:o:r:n:v:h")) != -1)
    {
        switch (opt)
        {
            case 'i':
                in_file_path = optarg;
                break;
            case 'o':
                out_file_path = optarg;
                break;
            case 'r':
                reaction_file_path = optarg;
                break;
            case 'n':
                file_num_heavy_atoms = std::stoi(optarg);
                break;    
            case 'v':
                verbose = true;
                break;
            case 'h':
                std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -r <reaction_file_path> -n <num_heavy_atoms>" << std::endl;
                return 0;
            default:
                std::cerr << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -r <reaction_file_path> -n <num_heavy_atoms>" << std::endl;
                return 1;
        }
    }

    // Verify path
    if (in_file_path.empty())
    {
        std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -r <reaction_file_path> -n <num_heavy_atoms>" << std::endl;
        return 1;
    }

    // Verify path
    if (reaction_file_path.empty())
    {
        std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -r <reaction_file_path> -n <num_heavy_atoms>" << std::endl;
        return 1;
    }
    
    // 
    if (file_num_heavy_atoms == 0)
    {
        std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -r <reaction_file_path> -n <num_heavy_atoms>" << std::endl;
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

    // Extract the directory path from the file path of "foo"
    fs::path out_base_path = fs::path(out_file_path).parent_path();
    fs::path out_file_name = fs::path(out_file_path).filename();

    // Handles hierarchical folders
    if (out_base_path.string() == "")
    {
        // plus one heavy atom stream
        plus_out_file_path.append("plus1." + out_file_name.string());
    }
    else
    {
        plus_out_file_path.append(out_base_path.string() + "/plus1." + out_file_name.string());
    }

    // plus one heavy atom stream
    std::ofstream plus_output_stream(plus_out_file_path);
    if (!plus_output_stream)
    {
        std::cerr << "Error opening output file: " << plus_out_file_path << std::endl;
        return 1;
    }

    // Create a vector of reaction patterns
    std::vector<std::pair<std::string,RDKit::ChemicalReaction*>> reactions;

    // Read the reaction pattern file into the vector
    read_reaction_file(reaction_file_path, reactions);

    // Ensure uniqueness of the output via string hashing
    std::unordered_set<size_t> seen_hashes;

    // Create a hashing object
    std::hash<std::string> stringHasher;

    // String for reading lines
    std::string line;

    // Read the input file line by line
    while (std::getline(input_stream, line))
    {
        // Turn the line (SMILES) into a molecule
        std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(line));
        
        // Combinatorially functional groups for the current input molecule 
        std::set<std::string> decorated_smiles = generate_decorations(*mol, reactions, verbose);

        // Standardize the tautomers of the decorated molecules
        std::set<std::string> standardized_smiles = standardize_tautomers(decorated_smiles);

        // Iterate over all canonical tautomers of decorated molecules
        for (const std::string& smiles : standardized_smiles)
        {
            // Hash the SMILES string
            size_t hashed_smiles = stringHasher(smiles);

            // If you have not seen this hash
            if (seen_hashes.find(hashed_smiles) == seen_hashes.end()) 
            {
                // Count the number of heavy atoms
                unsigned int num_heavy_atoms = count_heavy_atoms(smiles);

                // Decide where to write the output SMILES to
                if (num_heavy_atoms == file_num_heavy_atoms)
                {
                    // Write the unsaturated SMILES string to output file
                    output_stream << smiles + "\n";

                    // Store the hashed SMILES string
                    seen_hashes.insert(hashed_smiles);
                }
                else
                {   
                    // Write the unsaturated SMILES string to output file
                    plus_output_stream << smiles + "\n";
                    
                    // Store the hashed SMILES string
                    seen_hashes.insert(hashed_smiles);
                }
            }
        }
    }

    // Close the file streams
    input_stream.close();
    output_stream.close();
    plus_output_stream.close();

    // Signal success 
    return 0;
}
