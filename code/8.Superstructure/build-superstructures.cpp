/*
    Author: Andreas Luttens
    Contact: andreas.luttens@gmail.com
    Date: June 22, 2023
    Description: Attach substituents to activated scaffolds
*/

// Basic utilities
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>
#include <utility>
#include <iterator>
#include <filesystem>
#include <functional> 

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

// Remove intermediate file containing half-reacted molecules
void remove_intermediate_file(const std::string& remove_file_path)
{
    const std::filesystem::path file_path = remove_file_path;

    // Remove the file using std::filesystem
    try 
    {
        // Sometimes the intermediate file gets remove too quickly
        // and a sleep helps prevent that
        //if (std::filesystem::remove(file_path) != 0)
        //{
            // Wait for 3 seconds
            //std::this_thread::sleep_for(std::chrono::seconds(2));
        //}
        std::filesystem::remove(file_path);
    } 
    catch (const std::filesystem::filesystem_error& e) 
    {
        std::cerr << "Error removing file: " << e.what() << std::endl;
    }
}

// For a given scaffold molecule, and a specific reaction, react all substituent molecules in a file 
void connect_scaffold_with_substituent(RDKit::ChemicalReaction& rxn,const std::string& substituent_file_path,\
                                       const unsigned int& counter,std::vector<std::string> intermediate_file_names,\
                                       const unsigned int& compound_number, int& sym_number, bool& verbose)
{
    // Initialize the reaction
    rxn.initReactantMatchers();

    // Name of the file containing the scaffolds
    std::string scaffold_file_path = intermediate_file_names[counter];

    // Debugging purposes
    if (verbose)
    {
        std::cout << "we are reacting to:\t" << scaffold_file_path << std::endl;
    }

    // Open the file containing the scaffolds
    std::ifstream scaffold_stream(scaffold_file_path);

    // Open the file containing the substituents
    std::ifstream substituent_stream(substituent_file_path);

    // Format the name of the intermediate file
    std::ostringstream intermediate_file_path; 
    intermediate_file_path << std::to_string(compound_number)\
                           << "."\
                           << std::to_string(counter + 1)\
                           << +".intermediate.ism";
    
    // Debug purposes
    if (verbose)
    {
        std::cout << "Output reactionfile:\t" << intermediate_file_path.str() << std::endl;
    }

    // Open a file containing intermediate molecules
    std::ofstream intermediate_stream(intermediate_file_path.str());

    // To store unique molecules in case there are multiple equivalent reaction sites
    std::unordered_set<std::string> unique_molecules;

    // Create empty environment for reaction to take place 
    RDKit::MOL_SPTR_VECT reaction_flask;

    // Strings for reading lines
    std::string scaffold_line;
    std::string substituent_line;

    // Read the scaffold file line by line
    while (std::getline(scaffold_stream, scaffold_line))
    {
        // Debug purposes
        if (verbose)
        {
            std::cout << "current scaffold:\t" << scaffold_line << std::endl;
        }

        // Read the substituent file line by line
        while (std::getline(substituent_stream, substituent_line))
        {
            // Debug purposes
            if (verbose)
            {
                std::cout << "current substituent:\t" << substituent_line << std::endl;
            }
            
            // Turn the substituent smiles into an RDKit mol object
            RDKit::ROMol *substituent(RDKit::SmilesToMol(substituent_line));

            // Add hydrogens to the substituent
            RDKit::ROMol *protonated_substituent(RDKit::MolOps::addHs(*substituent));

            // Turn the scaffold smiles into an RDKit mol object
            // This has to happen for each substituent to connect them individually
            RDKit::ROMol *scaffold(RDKit::SmilesToMol(scaffold_line));

            // Add hydrogens to the scaffold
            RDKit::ROMol *protonated_scaffold(RDKit::MolOps::addHs(*scaffold));

            // Add the scaffold as the first reagent
            reaction_flask.emplace_back(RDKit::ROMOL_SPTR(protonated_scaffold));

            // Add the substituent as the second reagent
            reaction_flask.emplace_back(RDKit::ROMOL_SPTR(protonated_substituent));

            // Vector to store the converted molecules
            std::vector<RDKit::MOL_SPTR_VECT> connected_molecules;

            // Run the reaction at 1 or 2 reaction sites at a time
            // depending on symmetry factor, e.g. fusion two orientations
            connected_molecules = rxn.runReactants(reaction_flask, sym_number);

            // Iterate over all formed products
            for (const auto& product_vector : connected_molecules)
            {
                // RDKit has this as two nested loops 
                for (const auto& product : product_vector) 
                {  
                    // Obtain product SMILES
                    const std::string connected_smiles = RDKit::MolToSmiles(*product);

                    // Debug purposes
                    if (verbose)
                    {
                        std::cout << connected_smiles << std::endl;
                    }

                    // If there is more than 1 reaction site
                    if (sym_number > 1u)
                    {   
                        // Make sure the output molecules are unique
                        // Probably better ways to circumvent this potential duplication
                        if (unique_molecules.count(connected_smiles) > 0)
                        {
                            continue;
                        }
                        else
                        {
                            unique_molecules.emplace(connected_smiles);
                        }
                    }
                    
                    // Write the SMILES to the intermediate stream
                    intermediate_stream << connected_smiles + "\n";
                }
            }

            // Empty the flasks
            reaction_flask.clear();
        }

        // Go back to beginning of substituent file for the next activated scaffold
        substituent_stream.clear();
        substituent_stream.seekg(0, std::ios::beg);
    
    }
    // Close the file handles
    scaffold_stream.close();
    substituent_stream.close();
    intermediate_stream.close();
}

// Pull all molecules from input file through a chemical pattern filter
void pattern_filter(const unsigned int& counter, const std::string& input_file_path,\
                    const std::string& output_file_path, std::vector<RDKit::RWMol *> filters)
{
    // Open the input and output file streams
    std::ifstream input_stream(input_file_path);
    std::ofstream output_stream(output_file_path, std::ios::app); // Append mode

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

        // Bool to monitor molecule's validity
        bool is_valid = true;

        // Iterate over all patterns in the filter
        for (const auto& pattern : filters)
        {
            // Store pattern matches
            RDKit::MatchVectType hits;

            // If the mol hits a pattern
            if(RDKit::SubstructMatch(*mol, *pattern, hits)) 
            {
                // Is not valid molecule, break loop
                is_valid = false;
                break;
            }
        }

        // If the molecule survived all patterns
        if (is_valid)
        {
            // Deprotonate the molecule after pattern matching
            RDKit::ROMol *deprotonated_mol = RDKit::MolOps::removeHs(*mol);

            // Obtain product SMILES
            const std::string smiles = RDKit::MolToSmiles(*deprotonated_mol);

            // Hash the SMILES string
            size_t hashed_smiles = stringHasher(smiles);

            // If you have not seen this hash
            if (seen_hashes.find(hashed_smiles) == seen_hashes.end()) 
            {
                // Write to final output file
                output_stream << smiles + "\n";

                // Store the hashed SMILES string
                seen_hashes.insert(hashed_smiles);
            }
        }   
    }

    // Close the file streams
    input_stream.close();
    output_stream.close();
    
    // Remove the last intermediate file
    remove_intermediate_file(input_file_path);   
}

// Format a valid SMIRKS string to create a chemical reaction object later on
std::string format_smirks_string(std::string& input_string, bool& fuse_hack, bool& verbose)
{
    // All possible activation tags
    std::vector<std::string> activation_tags = {"Lu", "Hf", "Ta", "W", "Re", "Os"};

    // Container string to store current activation type
    std::string activation_tag;

    // Index corresponding to activation type
    int tag_idx = 0;
    
    // Iterate over all activation types and determine what this reaction is
    for (int t = 0; t < activation_tags.size(); ++t)
    {
        // Find which activation type it is
        size_t position = input_string.find(activation_tags[t]);

        // Should always be one of these six
        if (position != std::string::npos)
        {
            // Store the information
            activation_tag = activation_tags[t];
            tag_idx = t;
            break;
        }
    }

    // Variable to hold the formatted SMIRKS string
    std::ostringstream formatted_string;

    // Find the position of the substring to be replaced
    size_t bgn_position = input_string.find("[") + 1;
    size_t end_position = input_string.find("]");
    size_t length = end_position - bgn_position - activation_tag.size();
    std::string isotope_string = input_string.substr(bgn_position, length);

    // Container string for formatting the SMIRKS based on bond type
    std::string bond_string;

    std::string sym_char = "1";

    // Depending on what activation type
    switch (tag_idx)
    {
        case 0: // Lu : single bond
            bond_string = "-";

            // Form a valid SMIRKS string
            formatted_string << "[*:1]"\
                             << bond_string\
                             << "["\
                             << isotope_string\
                             << activation_tag\
                             << ":2].[*:3]"\
                             << bond_string\
                             << "["\
                             << activation_tag\
                             << ":4]>>[*:1]"\
                             << bond_string\
                             << "[*:3]";
            break;

        case 1: // Hf : double bond
            bond_string = "=";

            // Form a valid SMIRKS string
            formatted_string << "[*:1]"\
                             << bond_string\
                             << "["\
                             << isotope_string\
                             << activation_tag\
                             << ":2].[*:3]"\
                             << bond_string\
                             << "["\
                             << activation_tag\
                             << ":4]>>[*:1]"\
                             << bond_string\
                             << "[*:3]";
            break;

        case 2: // Ta : triple bond
            bond_string = "#";

            // Form a valid SMIRKS string
            formatted_string << "[*:1]"\
                             << bond_string\
                             << "["\
                             << isotope_string\
                             << activation_tag\
                             << ":2].[*:3]"\
                             << bond_string\
                             << "["\
                             << activation_tag\
                             << ":4]>>[*:1]"\
                             << bond_string\
                             << "[*:3]";
            break;

        case 3: // W : spiro
            bond_string = "-";

            // Form a valid SMIRKS string
            formatted_string << "[*:0][A:1](["\
                             << isotope_string\
                             << activation_tag\
                             << ":2])(["\
                             << isotope_string\
                             << activation_tag\
                             << ":3])[*:4].[A:5](["\
                             << activation_tag\
                             << ":6])["\
                             << activation_tag\
                             << ":7]>>[*:0][A:5][*:4]";
            break;    

        case 4: // Re : single fusion
            bond_string = "-";

            sym_char = "2"; // two ways to fuse two rings

            if (!fuse_hack)
            {
                // Form a valid SMIRKS string
                formatted_string << "[*:0][*:1](["\
                                << isotope_string\
                                << activation_tag\
                                << ":2])[*:3](["\
                                << isotope_string\
                                << activation_tag\
                                << ":4])[*:5].[*:6]@[*:7](["\
                                << activation_tag\
                                << ":8])[*:9](["\
                                << activation_tag\
                                << ":10])@[*:11]>>[*:0][*:1]([*:6])[*:3]([*:11])[*:5]"; 
            }
            else
            {
                // Form a valid SMIRKS string
                formatted_string << "[*:0][*:1](["\
                                << isotope_string\
                                << activation_tag\
                                << ":2])[*:3](["\
                                << isotope_string\
                                << activation_tag\
                                << ":4])[*:5].[*:6]1[*:7](["\
                                << activation_tag\
                                << ":8])[*:9]1(["\
                                << activation_tag\
                                << ":10])>>[*:0][*:1]1[*:6][*:3]1[*:5]";
            }
            break;

        case 5: // Os : double fusion
            bond_string = "";

            sym_char = "2"; // two ways to fuse two rings

            // Form a valid SMIRKS string
            formatted_string << "[*:0][*:1](["\
                             << isotope_string\
                             << activation_tag\
                             << ":2])=[*:3](["\
                             << isotope_string\
                             << activation_tag\
                             << ":4])[*:5].[*:6][*:7](["\
                             << activation_tag\
                             << ":8])[*:9](["\
                             << activation_tag\
                             << ":10])[*:11]>>[*:0][*:1](:[*:6]):[*:3](:[*:11])[*:5]";
            break;

        default:
            break;
    }
    
    if (verbose)
    {
        std::cout << formatted_string.str() << std::endl;
    }

    // Return as a clean c string + symmetry factor
    return formatted_string.str() + sym_char;
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

    // String to inform user
    std::string usage_string = " -i <input_file_path>\n\
                               -o <output_file_path>\n\
                               -f <filter_file_path>\n\
                               -v <verbose>";

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

    // Verify path
    if (filter_file_path.empty())
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

    // Create a vector with SMARTS pattern filters
    // from cheapest to most expensive calculation
    std::vector<RDKit::RWMol*> filters;

    // Read in the filter file, molecules need to match here
    read_filter_file(filter_file_path, filters);

    // String for reading lines
    std::string line;

    // Count which molecule we're finding superstructures for
    unsigned int compound_number = 0;

    // Read the input file line by line, each line is a molecule
    // with a recipe of instructions to build superstructures
    while (std::getline(input_stream, line))
    {
        // Debugging purposes
        if (verbose)
        {
            std::cout << compound_number << " + " << line << std::endl;
        }
        
        // Vector to store the tokens (substrings) of the current line
        std::vector<std::string> tokens;

        // Stream to tokenize the current line
        std::istringstream iss(line);

        // String to hold token substrings
        std::string token;

        // Tokenize the line
        while (std::getline(iss, token, ' ')) // delimiter is a white space
        {
            // Add the string to the vector
            tokens.emplace_back(token);
        }

        // Get the number of sites to connect onto this scaffold
        unsigned int num_tokens = tokens.size();

        // Debugging purposes
        if (verbose)
        {
            std::cout << "Num line tokens:\t" << num_tokens << std::endl;
        }
    
        // Vector to store names of files containing intermediate molecules
        std::vector<std::string> intermediate_file_names;
        
        // Intermediate file counter
        int intermediate_counter = 0;

        // Name of file containing intermediate molecules
        std::ostringstream intermediate_file_path;

        // Initialize a chemical reaction object
        RDKit::ChemicalReaction *rxn;

        // TODO
        int sym_number = 1;

        // Hack for three-mem ring fusions
        bool fuse_hack = false;

        // Iterate over all connection sites and corresponding reactions
        for (unsigned int token_idx = 1 ; token_idx < num_tokens ; ++token_idx)
        {
            // Current string
            std::string token_string = tokens[token_idx];
            
            // If the current token index is odd, its a activation type + location + isotope label
            if (token_idx%2 != 0)
            {
                // Debugging purposes
                if (verbose)
                {
                    std::cout << "activation type:\t" << token_string << std::endl;
                }

                // Hack for three-membered ring fusions
                if (tokens[token_idx + 1] == "1.single_fusion.ism")
                {
                    fuse_hack = true;
                }
                else
                {
                    fuse_hack = false;
                }
                
                // Format the reaction SMIRKS based activation type and isotope label
                std::string smirks_string = format_smirks_string(token_string, fuse_hack, verbose); 

                // The symmetry number will be the last char of the formatted string
                sym_number = std::atoi(&smirks_string.back());

                // Remove the last character
                smirks_string.resize(smirks_string.size() - 1);

                // Turn the reaction SMIRKS into a chemical reaction object
                rxn = RDKit::RxnSmartsToChemicalReaction(smirks_string);
            }

            // If the current token index is even, its a substituent library
            else
            {
                // Debugging purposes
                if (verbose)
                {
                    std::cout << "substituent library:\t" << token_string << std::endl;
                }

                // Clear the previous file name
                intermediate_file_path.str("");

                // Format the name of the file
                intermediate_file_path << std::to_string(compound_number)\
                                       << "."\
                                       << std::to_string(intermediate_counter)\
                                       << ".intermediate.ism";

                // Store the name of the file in vector
                intermediate_file_names.emplace_back(intermediate_file_path.str());

                // Write the initial blueprint scaffold to this intermediate file
                if (intermediate_counter == 0)
                {
                    std::ofstream intermediate_stream(intermediate_file_path.str());
                    intermediate_stream << tokens[0] + "\n";
                    intermediate_stream.close();
                }
                
                // Perform the last reaction with the current substituent library file
                connect_scaffold_with_substituent(*rxn, token_string, intermediate_counter, intermediate_file_names, compound_number, sym_number, verbose);

                // Fetch name of intermediate file name
                remove_intermediate_file(intermediate_file_path.str());

                // Sleep a few seconds ?
                
                // Increment the number of intermediate files
                ++intermediate_counter;
            }
        }

        // Clear the previous file name
        intermediate_file_path.str("");

        // Format the name of the file
        intermediate_file_path << std::to_string(compound_number)\
                               << "."\
                               << std::to_string(intermediate_counter)\
                               << ".intermediate.ism";

        // The last intermediate file will be pulled through a pattern filter
        // whose output is the true final output file
        pattern_filter(intermediate_counter, intermediate_file_path.str(), out_file_path, filters);

        // Increment the number of compounds processed
        ++compound_number;
    }

    // Close the input stream
    input_stream.close();

    // Signal success
    return 0;
}