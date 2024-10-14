/*
    Author: Andreas Luttens
    Contact: andreas.luttens@gmail.com
    Date: June 22, 2023
    Description: Generate all combinations you can activate a scaffold for superstructure generation
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>
#include <utility>
#include <unistd.h> 

//RDKit 
#include <GraphMol/GraphMol.h>
#include <GraphMol/Canon.h> // Check for symmetry later
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>


using namespace std;

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

    //
    RDKit::MolOps::sanitizeMol(*mol); 

    // Can I skip this?
    RDKit::ROMol sanitized_mol(*mol);

    // Transform the molecule into a valid SMILES string
    const std::string curated_smiles = RDKit::MolToSmiles(sanitized_mol);

    return curated_smiles;
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
    
    std::ifstream input_stream(in_file_path);
    if (!input_stream)
    {
        std::cerr << "Error opening input file: " << in_file_path << std::endl;
        return 1;
    }

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

    // Read the input file line by line
    while (std::getline(input_stream, line))
    {  
        std::string curated_smiles = clear_isotope_info(line);
        
        output_stream << curated_smiles + "\n";
    }

    // Close the file streams
    input_stream.close();
    output_stream.close();

    // Signal success 
    return 0;
}