/*
    Author: Andreas Luttens
    Contact: andreas.luttens@gmail.com
    Date: June 22, 2023
    Description: Convert adjacency matrices into hydrocarbon SMILES
*/

// Basic utilities
#include <iostream>
#include <vector>
#include <string>

// RDKit functions
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

// Function to convert adjacency matrix to a SMILES
std::string adjmat2smiles(const std::vector<std::vector<unsigned int>> &adjmat) 
{
    // Get the number of rows (or columns) in the adjancency matrix
    int num_atoms = adjmat.size();

    // Create a new RDKit work molecule that we'll give atoms and bonds
    // Hydrocarbon based on adjacency matrix, all single bonds, fill in
    // hydrogens based on valencies of carbons

    // Create a string by repeating the substring N times
    std::string disconnect_smiles;
    for (int i = 0; i < num_atoms; ++i) 
    {
        disconnect_smiles += "C.";
    }

    // Remove the last redundant dot
    disconnect_smiles.resize(disconnect_smiles.size() - 1);

    // Convert that SMILES into a molecule of N carbons
    RDKit::ROMol *molecule = RDKit::SmilesToMol(disconnect_smiles);

    // Create a work molecule that we can add new bonds into
    std::unique_ptr<RDKit::RWMol> work_mol(new RDKit::RWMol(*molecule));

    // Scan upper triangular matrix
    for (int a = 0; a < num_atoms - 1; ++a) 
    {
        for (int b = a + 1; b < num_atoms; ++b) 
        {
            // If adjacency matrix element is not zero
            if (adjmat[a][b] == 1)
            {
                // Make an RDKit object of this atom
                RDKit::Atom *atom_a = work_mol->getAtomWithIdx(a);

                // Make an RDKit object of this atom
                RDKit::Atom *atom_b = work_mol->getAtomWithIdx(b);
                
                // Add a single sigma bond in the molecule
                work_mol->addBond(atom_a,atom_b,RDKit::Bond::SINGLE);
            }
        }
    }

    // Helps fix hydrogens
    RDKit::MolOps::sanitizeMol(*work_mol); 

    // Convert the molecule to its corresponding SMILES string
    const std::string smiles = RDKit::MolToSmiles(RDKit::ROMol(*work_mol));

    return smiles;
}

// Driver function
int main(int argc, char* argv[]) 
{    
    // File names
    std::string out_file_path;

    // Argument parser
    int opt;
    while ((opt = getopt(argc, argv, "o:h")) != -1) 
    {
        switch (opt) 
        {
            case 'o':
                out_file_path = optarg;
                break;
            case 'h':
                std::cout << "Usage: " << argv[0] << "-o <output_file_path>" << std::endl;
                return 0;
            default:
                std::cerr << "Usage: " << argv[0] << "-o <output_file_path>" << std::endl;
                return 1;
        }
    }

    // Verify path
    if (out_file_path.empty()) 
    {
        std::cerr << "Error: Output file path not provided. Use the -o option." << std::endl;
        std::cout << "Usage: " << argv[0] << "-o <output_file_path>" << std::endl;
        return 1;
    }

    // Verify you can open output file
    std::ofstream output_file(out_file_path);
    if (!output_file) 
    {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    // String to store a line from the input file
    std::string line;

    // Placeholder for an adjancency matrix
    std::vector<std::vector<unsigned int>> adjmat;

    //
    bool grow = false;

    // Read the input file line by line
    while (std::getline(std::cin, line)) 
    {
        // If the line isn't empty
        if (!line.empty())
        {
            // If the line starts with a G and the matrix has data
            if (line.rfind("G",0) == 0)
            {
                grow = false;

                if (adjmat.size() != 0)
                {
                    // Convert adjacency matrix to valid SMILES of hydrocarbon
                    std::string smiles = adjmat2smiles(adjmat);

                    // Write the SMILES to the output file
                    output_file << smiles + "\n";
                }

                // Placeholder for an adjancency matrix
                adjmat.clear();
            }

            // Toggle to grow the matrix with bond orders
            else
            {
                grow = true;
            }

            if (grow)
            {
                // Initialize a new row
                std::vector<unsigned int> row;

                // Iterate over the bond orders in the line
                int size = line.size();
                for (int c = 0; c < size; ++c)
                {
                    // Add the bond order to the vector
                    row.emplace_back(line.at(c) - '0');
                }

                // Add the row to the matrix
                adjmat.emplace_back(row);
            }
        }
    }

    // Close the file stream
    output_file.close();

    // Signal success
    return 0;
}