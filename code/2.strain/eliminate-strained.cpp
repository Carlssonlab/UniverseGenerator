/*
    Author: Andreas Luttens
    Contact: andreas.luttens@gmail.com
    Date: June 22, 2023
    Description: Filter molecules based on the smallest atomic volume their all-carbon tetrahedrons form
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
#include <unistd.h> 
#include <algorithm>
#include <cmath>

// RDKit functions
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>


// Helper function to calculate the magnitude of a 3D vector.
double magnitude(const std::vector<double>& vec)
{
    double sum_squares = 0.0;
    for (double value : vec) 
    {
        sum_squares += value * value;
    }
    return std::sqrt(sum_squares);
}

// Helper function to calculate the subtraction of two 3D vectors.
std::vector<double> vector_subtraction(const std::vector<double>& a, const std::vector<double>& b)
{
    std::vector<double> result(3);
    for (int i = 0; i < 3; ++i) 
    {
        result[i] = a[i] - b[i];
    }
    return result;
}

// Helper function to calculate the division of a 3D vector by its magnitude.
std::vector<double> unit_vector(const std::vector<double>& vec)
{
    double mag = magnitude(vec);
    std::vector<double> result(3);
    for (int i = 0; i < 3; ++i)
    {
        result[i] = vec[i] / mag;
    }
    return result;
}

// Function to create a unit vector based on two input vectors.
std::vector<double> update_coordinate(const std::vector<double>& ori, const std::vector<double>& dst) 
{
    std::vector<double> v = vector_subtraction(dst, ori);
    std::vector<double> v_hat = unit_vector(v);
    std::vector<double> new_coord(3);
    for (int i = 0; i < 3; ++i) 
    {
        new_coord[i] = ori[i] + v_hat[i];
    }
    return new_coord;
}

// Helper function to calculate the cross product of two 3D vectors.
std::vector<double> cross_product(const std::vector<double>& a, const std::vector<double>& b) 
{
    std::vector<double> result(3);
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    return result;
}

// Helper function to calculate the dot product of two 3D vectors.
double dot_product(const std::vector<double>& a, const std::vector<double>& b)
{
    double result = 0.0;
    for (int i = 0; i < 3; ++i) 
    {
        result += a[i] * b[i];
    }
    return result;
}

// Function to calculate the volume of a tetrahedron given three input vectors.
double calculate_tetrahedron(const std::vector<double>& v_a, const std::vector<double>& v_b, const std::vector<double>& v_c)
{
    std::vector<double> cross_ab = cross_product(v_a, v_b);
    double dot_abc = dot_product(cross_ab, v_c);
    return std::abs(dot_abc) / 6.0;
}

// Filter out molecule based on the smallest atomic value threshold
bool sav_filter(RDKit::ROMol& mol, double& threshold)
{
    // bool for tracking if the conformer has strained carbon atoms
    bool is_strained = false;
    
    // Get the molecule's conformer, it will be just one 'low' energy conformer
    RDKit::Conformer &conf = mol.getConformer();

    // Iterate over all tetravalent carbon atoms
    for(unsigned int i = 0; i < mol.getNumAtoms() ; ++i) 
    {
        // Get the atom with index i
        const RDKit::Atom *atom = mol.getAtomWithIdx(i);

        // Get the atomic number of this atom
        unsigned int elem_num = atom->getAtomicNum();

        // The atom has to be a carbon
        if (elem_num != 6)
        {
            continue;
        }

        // Create a map to store unit vectors
        std::map<int,std::vector<double>> connectivity = {};

        // Get the coordinate of this atom
        RDGeom::Point3D origin = conf.getAtomPos(i); 
        std::vector<double> ori_vec = {origin.x, origin.y, origin.z};
        
        // Count number of neighboring atoms
        int num_neighbors = 0;

        // Iterate over all neighboring atoms of this atom
        for(const int &neighbor_index: make_iterator_range(mol.getAtomNeighbors(atom))) 
        {    
            // Get the neighbor atom coordinate
            RDGeom::Point3D destination = conf.getAtomPos(neighbor_index);

            // Turn coordinate into vector
            std::vector<double> dest_vec = {destination.x, destination.y, destination.z};

            // Get unit vector coordinates
            std::vector<double> norm_vec = update_coordinate(ori_vec, dest_vec);

            // Update the neighboring atom counter
            ++num_neighbors;

            // Store the unit vector
            connectivity.emplace(num_neighbors, norm_vec);
        }

        // Check that there were at least four neighbors
        if (num_neighbors != 4)
        {
            continue;
        }
        
        else
        {
            // Define the vector of the tetrahedron
            std::vector<double> v_a = vector_subtraction(connectivity[2], connectivity[1]);
            std::vector<double> v_b = vector_subtraction(connectivity[3], connectivity[1]);
            std::vector<double> v_c = vector_subtraction(connectivity[4], connectivity[1]);
            
            // Calculate the smallest atomic volume
            double sav = std::abs(calculate_tetrahedron(v_a, v_b, v_c));

            // Compare smallest atomic volume with the threshold value
            if (sav < threshold)
            {
                is_strained = true;
            }
        }
    }

    return is_strained;
}

// Driver function
int main(int argc, char* argv[]) 
{
    // File names
    std::string in_file_path;
    std::string out_file_path;

    // Threshold for smallest atomic volume
    double threshold = 0.345;

    // Argument parser
    int opt;
    while ((opt = getopt(argc, argv, "i:o:t:h")) != -1) 
    {
        switch (opt) 
        {
            case 'i':
                in_file_path = optarg;
                break;
            case 'o':
                out_file_path = optarg;
                break;
            case 't':
                threshold = atof(optarg);
                break;
            case 'h':
                std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -t <threshold>" << std::endl;
                return 0;
            default:
                std::cerr << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -t <threshold>" << std::endl;
                return 1;
        }
    }

    // Verify you can open input file
    if (in_file_path.empty()) 
    {
        std::cerr << "Error: Input file path not provided. Use the -i option." << std::endl;
        std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -t <threshold>" << std::endl;
        return 1;
    }

    // Verify you can open output file
    if (out_file_path.empty()) 
    {
        std::cerr << "Error: Output file path not provided. Use the -o option." << std::endl;
        std::cout << "Usage: " << argv[0] << " -i <input_file_path> -o <output_file_path> -t <threshold>" << std::endl;
        return 1;
    }

    std::ifstream input_file(in_file_path);
    if (!input_file) 
    {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    std::ofstream output_file(out_file_path);
    if (!output_file) 
    {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    std::string line;
    std::string current_mol2_block = "";

    //RDKit::Mol2ParserParams params;

    // Read the input file line by line
    while (std::getline(input_file, line)) 
    {
        if (line.find("@<TRIPOS>MOLECULE") == 0) 
        {
            // If the line starts with the provided start string,
            // push the current_mol2_block into the vector and start a new block.
            if (!current_mol2_block.empty()) // and current_mol2_block.find("@<TRIPOS>MOLECULE") == 0 
            {

                // Create RDKit mol object from mol2block text
                std::shared_ptr<RDKit::ROMol> mol(RDKit::Mol2BlockToMol(current_mol2_block, true, false));

                // If molecule passes SAV filter
                if (!sav_filter(*mol, threshold))
                {
                    // Create a work mol object we can transform with activations
                    std::shared_ptr<RDKit::RWMol> protonated_mol(new RDKit::RWMol(*mol));

                    // Remove explicit hydrogens again
                    RDKit::MolOps::removeHs(*protonated_mol);

                    // Write the SMILES to the output file
                    output_file << RDKit::MolToSmiles(*protonated_mol) + "\n";
                }
                
                // Empty block of text
                current_mol2_block.clear();
            }
        }
        current_mol2_block += line + "\n";
    }

    // Analyse the last block of text
    if (!current_mol2_block.empty()) {
        std::shared_ptr<RDKit::ROMol> mol(RDKit::Mol2BlockToMol(current_mol2_block, true, false));

        // If molecule passes SAV filter
        if (!sav_filter(*mol, threshold))
        {
            // Create a work mol object we can transform with activations
            std::shared_ptr<RDKit::RWMol> protonated_mol(new RDKit::RWMol(*mol));

            // Remove explicit hydrogens again
            RDKit::MolOps::removeHs(*protonated_mol);

            // Write the SMILES to the output file
            output_file << RDKit::MolToSmiles(*protonated_mol) + "\n";
        }
    }

    // Close file streams
    input_file.close();
    output_file.close();

    // Signal success
    return 0;

}
