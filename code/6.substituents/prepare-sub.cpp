/*
    Author: Andreas Luttens
    Date: June 22, 2023
    Description: 
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
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

using namespace std;

bool is_activated_bond(const RDKit::ROMol &mol, const RDKit::Bond &bond)
{
    
    std::unordered_set<unsigned int> activation_tags = {74u,75u,76u,85u};
        
    // Get the atoms
    const RDKit::Atom *bgn = bond.getBeginAtom();
    
    //
    const RDKit::Atom *end = bond.getEndAtom();

    const unsigned int bgn_elem = bgn->getAtomicNum();
    const unsigned int end_elem = end->getAtomicNum();

    if ((bgn_elem != 6u) || (end_elem != 6u))
    {
        return false;
    }

    const unsigned int bgn_hcount = bgn->getNumImplicitHs();
    const unsigned int end_hcount = end->getNumImplicitHs();

    if ((bgn_hcount == 0) || (end_hcount == 0))
    {
        return false;
    }

    //
    std::vector<unsigned int> bond_tags;

    // Iterate over all neighboring atoms
    for(const auto &neighbor_atom: mol.atomNeighbors(bgn))
    {
        const unsigned int neighbor_elem_num = neighbor_atom->getAtomicNum();

        if (activation_tags.count(neighbor_elem_num) > 0)
        {
            bond_tags.emplace_back(neighbor_elem_num);
        }
    }

    // Iterate over all neighboring atoms
    for(const auto &neighbor_atom: mol.atomNeighbors(end))
    {
        const unsigned int neighbor_elem_num = neighbor_atom->getAtomicNum();

        if (activation_tags.count(neighbor_elem_num) > 0)
        {
            bond_tags.emplace_back(neighbor_elem_num);
        }
    }

    for (const unsigned int& tag : bond_tags)
    {
        if (count(bond_tags.begin(), bond_tags.end(), tag) > 1)
        {
            return true;
        }   
    }

    return false;
}

bool is_activated_atom(const RDKit::ROMol &mol, const RDKit::Atom &atom)
{
    
    std::unordered_set<unsigned int> activation_tags = {74u,75u,76u,85u};
        
    const unsigned int elem = atom.getAtomicNum();

    if (elem != 6u)
    {
        return false;
    }

    //const unsigned int hcount = atom.getNumImplicitHs();

    //if (hcount == 0)
    //{
    //    return false;
    //}

    //
    std::vector<unsigned int> bond_tags;

    // Iterate over all neighboring atoms
    for(const auto &neighbor_atom: mol.atomNeighbors(&atom))
    {
        const unsigned int neighbor_elem_num = neighbor_atom->getAtomicNum();

        if (activation_tags.count(neighbor_elem_num) > 0)
        {
            bond_tags.emplace_back(neighbor_elem_num);
        }
    }

    for (const unsigned int& tag : bond_tags)
    {
        if (count(bond_tags.begin(), bond_tags.end(), tag) > 1)
        {
            return true;
        }   
    }

    return false;
}

// Driver function
int main(int argc, char *argv[])
{
    
    // File names
    std::string in_file_path;
    std::string out_file_path;

    // Real isotope name
    std::string old_label;

    // Boolean to control verbose output
    bool verbose = false;

    bool aromatic = false;

    const std::string usage_string = " -i <input_file_path> -o <output_file_path> -l <label> -a <aromatic>";

    // Argument parser
    int opt;
    while ((opt = getopt(argc, argv, "i:o:l:a:v:h")) != -1)
    {
        switch (opt)
        {
            case 'i':
                in_file_path = optarg;
                break;
            case 'o':
                out_file_path = optarg;
                break;
            case 'l':
                old_label = optarg;
                break;
            case 'a':
                aromatic = true;
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

    // Verify you can open input file
    if (in_file_path.empty())
    {
        std::cout << "Usage: " << argv[0] << usage_string << std::endl;
        return 1;
    }

    // Verify you can open output file
    if (out_file_path.empty())
    {
        std::cout << "Usage: " << argv[0] << usage_string << std::endl;
        return 1;
    }

    std::ifstream input_stream(in_file_path);
    if (!input_stream)
    {
        std::cout << "Error opening input file: " << in_file_path << std::endl;
        return 1;
    }

    std::ofstream output_stream(out_file_path);
    if (!output_stream)
    {
        std::cout << "Error opening output file: " << out_file_path << std::endl;
        return 1;
    }

    // String for reading lines
    std::string line;

    // Read the input file line by line
    while (std::getline(input_stream, line))
    {
        // Turn the line (SMILES) into a molecule
        std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(line));

        // Perform ring analysis
        RDKit::MolOps::findSSSR(*mol);

        // Store ring info, can I do a const here? This molecule won't change
        RDKit::RingInfo* ring_info = mol->getRingInfo();

        // Iterate over all bonds in the molecule
        const unsigned int num_bonds = mol->getNumBonds();
        for (unsigned int bond_idx = 0 ; bond_idx < num_bonds ; ++bond_idx) 
        {
            // Retrieve the molecule's bond
            const RDKit::Bond *bond = mol->getBondWithIdx(bond_idx);

            const unsigned int type = bond->getBondType();

            if (is_activated_bond(*mol, *bond) && ring_info->numBondRings(bond_idx)); // bond->getBondType() == RDKit::Bond::AROMATIC && 
            {
                output_stream << line + "\n";
            }
        }        

        /*
        // Iterate over all atoms in the molecule
        const unsigned int num_atoms = mol->getNumAtoms();
        for (unsigned int atom_idx = 0 ; atom_idx < num_atoms ; ++atom_idx) 
        {
            // Retrieve the molecule's atom
            const RDKit::Atom *atom = mol->getAtomWithIdx(atom_idx);

            if (is_activated_atom(*mol, *atom) && ring_info->numAtomRings(atom_idx))
            {
                output_stream << line + "\n";
            }
        }
        */
        
    }

    input_stream.close();
    output_stream.close();

    // signal success
    return 0;
}
