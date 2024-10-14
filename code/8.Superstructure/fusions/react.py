from rdkit import Chem
from rdkit.Chem import AllChem

# Define the SMIRKS pattern for the reaction (adjust this to your specific reaction)
smirks_pattern = "[*:0][*:1]([0Re:2])[*:3]([0Re:4])[*:5].[*:6]@[*:7]([Re:8])[*:9]([Re:10])@[*:6]>>[*:0]1[*:1]([*:6])[*:3]1[*:5]"

# Define the SMILES strings for the two molecules
smiles1 = "[Re]C1C([Re])CCC1"  # Cyclopentane
smiles2 = "[Re]C1NC1[Re]"    # Aziridine

# Convert SMILES to RDKit molecules
mol1 = Chem.MolFromSmiles(smiles1)
mol2 = Chem.MolFromSmiles(smiles2)

# Create the reaction object from the SMIRKS pattern
reaction = AllChem.ReactionFromSmarts(smirks_pattern)

# Perform the reaction (using the reactants as input)
reactants = [mol1, mol2]
products = reaction.RunReactants(reactants)

# Display the products
for product_set in products:
    for product in product_set:
        print(Chem.MolToSmiles(product))

