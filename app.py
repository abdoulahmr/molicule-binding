from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, DataStructs
from rdkit.Chem.rdMolDescriptors import CalcTPSA
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

# Function to calculate Tanimoto similarity using Morgan Fingerprints
def calculate_similarity(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)

    # Calculate the Tanimoto similarity between the two fingerprints
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity

# Function to calculate molecular descriptors (LogP, Molecular Weight, etc.)
def calculate_molecular_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)

    descriptors = {
        'MolWt': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'NumHDonors': Lipinski.NumHDonors(mol),
        'NumHAcceptors': Lipinski.NumHAcceptors(mol),
        'TPSA': CalcTPSA(mol),  # Topological Polar Surface Area
        'HeavyAtomCount': Descriptors.HeavyAtomCount(mol),
    }

    return descriptors

# Function to compare functional groups using SMARTS patterns
def compare_functional_groups(smiles1, smiles2):
    patterns = {
        'Carboxylic Acid': '[CX3](=O)[O;H]',
        'Hydroxyl Group': '[OX2H]',
        'Aromatic Ring': 'a',
        'Amine': '[NX3;H2,H1;!$(NC=O)]',
    }

    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    functional_groups1 = {group: len(mol1.GetSubstructMatches(Chem.MolFromSmarts(pattern))) for group, pattern in patterns.items()}
    functional_groups2 = {group: len(mol2.GetSubstructMatches(Chem.MolFromSmarts(pattern))) for group, pattern in patterns.items()}

    # Count common functional groups
    common_groups = sum(min(functional_groups1[group], functional_groups2[group]) for group in patterns)
    return common_groups

# Main function to evaluate molecular binding potential
def evaluate_binding(smiles1, smiles2):
    # Generate 2D structures
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    # Draw and display molecules
    img1 = Draw.MolToImage(mol1, size=(200, 200))
    img2 = Draw.MolToImage(mol2, size=(200, 200))

    # Display the images
    display(img1)
    display(img2)

    # Step 1: Calculate similarity using Tanimoto
    similarity = calculate_similarity(smiles1, smiles2)
    print(f"Similarity (Tanimoto): {similarity:.3f}")

    # Step 2: Calculate molecular descriptors for both molecules
    descriptors1 = calculate_molecular_descriptors(smiles1)
    descriptors2 = calculate_molecular_descriptors(smiles2)

    print("\nMolecular Descriptors for Molecule X:")
    for key, value in descriptors1.items():
        print(f"{key}: {value}")

    print("\nMolecular Descriptors for Molecule Y:")
    for key, value in descriptors2.items():
        print(f"{key}: {value}")

    # Step 3: Compare functional groups
    functional_group_match = compare_functional_groups(smiles1, smiles2)
    print(f"\nCommon Functional Groups: {functional_group_match}")

    # Step 4: Calculate "binding affinity-like" score using descriptors
    # Note: This is a simplified scoring, based on molecular properties
    score1 = descriptors1['LogP'] - descriptors1['TPSA']  # Example heuristic
    score2 = descriptors2['LogP'] - descriptors2['TPSA']

    print(f"\nBinding Affinity-Like Score for Molecule X: {score1:.3f}")
    print(f"Binding Affinity-Like Score for Molecule Y: {score2:.3f}")

    # Step 5: Make prediction based on similarity and descriptors
    if similarity > 0.5 and functional_group_match > 0:
        if score1 > score2:  # Molecule X likely has a better binding profile
            print("\nMolecule X is likely to bind better to the receptor.")
        else:
            print("\nMolecule Y is likely to bind better to the receptor.")
    else:
        print("\nBoth molecules are too different in structure to bind similarly to the receptor.")


if __name__ == "__main__":
    # Define SMILES for two molecules
    smiles1 = str(input("Enter the first smiles and press enter: "))
    smiles2 = str(input("Enter the second smiles and press enter: "))

    # Run the evaluation
    try:
        evaluate_binding(smiles1, smiles2)
    except ModuleNotFoundError:
        print("run the first code to install the package")
    except:
        print("smiles not valid")