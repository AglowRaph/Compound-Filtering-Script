# Compound-Filtering-Script

# Example Usage
if __name__ == "__main__":
    # Input: List of SMILES strings
    smiles_list = ["CCO", "CC(=O)O", "CC(C)O"]  # Replace with your dataset

    # Step 1: Lipinski Filtering
    lipinski_results = lipinski_filter(smiles_list)
    print("Lipinski's Rule of Five Filter Results:")
    print(lipinski_results)

    # Step 2: Advanced Filtering
    advanced_results = advanced_filter(smiles_list)
    print("Advanced Filtering Results:")
    print(advanced_results)
