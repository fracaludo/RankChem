from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from morfeus import read_xyz, XTB

def mol_to_xyz(smiles_list):
    """
    Convert a list of SMILES strings to their XYZ format.

    Parameters
    ----------
    smiles_list : list of str
        A list of text strings where each string is a SMILES (Simplified
        Molecular Input Line Entry System) notation representing a molecule.

    Returns
    -------
    list of str
        A list of XYZ formatted strings. Each string represents the 3D coordinates
        of a molecule corresponding to the input SMILES strings.

    Examples
    --------
    >>> mol_to_xyz(["C(=O)O", "CCO"])
    ['3\nC(=O)O\nC 0.0000 0.0000 0.0000\nO 1.2093 0.0000 0.0000\nO -1.2093 0.0000 0.0000', 
     '6\nCCO\nC 0.0000 0.0000 0.0000\nC 1.2093 0.0000 0.0000\nO 2.4186 0.0000 0.0000']

    """
    results = [] 
    
    # Iterate over each SMILES string in the list
    for smiles in smiles_list:
        # SMILES to RDKit
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, rdDistGeom.ETKDG()) # 3D coordinates 
        num_atoms = mol.GetNumAtoms() # RDKit XYZ string
        xyz_lines = [str(num_atoms)]
        xyz_lines.append(smiles)
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            conf = mol.GetConformer()
            pos = conf.GetAtomPosition(atom.GetIdx())
            xyz_lines.append(f"{symbol} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}")
        
        # Append the XYZ to results list
        results.append("\n".join(xyz_lines))

    return results


def get_molecules(smiles_list, N):
        """
    Compute Fukui indices for electrophilicity and nucleophilicity for a list of molecules.

    Parameters
    ----------
    smiles_list : list of str
        A list of text strings where each string is a SMILES (Simplified Molecular Input Line Entry System) notation representing a molecule.
    N : int
        The number of times to calculate Fukui indices for each molecule to obtain an average.

    Returns
    -------
    None
        This function prints the average maximum values of Fukui indices for electrophilicity or nucleophilicity for each molecule in the input SMILES list.

    Examples
    --------
    >>> get_molecules(["C(=O)O", "CCO"], 10)
    {'C(=O)O': 0.1234, 'CCO': 0.5678}
    
    Notes
    -----
    This function relies on external functions and modules such as `mol_to_xyz`, `read_xyz`, and the `XTB` class to perform its operations. 
    """
    xyz_contents = mol_to_xyz(smiles_list)# SMILES list to XYZ 
    n_nucs_max = []
    n_eles_max = []
    for i, xyz_content in enumerate(xyz_contents):
        fukui_electrophilicity = []
        fukui_nucleophilicity = []
        for _ in range(N):
            xyz_file = f"molecule_{i}.xyz"  # XYZ to a file
            with open(xyz_file, 'w') as file:
                file.write(xyz_content)
            
            elements, coordinates = read_xyz(xyz_file) # Morfeus to read XYZ file
            
            xtb = XTB(elements, coordinates) # Create  XTB object 
            
            fukui_electrophilicity.append(xtb.get_fukui("electrophilicity"))
            fukui_nucleophilicity.append(xtb.get_fukui("nucleophilicity"))

        # Fukui electrophilicity max value
        max_values_electrophilicity = [max(f.values()) for f in fukui_electrophilicity]
        n_eles_max.append(sum(max_values_electrophilicity) / len(max_values_electrophilicity))

        # Fukui nucleophilicity max value
        max_values_nucleophilicity = [max(f.values()) for f in fukui_nucleophilicity]
        n_nucs_max.append(sum(max_values_nucleophilicity) / len(max_values_nucleophilicity))

    n_eles_max = dict(zip(smiles_list, n_eles_max))
    n_nucs_max = dict(zip(smiles_list, n_nucs_max))

    if nucl_elec == "E":
        print(n_nucs_max) #nucleophilicity max values
    else:
        print(n_eles_max) #electrophilicity max values

user_input = input("Enter a list of SMILES strings separated by commas: ")
smiles_list_input = user_input.split(',')
smiles_list_input = [smiles.strip() for smiles in smiles_list_input]

nucl_elec = input("Should we be ranking according to electrophilicity(E) or nucleophilicity(N)?")

get_molecules(smiles_list_input, 400)
