from rdkit import Chem
from rdkit.Chem import rdchem, Draw
import random
from data.molecules import Smiles, Bridge




def process_random_molecule(Smiles_dict, modify_H=True):
    """
    Process a random molecule from the SMILES dictionary.
    Args:
        Smiles_dict (dict): Dictionary of names and SMILES
        modify_H (bool): If True, removes a random hydrogen
    Returns:
        tuple: (original_mol, modified_mol, atom_info, C_idx)
               where C_idx is the index of the carbon that was bonded to the removed H
               (None if no H was removed)
    """
    # 1. Randomly select a molecule
    smiles = random.choice(list(Smiles_dict.values()))
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")

    # 2. Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # 3. Get atomic information
    atom_info = [{
        "index": atom.GetIdx(),
        "symbol": atom.GetSymbol(),
        "hybridization": str(atom.GetHybridization()),
        "num_H": atom.GetTotalNumHs(),
        "is_aromatic": atom.GetIsAromatic()
    } for atom in mol.GetAtoms()]

    # 4. If requested, modify a random hydrogen
    modified_mol = None
    atom_idx = None
    if modify_H:
        rw_mol = Chem.RWMol(mol)
        h_atoms = [atom for atom in rw_mol.GetAtoms() if atom.GetSymbol() == 'H']

        if h_atoms:
            h_atom = random.choice(h_atoms)
            h_idx = h_atom.GetIdx()
            neighbor = h_atom.GetNeighbors()[0]
            atom_idx = neighbor.GetIdx()
            print(f"Removed hydrogen {h_idx} bonded to {neighbor.GetSymbol()}{atom_idx}")
            rw_mol.RemoveAtom(h_idx)
            modified_mol = rw_mol.GetMol()

    return mol, modified_mol, atom_info, atom_idx

def combine_molecules(mol1, mol2, atom_idx1=None, atom_idx2=None, bond_type='SINGLE', sanitize=True):
    """
    Combine two molecules into a single entity, optionally creating a bond between them.
    
    Args:
        mol1 (rdkit.Chem.Mol): First molecule
        mol2 (rdkit.Chem.Mol): Second molecule
        atom_idx1 (int, optional): Index of atom in first molecule for bonding
        atom_idx2 (int, optional): Index of atom in second molecule for bonding
        bond_type (str): Bond type ('SINGLE', 'DOUBLE', 'TRIPLE', 'AROMATIC')
        sanitize (bool): If True, performs RDKit sanitization
        
    Returns:
        rdkit.Chem.Mol: Combined molecule
    """
    # Create a combined molecule
    combined = Chem.RWMol(Chem.CombineMols(mol1, mol2))

    # If atoms are specified for bonding
    if atom_idx1 is not None and atom_idx2 is not None:
        # Convert bond type
        bond_type_dict = {
            'SINGLE': rdchem.BondType.SINGLE,
            'DOUBLE': rdchem.BondType.DOUBLE,
            'TRIPLE': rdchem.BondType.TRIPLE,
            'AROMATIC': rdchem.BondType.AROMATIC
        }
        
        if bond_type not in bond_type_dict:
            raise ValueError(f"Invalid bond type: {bond_type}")
            
        adjusted_atom_idx2 = atom_idx2 + mol1.GetNumAtoms()

        # Verify atoms exist
        if atom_idx1 >= mol1.GetNumAtoms():
            raise ValueError(f"Index {atom_idx1} invalid for first molecule")
        if atom_idx2 >= mol2.GetNumAtoms():
            raise ValueError(f"Index {atom_idx2} invalid for second molecule")

        # Add the bond
        combined.AddBond(atom_idx1, adjusted_atom_idx2, bond_type_dict[bond_type])

        # Remove hydrogens if needed (to maintain correct valence)
        atom1 = combined.GetAtomWithIdx(atom_idx1)
        atom2 = combined.GetAtomWithIdx(adjusted_atom_idx2)

        if atom1.GetSymbol() != 'H':
            atom1.SetNumExplicitHs(max(0, atom1.GetNumExplicitHs() - 1))
        if atom2.GetSymbol() != 'H':
            atom2.SetNumExplicitHs(max(0, atom2.GetNumExplicitHs() - 1))

    # Convert to regular molecule
    result = combined.GetMol()

    # Sanitization (optional but recommended)
    if sanitize:
        try:
            Chem.SanitizeMol(result)
        except:
            print("Warning: Sanitization failed - the resulting molecule might be unstable")

    return result

def process_bridge_molecule(combined_with_bridge, modify_C=True):
    """
    Process a molecule with bridge, now formed
    Args: combined_with_bridge : the molecule formed by combining an aromatic molecule with a bridge
          modify_C (bool): If True, removes a random hydrogen from an SP3 carbon atom
    Returns:
        tuple: (modified_mol, atom_info, C_idx)
               where C_idx is the index of the carbon that was bonded to the removed H
               (None if no H was removed)
    """
    mol = Chem.AddHs(combined_with_bridge)

    # 1. Get atomic information
    atom_info = [{
        "index": atom.GetIdx(),
        "symbol": atom.GetSymbol(),
        "hybridization": str(atom.GetHybridization()),
        "num_H": atom.GetTotalNumHs(),
        "is_aromatic": atom.GetIsAromatic()
    } for atom in mol.GetAtoms()]

    # 2. If requested, find a SP3 Carbon atom and remove one of its hydrogens
    modified_mol = None
    C_idx = None
    if modify_C:
        rw_mol = Chem.RWMol(mol)
        SP3_C = [atom for atom in rw_mol.GetAtoms() 
                 if str(atom.GetHybridization()) == "SP3" 
                 and atom.GetSymbol() == 'C' 
                 ]

        print("SP3 Carbon atoms found:", [(atom.GetIdx(), atom.GetSymbol(), atom.GetTotalNumHs()) for atom in SP3_C])

        if SP3_C:
            C_atom = random.choice(SP3_C)
            C_idx = C_atom.GetIdx()
            # Find a hydrogen bonded to this carbon
            h_atoms = [n for n in C_atom.GetNeighbors() if n.GetSymbol() == 'H']
            if h_atoms:
                h_atom = random.choice(h_atoms)
                h_idx = h_atom.GetIdx()
                print(f"Removed hydrogen {h_idx} bonded to {C_atom.GetSymbol()}{C_idx}")
                rw_mol.RemoveAtom(h_idx)
                modified_mol = rw_mol.GetMol()

    return modified_mol if modified_mol is not None else mol, atom_info, C_idx

