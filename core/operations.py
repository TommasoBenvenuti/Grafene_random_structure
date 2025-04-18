from core import process_random_molecule, combine_molecules, process_bridge_molecule
from data.molecules import Smiles, Bridge 



## ============================== BRIDGE =============================== ##
## Here i create a molecule with bridge. Two condensate ring, and bridge ##
## ===================================================================== ##


def final_with_bridge():
    ## Process molecules. The program takes a random molecule from the list, remove H, label the position. The same does for the bridge molecule
    mol1, mol1_modified, info1, atom_idx1 = process_random_molecule(Smiles)
    bridge, bridge_modified, info_bridge, bridge_idx = process_random_molecule(Bridge)

    ## Use modified molecules if available, if everything is working in core.py
    mol1_to_combine = mol1_modified if mol1_modified is not None else mol1
    bridge_to_combine = bridge_modified if bridge_modified is not None else bridge

    ## Combine with bond if hydrogens were removed
    if atom_idx1 is not None and bridge_idx is not None:
        combined_with_bridge = combine_molecules(mol1_to_combine, bridge_to_combine, atom_idx1, bridge_idx, 'SINGLE')
        print("First combination successful!")
    else:
        print("Couldn't combine with bond - no hydrogens were removed from one or both molecules")
        combined_with_bridge = combine_molecules(mol1_to_combine, bridge_to_combine)

    ## Process the combined molecule. This molecule is clearly an intermediate (here you have a ring condensate molecule, bonded to a "bridge"). C_idx is the position labelled of Carbon (for the bond!!)
    Mol_with_bridge, info_intermol, C_idx = process_bridge_molecule(combined_with_bridge)

    ## Process second molecule. Another condensated ring mol.
    mol2, mol2_modified, info2, atom_idx2 = process_random_molecule(Smiles)
    mol2_to_combine = mol2_modified if mol2_modified is not None else mol2

    # Final combination
    if atom_idx2 is not None and C_idx is not None:
        combined_bridge = combine_molecules(mol2_to_combine, Mol_with_bridge, atom_idx2, C_idx, 'SINGLE')
        print("Final combination successful!")
    else:
        print("Couldn't combine with bond - no hydrogens were removed from one or both molecules")
        combined_bridge = combine_molecules(mol2_to_combine, Mol_with_bridge)

    return combined_bridge
## ======================= NO BRIDGE  =============================== ##
## Here i create a molecule without bridge. Just two condensate ring. ##
## ================================================================== ##


def final_no_bridge():
## Process molecules. The program takes a random molecule from the list, remove H, label the position. The same does for the bridge molecule
    mol1, mol1_modified, info1, atom_idx1 = process_random_molecule(Smiles)
    mol2, mol2_modified, info2, atom_idx2 = process_random_molecule(Smiles)

   
## Use modified molecules if available, if everything is working in core.py
    mol1_to_combine = mol1_modified if mol1_modified is not None else mol1
    mol2_to_combine = mol2_modified if mol2_modified is not None else mol2

## Combine with bond if hydrogens were removed
    if atom_idx1 is not None and atom_idx2 is not None:
       combined = combine_molecules(mol1_to_combine, mol2_to_combine, atom_idx1, atom_idx2, 'SINGLE')
       print("First combination successful!")
    else:
       print("Couldn't combine with bond - no hydrogens were removed from one or both molecules")
       combined = combine_molecules(mol1_to_combine, mol2_to_combine) 

    return combined 
