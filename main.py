import os
from core.operations import final_with_bridge, final_no_bridge
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles, Draw


def save_as_xyz(mol, folder, filename):
    """Salva la molecola in formato .xyz nella cartella specificata."""
    mol = Chem.AddHs(mol)
    success = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if success == -1:
        print(f"Errore: generazione 3D fallita per {filename}")
        return False
    filepath = os.path.join(folder, filename)
    rdmolfiles.MolToXYZFile(mol, filepath)
    print(f"Salvato: {filepath}")
    return True

def generate_structures():
    """Genera strutture e le salva in una cartella dedicata."""
    output_folder = "output_structures"  # Cartella di output
    
    for i in range(1, 6):
        # Struttura CON bridge
        mol_with = final_with_bridge()
        if mol_with:
            save_as_xyz(mol_with, output_folder, f"molecola_con_bridge_{i}.xyz")
            Draw.MolToFile(mol_with, os.path.join(output_folder, f"molecola_con_bridge_{i}.png"))
        
        # Struttura SENZA bridge
        mol_without = final_no_bridge()
        if mol_without:
            save_as_xyz(mol_without, output_folder, f"molecola_no_bridge_{i}.xyz")
            Draw.MolToFile(mol_without, os.path.join(output_folder, f"molecola_no_bridge_{i}.png"))

def main():
    generate_structures()
    print("\nGenerazione completata! Controlla la cartella 'output_structures'.")

if __name__ == "__main__":
    main()
