from openeye import oechem

input_file = "C:/Users/sabri/OneDrive/Documentos/UFG/project-BRICS/HsDHODH/Paper 9_ H3D virtual screening/Shape/combined.oeb.gz"
output_file = "C:/Users/sabri/OneDrive/Documentos/UFG/project-BRICS/HsDHODH/Paper 9_ H3D virtual screening/Shape/combined_fixed.oeb.gz"

ifs = oechem.oemolistream()
if not ifs.open(input_file):
    print(f"Error: Cannot open {input_file}")
    exit()

ofs = oechem.oemolostream()
if not ofs.open(output_file):
    print(f"Error: Cannot open {output_file}")
    exit()

reference_query_title = "1D3G" # THIS MUST MATCH YOUR REFERENCE'S CONFORMATION TITLE

print(f"Fixing ROCS_ShapeQuery in {input_file} and saving to {output_file}...")
for mol in ifs.GetOEGraphMols():
    oechem.OESetSDData(mol, "ROCS_ShapeQuery", reference_query_title)
    oechem.OEWriteMolecule(ofs, mol)

ifs.close()
ofs.close()
print("Fixing complete. Try running your main script with the new output_file.")