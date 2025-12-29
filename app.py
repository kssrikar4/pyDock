import streamlit as st
import os
import sys
import subprocess
import tempfile
import shutil
import traceback
import re
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

DEFAULT_PADDING = 5.0
VINA_EXHAUSTIVENESS = 8
NUM_POSES = 5
MIN_SEARCH_SPACE = 20.0

def parse_pdb_coordinates(pdb_content):
    coordinates = []
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coordinates.append((x, y, z))
            except (ValueError, IndexError):
                continue
    return coordinates

def calculate_binding_box(pdb_content, padding=DEFAULT_PADDING):
    coordinates = parse_pdb_coordinates(pdb_content)
    if not coordinates:
        raise ValueError("No valid atom coordinates found in PDB file")
    coords_array = np.array(coordinates)
    min_coords = coords_array.min(axis=0)
    max_coords = coords_array.max(axis=0)
    center = (min_coords + max_coords) / 2.0
    size = max_coords - min_coords + (2 * padding)
    size = np.maximum(size, MIN_SEARCH_SPACE)
    return {
        'center_x': float(center[0]),
        'center_y': float(center[1]),
        'center_z': float(center[2]),
        'size_x': float(size[0]),
        'size_y': float(size[1]),
        'size_z': float(size[2]),
    }

def read_ligand(ligand_path):
    ligand_path = str(ligand_path)
    if ligand_path.lower().endswith('.sdf'):
        mol = next(Chem.SDMolSupplier(ligand_path, removeHs=False))
        if mol is None:
            raise ValueError("Failed to read SDF file or molecule is invalid")
    elif ligand_path.lower().endswith('.mol2'):
        mol = Chem.MolFromMol2File(ligand_path, removeHs=False)
        if mol is None:
            raise ValueError("Failed to read MOL2 file or molecule is invalid")
    else:
        raise ValueError("Ligand format must be SDF or MOL2")
    return mol

def validate_ligand(mol):
    if mol is None:
        raise ValueError("Invalid molecule object")
    mw = Descriptors.MolWt(mol)
    if mw < 50 or mw > 5000:
        raise ValueError(f"Molecular weight {mw:.1f} may be problematic for docking")
    n_atoms = mol.GetNumAtoms()
    if n_atoms < 1:
        raise ValueError("Ligand has no atoms")
    return True

def prepare_protein_meeko(pdb_path, output_pdbqt_path):
    try:
        python_dir = Path(sys.executable).parent
        script_path = python_dir / "Scripts" / "mk_prepare_receptor.py"
        if not script_path.exists():
            script_path = python_dir / "Scripts" / "mk_prepare_receptor"
        if script_path.exists():
            cmd = [sys.executable, str(script_path), "--read_pdb", str(pdb_path), "-o", str(output_pdbqt_path)]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            if result.returncode == 0 and os.path.exists(output_pdbqt_path):
                clean_pdbqt_file(output_pdbqt_path)
                return output_pdbqt_path
        with open(pdb_path, 'r') as f:
            pdb_content = f.read()
        pdbqt_content = convert_pdb_to_pdbqt_simple(pdb_content)
        with open(output_pdbqt_path, 'w') as f:
            f.write(pdbqt_content)
        if not os.path.exists(output_pdbqt_path):
            raise RuntimeError("Receptor PDBQT file was not created")
        return output_pdbqt_path
    except Exception as e:
        raise RuntimeError(f"Protein preparation failed: {str(e)}")

def clean_pdbqt_file(pdbqt_path):
    with open(pdbqt_path, 'r') as f:
        lines = f.readlines()
    cleaned_lines = []
    for line in lines:
        if not line.startswith('MODEL') and not line.startswith('ENDMDL'):
            cleaned_lines.append(line)
    with open(pdbqt_path, 'w') as f:
        f.writelines(cleaned_lines)

def convert_pdb_to_pdbqt_simple(pdb_content):
    pdbqt_lines = []
    atom_types = {
        'C': 'C', 'N': 'N', 'O': 'O', 'S': 'S',
        'P': 'P', 'H': 'HD', 'F': 'F', 'Cl': 'Cl',
        'Br': 'Br', 'I': 'I', 'Fe': 'Fe', 'Mg': 'Mg',
        'Ca': 'Ca', 'Zn': 'Zn', 'Mn': 'Mn'
    }
    for line in pdb_content.split('\n'):
        if line.startswith('MODEL') or line.startswith('ENDMDL'):
            continue
        if line.startswith('ATOM') or line.startswith('HETATM'):
            if len(line) > 77 and line[76:78].strip():
                element = line[76:78].strip()
            else:
                atom_name = line[12:16].strip()
                element = ''.join([c for c in atom_name if c.isalpha()])
                if len(element) >= 2 and element[1].islower():
                    element = element[0].upper() + element[1]
                elif len(element) >= 1:
                    element = element[0].upper()
                else:
                    element = 'C'
            ad_type = atom_types.get(element, 'C')
            if len(line) >= 66:
                base_line = line[:66]
            else:
                base_line = line.ljust(66)
            pdbqt_line = f"{base_line}    {0.000:6.3f}  {ad_type:>2s}"
            pdbqt_lines.append(pdbqt_line)
        elif line.startswith('TER'):
            pdbqt_lines.append(line)
        elif line.startswith('END') and not line.startswith('ENDMDL'):
            pdbqt_lines.append(line)
    if pdbqt_lines and not pdbqt_lines[-1].startswith('END'):
        pdbqt_lines.append('END')
    return '\n'.join(pdbqt_lines)

def prepare_ligand_meeko(sdf_path, output_pdbqt_path):
    try:
        ligand_path = str(sdf_path)
        if ligand_path.lower().endswith('.sdf'):
            mol = next(Chem.SDMolSupplier(ligand_path, removeHs=False))
        elif ligand_path.lower().endswith('.mol2'):
            mol = Chem.MolFromMol2File(ligand_path, removeHs=False)
        else:
            raise ValueError("Ligand format must be SDF or MOL2")
        if mol is None:
            raise RuntimeError(f"Failed to read ligand file: {sdf_path}")
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)
        if not mol_setups:
            raise RuntimeError("Meeko ligand preparation returned no results")
        setup = mol_setups[0]
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
        if not is_ok:
            raise RuntimeError(f"Meeko ligand preparation failed: {error_msg}")
        with open(output_pdbqt_path, 'w') as f:
            f.write(pdbqt_string)
        if not os.path.exists(output_pdbqt_path):
            raise RuntimeError("Meeko did not generate ligand PDBQT file")
        return output_pdbqt_path
    except ImportError as e:
        raise RuntimeError(f"Meeko not properly installed: {e}")

def run_vina_docking(vina_exe, protein_pdbqt, ligand_pdbqt, grid_params, output_pdbqt):
    try:
        cmd = [
            vina_exe,
            '--receptor', str(protein_pdbqt),
            '--ligand', str(ligand_pdbqt),
            '--center_x', f"{grid_params['center_x']:.3f}",
            '--center_y', f"{grid_params['center_y']:.3f}",
            '--center_z', f"{grid_params['center_z']:.3f}",
            '--size_x', f"{grid_params['size_x']:.3f}",
            '--size_y', f"{grid_params['size_y']:.3f}",
            '--size_z', f"{grid_params['size_z']:.3f}",
            '--exhaustiveness', str(VINA_EXHAUSTIVENESS),
            '--out', str(output_pdbqt),
            '--num_modes', str(NUM_POSES)
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            raise RuntimeError(f"Vina docking failed:\n{result.stderr}")
        if not os.path.exists(output_pdbqt):
            raise RuntimeError("Vina did not generate output PDBQT")
        return output_pdbqt
    except FileNotFoundError:
        raise RuntimeError(
            f"AutoDock Vina executable not found at: {vina_exe}\n"
            "Download from: https://vina.scripps.edu/download/"
        )

def parse_vina_output(pdbqt_file):
    poses = []
    current_model = []
    current_affinity = None
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith('MODEL'):
                current_model = []
                current_affinity = None
            elif line.startswith('REMARK VINA RESULT'):
                parts = line.split()
                try:
                    current_affinity = float(parts[3])
                except (IndexError, ValueError):
                    current_affinity = None
            elif line.startswith('ENDMDL'):
                if current_model and current_affinity is not None:
                    poses.append({
                        'affinity': current_affinity,
                        'pdblines': current_model
                    })
                current_model = []
                current_affinity = None
            elif line.startswith(('ATOM', 'HETATM', 'CONECT')):
                current_model.append(line)
    poses.sort(key=lambda x: x['affinity'])
    return poses[:NUM_POSES]

def extract_top_poses(vina_output_pdbqt, protein_pdb, output_dir):
    poses = parse_vina_output(vina_output_pdbqt)
    with open(protein_pdb, 'r') as f:
        protein_lines = f.readlines()
    protein_atoms = [line for line in protein_lines if line.startswith(('ATOM', 'HETATM'))]
    complex_files = []
    for idx, pose in enumerate(poses, 1):
        complex_pdb = os.path.join(output_dir, f"complex_pose{idx}.pdb")
        with open(complex_pdb, 'w') as f:
            f.write(f"REMARK   Blind Docking Complex - Pose {idx}\n")
            f.write(f"REMARK   Binding Affinity: {pose['affinity']:.2f} kcal/mol\n")
            f.writelines(protein_atoms)
            f.write("TER\n")
            ligand_lines = convert_ligand_to_hetatm(pose['pdblines'])
            f.writelines(ligand_lines)
            f.write("END\n")
        complex_files.append({
            'path': complex_pdb,
            'pose': idx,
            'affinity': pose['affinity']
        })
    return complex_files

def convert_ligand_to_hetatm(pdbqt_lines):
    hetatm_lines = []
    for line in pdbqt_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            if len(line) >= 66:
                pdb_line = "HETATM" + line[6:66] + "\n"
            else:
                pdb_line = "HETATM" + line[6:].rstrip() + "\n"
            hetatm_lines.append(pdb_line)
        elif line.startswith('CONECT'):
            hetatm_lines.append(line)
    return hetatm_lines

def analyze_interactions(complex_pdb_path):
    try:
        with open(complex_pdb_path, 'r') as f:
            pdb_lines = f.readlines()
        protein_atoms = []
        ligand_atoms = []
        for line in pdb_lines:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                residue = line[17:20].strip()
                resnum = line[22:26].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                protein_atoms.append({
                    'name': atom_name,
                    'residue': f"{residue}{resnum}",
                    'coords': np.array([x, y, z]),
                    'element': element
                })
            elif line.startswith('HETATM'):
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                ligand_atoms.append({
                    'name': atom_name,
                    'coords': np.array([x, y, z]),
                    'element': element
                })
        if not ligand_atoms:
            st.warning("⚠️ No ligand atoms found in complex")
            return None
        interactions = {
            'hydrogen_bonds': [],
            'hydrophobic': [],
            'salt_bridges': [],
            'pi_stacking': [],
            'pi_cation': [],
            'halogen': []
        }
        residue_interactions = {}
        for lig_atom in ligand_atoms:
            lig_pos = lig_atom['coords']
            lig_elem = lig_atom['element']
            for prot_atom in protein_atoms:
                prot_pos = prot_atom['coords']
                prot_elem = prot_atom['element']
                residue = prot_atom['residue']
                distance = np.linalg.norm(lig_pos - prot_pos)
                if distance < 3.5:
                    if (lig_elem in ['N', 'O'] and prot_elem in ['N', 'O']):
                        key = f"hbond_{residue}"
                        if key not in residue_interactions or distance < residue_interactions[key]['data']['distance']:
                            residue_interactions[key] = {
                                'type': 'hydrogen_bonds',
                                'data': {
                                    'donor': tuple(prot_pos),
                                    'acceptor': tuple(lig_pos),
                                    'residue': residue,
                                    'distance': distance
                                }
                            }
                if distance < 4.5:
                    if lig_elem == 'C' and prot_elem == 'C':
                        key = f"hydrophobic_{residue}"
                        if key not in residue_interactions or distance < residue_interactions[key]['data']['distance']:
                            residue_interactions[key] = {
                                'type': 'hydrophobic',
                                'data': {
                                    'ligand_coords': tuple(lig_pos),
                                    'protein_coords': tuple(prot_pos),
                                    'residue': residue,
                                    'distance': distance
                                }
                            }
                if distance < 3.5:
                    if lig_elem in ['F', 'Cl', 'Br', 'I'] and prot_elem in ['O', 'N']:
                        key = f"halogen_{residue}"
                        if key not in residue_interactions or distance < residue_interactions[key]['data']['distance']:
                            residue_interactions[key] = {
                                'type': 'halogen',
                                'data': {
                                    'donor': tuple(lig_pos),
                                    'acceptor': tuple(prot_pos),
                                    'residue': residue,
                                    'distance': distance
                                }
                            }
        for inter in residue_interactions.values():
            interactions[inter['type']].append(inter['data'])
        return interactions
    except Exception as e:
        st.warning(f"⚠️ Interaction analysis failed: {str(e)}")
        st.code(traceback.format_exc())
        return None

def generate_hbond_js(hbonds):
    js_lines = []
    for hb in hbonds:
        donor = hb['donor']
        acceptor = hb['acceptor']
        js_lines.append(f"""
        viewer.addLine({{
            dashed: true,
            start: {{x: {donor[0]:.3f}, y: {donor[1]:.3f}, z: {donor[2]:.3f}}},
            end: {{x: {acceptor[0]:.3f}, y: {acceptor[1]:.3f}, z: {acceptor[2]:.3f}}},
            color: 'cyan',
            linewidth: 3
        }});
        """)
    return '\n'.join(js_lines)

def generate_hydrophobic_js(hydrophobic):
    js_lines = []
    for hc in hydrophobic:
        lig = hc['ligand_coords']
        prot = hc['protein_coords']
        js_lines.append(f"""
        viewer.addLine({{
            dashed: true,
            start: {{x: {lig[0]:.3f}, y: {lig[1]:.3f}, z: {lig[2]:.3f}}},
            end: {{x: {prot[0]:.3f}, y: {prot[1]:.3f}, z: {prot[2]:.3f}}},
            color: 'gray',
            linewidth: 3
        }});
        """)
    return '\n'.join(js_lines)

def generate_saltbridge_js(saltbridges):
    js_lines = []
    for sb in saltbridges:
        lig = sb['ligand_coords']
        prot = sb['protein_coords']
        js_lines.append(f"""
        viewer.addLine({{
            dashed: true,
            start: {{x: {lig[0]:.3f}, y: {lig[1]:.3f}, z: {lig[2]:.3f}}},
            end: {{x: {prot[0]:.3f}, y: {prot[1]:.3f}, z: {prot[2]:.3f}}},
            color: 'yellow',
            linewidth: 3
        }});
        """)
    return '\n'.join(js_lines)

def generate_pistacking_js(pistacking):
    js_lines = []
    for ps in pistacking:
        lig = ps['ligand_coords']
        prot = ps['protein_coords']
        js_lines.append(f"""
        viewer.addLine({{
            dashed: true,
            start: {{x: {lig[0]:.3f}, y: {lig[1]:.3f}, z: {lig[2]:.3f}}},
            end: {{x: {prot[0]:.3f}, y: {prot[1]:.3f}, z: {prot[2]:.3f}}},
            color: 'magenta',
            linewidth: 3
        }});
        """)
    return '\n'.join(js_lines)

def generate_pication_js(pication):
    js_lines = []
    for pc in pication:
        lig = pc['ligand_coords']
        prot = pc['protein_coords']
        js_lines.append(f"""
        viewer.addLine({{
            dashed: true,
            start: {{x: {lig[0]:.3f}, y: {lig[1]:.3f}, z: {lig[2]:.3f}}},
            end: {{x: {prot[0]:.3f}, y: {prot[1]:.3f}, z: {prot[2]:.3f}}},
            color: 'orange',
            linewidth: 3
        }});
        """)
    return '\n'.join(js_lines)

def generate_halogen_js(halogen):
    js_lines = []
    for halo in halogen:
        donor = halo['donor']
        acceptor = halo['acceptor']
        js_lines.append(f"""
        viewer.addLine({{
            dashed: true,
            start: {{x: {donor[0]:.3f}, y: {donor[1]:.3f}, z: {donor[2]:.3f}}},
            end: {{x: {acceptor[0]:.3f}, y: {acceptor[1]:.3f}, z: {acceptor[2]:.3f}}},
            color: 'green',
            linewidth: 3
        }});
        """)
    return '\n'.join(js_lines)

def generate_residue_labels(interactions):
    residues = set()
    positions = {}
    for interaction_type in interactions.values():
        for inter in interaction_type:
            if 'residue' in inter:
                residue = inter['residue']
                residues.add(residue)
                if 'protein_coords' in inter:
                    positions[residue] = inter['protein_coords']
                elif 'acceptor' in inter:
                    positions[residue] = inter['acceptor']
    js_lines = []
    for residue in residues:
        if residue in positions:
            pos = positions[residue]
            js_lines.append(f"""
            viewer.addLabel('{residue}', {{
                position: {{x: {pos[0]:.3f}, y: {pos[1]:.3f}, z: {pos[2]:.3f}}},
                backgroundColor: 'black',
                fontColor: 'white',
                fontSize: 12,
                backgroundOpacity: 0.7
            }});
            """)
    return '\n'.join(js_lines)

def read_pdb_to_string(pdb_path):
    with open(pdb_path, 'r') as f:
        return f.read()

def create_interactive_view_with_interactions(pdb_string, interactions):
    html_code = f"""
    <div id="container" style="height: 600px; width: 100%; position: relative;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        let viewer = $3Dmol.createViewer("container");
        let pdbData = `{pdb_string}`;
        viewer.addModel(pdbData, 'pdb');
        viewer.setStyle({{not: {{hetflag: true}}}}, {{cartoon: {{color: 'lightgray', radius: 0.2}}}});
        viewer.setStyle({{hetflag: true}}, {{stick: {{colorscheme: 'greenCarbon', radius: 0.3}}}});
        viewer.addStyle({{hetflag: true}}, {{sphere: {{colorscheme: 'greenCarbon', radius: 0.3}}}});
        {generate_hbond_js(interactions.get('hydrogen_bonds', []))}
        {generate_hydrophobic_js(interactions.get('hydrophobic', []))}
        {generate_saltbridge_js(interactions.get('salt_bridges', []))}
        {generate_pistacking_js(interactions.get('pi_stacking', []))}
        {generate_pication_js(interactions.get('pi_cation', []))}
        {generate_halogen_js(interactions.get('halogen', []))}
        {generate_residue_labels(interactions)}
        viewer.zoomTo();
        viewer.render();
    </script>
    """
    return html_code

def create_surface_view(pdb_string):
    html_code = f"""
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <div id="viewer" style="width: 100%; height: 600px; position: relative;"></div>
    <script>
        let viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor: "white"}});
        let data = `{pdb_string}`;
        viewer.addModel(data, "pdb");
        viewer.setStyle({{not: {{hetflag: true}}}}, {{cartoon: {{color: "spectrum"}}}});
        viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity: 0.8, color: "white"}}, {{not: {{hetflag: true}}}});
        viewer.setStyle({{hetflag: true}}, {{
            stick: {{colorscheme: "greenCarbon", radius: 0.3}},
            sphere: {{colorscheme: "greenCarbon", radius: 0.5}}
        }});
        viewer.zoomTo();
        viewer.render();
    </script>
    """
    return html_code

def create_cartoon_view(pdb_string):
    html_code = f"""
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <div id="viewer" style="width: 100%; height: 600px; position: relative;"></div>
    <script>
        let viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor: "white"}});
        let data = `{pdb_string}`;
        viewer.addModel(data, "pdb");
        viewer.setStyle({{not: {{hetflag: true}}}}, {{cartoon: {{color: "spectrum", thickness: 0.8}}}});
        viewer.setStyle({{hetflag: true}}, {{
            stick: {{colorscheme: "greenCarbon", radius: 0.3}},
            sphere: {{colorscheme: "greenCarbon", radius: 0.5}}
        }});
        viewer.zoomTo();
        viewer.render();
    </script>
    """
    return html_code

def create_stick_view(pdb_string):
    html_code = f"""
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <div id="viewer" style="width: 100%; height: 600px; position: relative;"></div>
    <script>
        let viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor: "white"}});
        let data = `{pdb_string}`;
        let model = viewer.addModel(data, "pdb");
        viewer.setStyle({{chain: ["A", "B", "C", "D", "E", "F"]}}, {{cartoon: {{color: "lightgray", opacity: 0.7}}}});
        viewer.setStyle({{hetflag: true}}, {{
            stick: {{colorscheme: "greenCarbon", radius: 0.3}},
            sphere: {{colorscheme: "greenCarbon", radius: 0.5}}
        }});
        viewer.zoomTo({{hetflag: true}});
        viewer.render();
    </script>
    """
    return html_code

def initialize_session_state():
    if 'vina_path' not in st.session_state:
        st.session_state.vina_path = None
    if 'docking_complete' not in st.session_state:
        st.session_state.docking_complete = False
    if 'complex_files' not in st.session_state:
        st.session_state.complex_files = []
    if 'temp_dir' not in st.session_state:
        st.session_state.temp_dir = None

def main():
    st.set_page_config(
        page_title="Protein-Ligand Blind Docking",
        page_icon="🧬",
        layout="wide"
    )
    initialize_session_state()
    st.title("🧬 Protein-Ligand Blind Docking Pipeline")
    st.markdown("""
    **Fully Automated Blind Docking Workflow**
    - Accepts protein (PDB) and ligand (SDF/MOL2) files
    - Automatically computes docking grid parameters
    - Uses AutoDock Vina for molecular docking
    - Outputs top 5 binding poses with interactive 3D visualization
    - Supports surface and cartoon representation modes
    **Reference:** CBDock2-style blind docking (covers entire protein surface)
    """)
    st.divider()
    st.header("Step 1: Configure AutoDock Vina")
    st.info(
        "⚠️ AutoDock Vina must be installed on your system.\n"
        "Download from: https://vina.scripps.edu/download/"
    )
    vina_path_input = st.text_input(
        "Enter the full path to AutoDock Vina executable:",
        value=st.session_state.vina_path or "",
        placeholder="/usr/bin/vina (Linux) or C:\\Program Files\\vina.exe (Windows)"
    )
    if vina_path_input:
        st.session_state.vina_path = vina_path_input
        if os.path.exists(vina_path_input):
            st.success("✅ Vina executable found!")
        else:
            st.error(f"❌ File not found at: {vina_path_input}")
    st.divider()
    st.header("Step 2: Upload Input Files")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Protein (PDB)")
        protein_file = st.file_uploader(
            "Upload protein PDB file:",
            type=['pdb'],
            key='protein_upload'
        )
    with col2:
        st.subheader("Ligand (SDF or MOL2)")
        ligand_file = st.file_uploader(
            "Upload ligand file:",
            type=['sdf', 'mol2'],
            key='ligand_upload'
        )
    st.divider()
    st.header("Step 3: Execute Docking")
    if st.button("▶️ Run Blind Docking", type="primary", use_container_width=True):
        if not st.session_state.vina_path:
            st.error("❌ Please specify AutoDock Vina executable path")
            st.stop()
        if not protein_file:
            st.error("❌ Please upload protein PDB file")
            st.stop()
        if not ligand_file:
            st.error("❌ Please upload ligand file (SDF or MOL2)")
            st.stop()
        if not os.path.exists(st.session_state.vina_path):
            st.error(f"❌ Vina executable not found at: {st.session_state.vina_path}")
            st.stop()
        temp_dir = tempfile.mkdtemp(prefix='docking_')
        st.session_state.temp_dir = temp_dir
        try:
            with st.spinner("🔬 Processing and preparing files..."):
                protein_path = os.path.join(temp_dir, 'protein.pdb')
                with open(protein_path, 'w') as f:
                    f.write(protein_file.getvalue().decode('utf-8'))
                ligand_path = os.path.join(temp_dir, f'ligand.{ligand_file.name.split(".")[-1]}')
                with open(ligand_path, 'wb') as f:
                    f.write(ligand_file.getvalue())
                st.write("✓ Files saved")
                with open(protein_path, 'r') as f:
                    protein_content = f.read()
                st.write("✓ Calculating blind docking grid parameters...")
                grid_params = calculate_binding_box(protein_content, padding=DEFAULT_PADDING)
                st.write(f"""
                **Grid Parameters (Automatic Calculation):**
                - Center: ({grid_params['center_x']:.2f}, {grid_params['center_y']:.2f}, {grid_params['center_z']:.2f})
                - Size: {grid_params['size_x']:.2f} × {grid_params['size_y']:.2f} × {grid_params['size_z']:.2f} Ų
                """)
                st.write("✓ Validating ligand structure (RDKit)...")
                mol = read_ligand(ligand_path)
                validate_ligand(mol)
                st.write(f"  Molecular weight: {Descriptors.MolWt(mol):.2f}")
                st.write(f"  Atom count: {mol.GetNumAtoms()}")
                st.write("✓ Preparing protein with Meeko...")
                protein_pdbqt = os.path.join(temp_dir, 'protein.pdbqt')
                prepare_protein_meeko(protein_path, protein_pdbqt)
                st.write("✓ Preparing ligand with Meeko...")
                ligand_pdbqt = os.path.join(temp_dir, 'ligand.pdbqt')
                prepare_ligand_meeko(ligand_path, ligand_pdbqt)
                st.write("✓ Running AutoDock Vina blind docking...")
                vina_output = os.path.join(temp_dir, 'docked_output.pdbqt')
                run_vina_docking(
                    st.session_state.vina_path,
                    protein_pdbqt,
                    ligand_pdbqt,
                    grid_params,
                    vina_output
                )
                st.write("✓ Extracting top 5 binding poses...")
                complex_files = extract_top_poses(vina_output, protein_path, temp_dir)
                st.session_state.complex_files = complex_files
                st.success("✅ Docking completed successfully!")
                st.session_state.docking_complete = True
        except Exception as e:
            st.error(f"❌ Error during docking:\n{str(e)}")
            if st.session_state.temp_dir and os.path.exists(st.session_state.temp_dir):
                shutil.rmtree(st.session_state.temp_dir)
            st.stop()
    st.divider()
    if st.session_state.docking_complete and st.session_state.complex_files:
        st.header("Step 4: Interactive 3D Visualization")
        pose_options = {
            f"Pose {cplx['pose']} (ΔG = {cplx['affinity']:.2f} kcal/mol)": cplx['path']
            for cplx in st.session_state.complex_files
        }
        selected_pose = st.selectbox(
            "Select binding pose to visualize:",
            options=pose_options.keys()
        )
        selected_pdb = pose_options[selected_pose]
        st.subheader("Visualization Mode")
        vizmode = st.radio(
            "Choose representation:",
            options=[
                "Surface Representation",
                "3D Cartoon Representation",
                "Interaction Analysis (with bonds & labels)"
            ],
            horizontal=False
        )
        pdb_content = read_pdb_to_string(selected_pdb)
        if vizmode == "Interaction Analysis (with bonds & labels)":
            with st.spinner("🔍 Analyzing protein-ligand interactions..."):
                interactions = analyze_interactions(selected_pdb)
            if interactions:
                st.subheader("Interaction Summary")
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Hydrogen Bonds", len(interactions['hydrogen_bonds']))
                    st.metric("Hydrophobic", len(interactions['hydrophobic']))
                with col2:
                    st.metric("Salt Bridges", len(interactions['salt_bridges']))
                    st.metric("π-Stacking", len(interactions['pi_stacking']))
                with col3:
                    st.metric("π-Cation", len(interactions['pi_cation']))
                    st.metric("Halogen Bonds", len(interactions['halogen']))
                st.divider()
                html_viz = create_interactive_view_with_interactions(pdb_content, interactions)
                st.components.v1.html(html_viz, height=650)
                st.info("""
                **Interaction Color Code:**
                - 🔵 **Cyan:** Hydrogen bonds
                - ⚪ **Gray:** Hydrophobic contacts
                - 🟡 **Yellow:** Salt bridges
                - 🟣 **Magenta:** π-π stacking
                - 🟠 **Orange:** π-Cation interactions
                - 🟢 **Green:** Halogen bonds
                **Residue Labels:** Black boxes show interacting amino acid residues
                """)
            else:
                st.warning("⚠️ Interaction analysis failed. Try another visualization mode.")
                html_viz = create_cartoon_view(pdb_content)
                st.components.v1.html(html_viz, height=650)
        elif vizmode == "Surface Representation":
            html_viz = create_surface_view(pdb_content)
            st.components.v1.html(html_viz, height=650)
        elif vizmode == "3D Cartoon Representation":
            html_viz = create_cartoon_view(pdb_content)
            st.components.v1.html(html_viz, height=650)
        st.info("""
        **Controls:**
        - **Rotate:** Click and drag with mouse
        - **Zoom:** Scroll wheel or pinch on trackpad
        - **Pan:** Right-click and drag
        - **Reset:** Double-click
        """)
        st.divider()
        st.subheader("Download Results")
        for idx, cplx in enumerate(st.session_state.complex_files):
            with st.expander(f"📥 Pose {cplx['pose']} (ΔG = {cplx['affinity']:.2f})"):
                with open(cplx['path'], 'r') as f:
                    pdb_content_dl = f.read()
                st.download_button(
                    label=f"Download Pose {cplx['pose']} PDB",
                    data=pdb_content_dl,
                    file_name=f"pose_{cplx['pose']}.pdb",
                    mime="text/plain"
                )
        st.subheader("Docking Summary")
        summary_data = {
            "Pose": [cplx['pose'] for cplx in st.session_state.complex_files],
            "Binding Affinity (kcal/mol)": [f"{cplx['affinity']:.2f}" for cplx in st.session_state.complex_files],
        }
        st.table(summary_data)
        if st.button("🗑️ Clear Results and Start Over"):
            if st.session_state.temp_dir and os.path.exists(st.session_state.temp_dir):
                shutil.rmtree(st.session_state.temp_dir)
            st.session_state.docking_complete = False
            st.session_state.complex_files = []
            st.rerun()

if __name__ == "__main__":
    main()