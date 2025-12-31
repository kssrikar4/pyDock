# pyDock

**Python Implementation for Protein-Ligand Blind Docking**

A fully automated blind docking pipeline that requires only three user inputs: protein PDB file, ligand file (SDF/MOL2), and AutoDock Vina executable path. Everything else grid parameters, binding box calculation, structure preparation is computed automatically.

## Short Video Tutorial

https://github.com/user-attachments/assets/3e3318de-2420-4542-a069-befae9d3e2ee

## What This Does

| Feature | Details |
|---------|---------|
| **Input** | Protein (PDB) + Ligand (SDF/MOL2) |
| **Processing** | Automatic grid calculation, structure validation, PDBQT preparation |
| **Docking Engine** | AutoDock Vina (blind docking across entire protein) |
| **Output** | Top 5 binding poses with binding affinity scores |
| **Visualization** | Interactive 3D viewer with surface, cartoon modes & Interaction Analysis (with bonds & labels) |
| **Interface** | Streamlit web UI |
| **Workflow Model** | Blind docking (unbiased, protein-wide search) |

###  Key Characteristics

| Component | Method | Details |
|-----------|--------|---------|
| **Grid Size** | Bounding box + 5Å padding | Entire protein covered, 20×20×20 Å minimum |
| **Protein Coverage** | No residue selection | All atoms included in search space |
| **Visualization Library** | py3Dmol | JavaScript-based, browser rendering, no system deps |
| **Framework** | Streamlit | Better state management, file handling, visualization integration |
| **RDKit Usage** | Ligand validation & properties | MW check, structure sanitization, coordinate extraction |
| **Meeko Usage** | PDBQT preparation | Charge assignment, atom type determination, format conversion |
| **Intermediate Files** | Temporary directory | PDBQT files, Vina config, docked output - cleaned after processing |
| **Vina Invocation** | subprocess.run() | Command-line execution with grid parameters |
| **Pose Extraction** | PDBQT parsing | Extract binding affinity, sort by energy, return top 5 |
| **Format Conversions** | PDB → PDBQT → PDB | RDKit reads SDF/MOL2, Meeko converts to PDBQT, output as PDB |
| **System Dependencies** | AutoDock Vina executable | User downloads and specifies path; Meeko as Python package |

## System Requirements

- **Operating System**: Windows, macOS, or Linux
- **Python**: 3.8 or higher
- **RAM**: Minimum 4 GB (8 GB recommended)
- **Disk Space**: ~500 MB (including dependencies)
- **Browser**: Modern browser with WebGL support (Chrome, Firefox, Safari)
- **AutoDock Vina**: 1.2.0 or higher


## Quick Start

### 1. Install System Dependencies

**Install git from [git-scm.com](https://git-scm.com/install/) & ensure git is installed**

```bash
git --version
```

**Install Python from python.org & ensure Python 3.8+ installed**

```bash
python3 --version
```

**Download AutoDock Vina: https://vina.scripps.edu/download/ & Extract**

### 2. Installation

#### Clone the repository:

```bash
git clone https://github.com/kssrikar4/pyDock.git
cd pyDock
```
##### Create virtual environment (recommended)

```bash
# Linux/macOS
python3 -m venv py
source py/bin/activate
```
or
```
# Windows
python3 -m venv py
py\Scripts\activate
```

#### Install dependencies
```
pip install -r requirements.txt
```

### 3. Verify Installation

```bash
python -c "import streamlit, rdkit, meeko; print('✓ All ready')"
```

### 4. Run Pipeline

```bash
streamlit run app.py
```

Browser opens to `http://localhost:8501`

### 5. Use the UI

1. **Step 1:** Paste AutoDock Vina path (e.g., `/usr/bin/vina`)
2. **Step 2:** Upload protein.pdb and ligand.sdf
3. **Step 3:** Click "Run Blind Docking"
4. **Step 4:** Wait 2-10 minutes for docking
5. **Step 5:** Select pose and visualization mode (Surface/Cartoon/Interaction Analysis)
6. **Step 6:** Download PDB files if needed

## File Structure

```
blind-docking-pipeline/
├── app.py                    # Main Streamlit application (full pipeline)
├── requirements.txt          # Python package dependencies
└── README.md                 # This file
```

## Architecture & Workflow

### Data Flow

```
┌─────────────────────────────────────────────────────────────┐
│ USER INPUT (via Streamlit UI)                               │
│ - Protein PDB file                                          │
│ - Ligand SDF/MOL2 file                                      │
│ - AutoDock Vina executable path                             │
└────────────────┬────────────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────────────────┐
│ AUTOMATIC GRID CALCULATION                                  │
│ - Parse PDB coordinates                                     │
│ - Calculate bounding box (min/max X,Y,Z)                    │
│ - Add 5Å padding for blind docking                          │
│ - Output: center (x,y,z) and size (sx,sy,sz)                │
└────────────────┬────────────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────────────────┐
│ LIGAND VALIDATION (RDKit)                                   │
│ - Read SDF/MOL2 format                                      │
│ - Check molecular weight (50-5000 Da)                       │
│ - Verify atom count and structure                           │
│ - Calculate properties                                      │
└────────────────┬────────────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────────────────┐
│ STRUCTURE PREPARATION (Meeko)                               │
│ - Convert protein PDB → PDBQT                               │
│ - Convert ligand SDF/MOL2 → PDBQT                           │
│ - Add partial charges (Meeko charge model)                  │
│ - Assign AutoDock atom types                                │
└────────────────┬────────────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────────────────┐
│ BLIND DOCKING (AutoDock Vina)                               │
│ - Invoke: vina --receptor ... --ligand ... --center ... ... │
│ - Search space: entire protein surface                      │
│ - Exhaustiveness: 8 (balanced accuracy/speed)               │
│ - Output: 5 models with lowest binding energies             │
└────────────────┬────────────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────────────────┐
│ POST-PROCESSING                                             │
│ - Parse Vina output PDBQT                                   │
│ - Extract binding affinity for each pose                    │
│ - Sort by affinity (most negative = best)                   │
│ - Create protein-ligand complex PDB files                   │
└────────────────┬────────────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────────────────┐
│ INTERACTIVE VISUALIZATION (py3Dmol)                         │
│ - Browser-based 3D viewer (no software installation)        │
│ - Surface Representation: Van der Waals surface             │
│ - Cartoon Representation: Protein ribbon + ligand atoms     │
│ - Interaction Analysis: H-bonds, hydrophobic, salt bridges, │
│   π-stacking, π-cation, halogen bonds with residue labels   │
│ - Controls: Rotate (click+drag), zoom (wheel), pan (right)  │
└────────────────┬────────────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────────────────┐
│ OUTPUT & DOWNLOAD                                           │
│ - Individual pose PDB files (downloadable)                  │
│ - Summary table with binding affinities                     │
└─────────────────────────────────────────────────────────────┘
```

## Grid Calculation Details

### Algorithm

```
Input: protein.pdb (all atoms)

1. Parse PDB → Extract coordinates
   atoms = [(x₁,y₁,z₁), (x₂,y₂,z₂), ..., (xₙ,yₙ,zₙ)]

2. Calculate bounding box
   bbox_min = (min(x), min(y), min(z))
   bbox_max = (max(x), max(y), max(z))

3. Add padding for blind docking
   padding = 5.0 Ų
   grid_min = bbox_min - padding
   grid_max = bbox_max + padding

4. Compute grid parameters
   center_x = (grid_min_x + grid_max_x) / 2
   size_x = grid_max_x - grid_min_x
   (repeat for y, z)

5. Enforce minimum size
   if size < 20 Ų: size = 20 Ų

Output: center (x,y,z) and size (sx,sy,sz)
```

### Why This is "Blind Docking"

- **No active site prediction:** Entire protein searched
- **No residue selection:** All atoms included
- **Unbiased:** Ligand can dock anywhere on surface
- **Complete coverage:** Bounding box + 5Å padding ensures all surface atoms sampled

## Streamlit UI

### Step 1: Configure Vina
```
┌─ AutoDock Vina Configuration ──────────────────────────────┐
│                                                            │
│ Enter path to Vina executable:                            │
│ [/usr/bin/vina                              ] ✓            │
│                                                            │
│ ✓ Vina executable found!                                  │
└────────────────────────────────────────────────────────────┘
```

### Step 2: Upload Files
```
┌─ Protein & Ligand Files ───────────────────────────────────┐
│                                                            │
│ Protein (PDB): [Choose File]                              │
│ ✓ protein.pdb uploaded (42 KB)                            │
│                                                            │
│ Ligand (SDF/MOL2): [Choose File]                          │
│ ✓ ligand.sdf uploaded (8 KB)                              │
└────────────────────────────────────────────────────────────┘
```

### Step 3: Execute Docking
```
┌─ Run Blind Docking ────────────────────────────────────────┐
│                                                            │
│ [▶️  RUN BLIND DOCKING]                                    │
│                                                            │
│ 🔬 Processing and preparing files...                      │
│ ✓ Files saved                                             │
│ ✓ Calculating grid parameters...                          │
│ Grid Center: (15.234, 20.156, 18.901)                     │
│ Grid Size: 28.5 × 32.1 × 26.3 Ų                           │
│ ✓ Validating ligand (RDKit)...                            │
│ MW: 382.41, Atoms: 31                                     │
│ ✓ Preparing protein (Meeko)...                            │
│ ✓ Preparing ligand (Meeko)...                             │
│ ✓ Running AutoDock Vina...                                │
│ ✓ Extracting top 5 poses...                               │
│ ✅ Docking completed!                                      │
└────────────────────────────────────────────────────────────┘
```

### Step 4: Visualization & Download
```
┌─ Interactive 3D Viewer ────────────────────────────────────┐
│                                                            │
│ Select Pose: [Pose 1 (ΔG = -8.23) ▼]                      │
│                                                            │
│ Representation:                                            │
│ ◉ Surface ○ Cartoon ○ Interaction Analysis                │
│                                                            │
│ ┌────────────────────────────────────────────────────────┐│
│ │  [3D rotating protein-ligand complex]                  ││
│ │  (py3Dmol viewer - interactive)                        ││
│ │  Rotate: click+drag | Zoom: scroll | Pan: right-click ││
│ └────────────────────────────────────────────────────────┘│
│                                                            │
│ Results Summary:                                           │
│ ┌────────────────────────────────────────────────────────┐│
│ │ Pose │ Binding Affinity (kcal/mol)                     ││
│ ├──────┼──────────────────────────────────────────────────││
│ │  1   │ -8.23                                            ││
│ │  2   │ -7.91                                            ││
│ │  3   │ -7.45                                            ││
│ │  4   │ -6.82                                            ││
│ │  5   │ -6.15                                            ││
│ └────────────────────────────────────────────────────────┘│
│                                                            │
│ Download: [Pose 1 PDB] [Pose 2 PDB] [All Poses]           │
└────────────────────────────────────────────────────────────┘
```

## Configuration Options

All configurable in code

```python
DEFAULT_PADDING = 5.0           # Ų beyond protein (blind docking)
VINA_EXHAUSTIVENESS = 8         # Sampling thoroughness (4-16)
NUM_POSES = 5                   # Output poses
MIN_SEARCH_SPACE = 20.0         # Minimum grid size (Ų)
```
## Troubleshooting

### "Failed to read ligand file"
- Ensure SDF/MOL2 has 3D coordinates (not SMILES)
- Regenerate from chemistry software

### "No atom coordinates in PDB"
- Verify PDB format
- Check ATOM/HETATM records present
- Download fresh PDB from RCSB/Alphafold

### Visualization not loading
- Check browser supports WebGL
- Try different browser (Chrome preferred)
- Verify internet connection (downloads 3Dmol.js)

## References & Citation

**If you use this pipeline, cite the underlying tools:**

### AutoDock Vina
```
Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed 
and accuracy of docking with a new scoring function, efficient 
optimization, and multithreading. Journal of Computational Chemistry, 
31(2), 455-461.
```

### RDKit
```
Landrum, G. RDKit: Open-source cheminformatics software. 
http://www.rdkit.org
```

### Meeko
```
Forli Lab. (2021). Meeko: Preparation of small molecules for AutoDock. 
https://github.com/forlilab/Meeko
```

### py3Dmol
```
Rego, S., & Koes, D. (2015). 3Dmol.js: molecular visualization with WebGL. 
Bioinformatics, 31(8), 1322-1324.
```
### NumPy
```
Harris, C. R., Millman, K. J., van der Walt, S. J., et al. (2020). 
Array programming with NumPy. Nature, 585(7825), 357-362.
```
### Streamlit
```
Streamlit Inc. (2024). Streamlit: A faster way to build and share data apps.
https://streamlit.io
```

## License 

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file 

### Third-Party Licenses

This software uses the following open-source packages:

| Package | License | Link |
|---------|---------|------|
| AutoDock Vina | Apache 2.0 | https://vina.scripps.edu |
| RDKit | BSD 3-Clause | http://www.rdkit.org |
| Meeko | LGPL v2.1 | https://github.com/forlilab/Meeko |
| Streamlit | Apache 2.0 | https://streamlit.io |
| NumPy | BSD 3-Clause | https://numpy.org |
| py3Dmol | BSD 3-Clause | https://3dmol.csb.pitt.edu |

All third-party licenses are respected and included in their respective packages.

## Contributing & Support

### Report Issues
Found a bug or have a feature request? Please report it:
- **Email**: kssrikar4@gmail.com
- **Include**: error messages, input files (if possible), system info

### Contributing
Contributions are welcome!

### Acknowledgments
Developed for academic research. Special thanks to the developers of AutoDock Vina, RDKit, Meeko, py3Dmol, Streamlit, and the open-source cheminformatics community for their invaluable tools.

**AI Assistance**: This pipeline was developed with assistance from Claude (Anthropic's AI assistant) for code structuring & optimization.
