import vtk
import ase.io as ase_io
import numpy as np
import json
from pymatgen.analysis.local_env import EconNN
import sys
from pymatgen.core import Lattice, Structure, Molecule



ecnn=EconNN()
#structure = Structure.from_file("./POSCAR3D/POSCAR")
# Load the POSCAR file using ASE
atoms = ase_io.read('./POSCAR3D/POSCAR')
structure=Structure(
    lattice=atoms.cell,
    species=atoms.get_chemical_symbols(),
    coords=atoms.get_scaled_positions()
)

# Loop over all atoms in the structure to add them to the renderer
positions = atoms.get_positions()
symbols = atoms.get_chemical_symbols()
lattice=atoms.cell
unique_atoms=set(symbols)
with open('atomic_data.json','r') as f:
    atomic_table=json.load(f)

atoms_colors={ele:[item/255 for item in atomic_table[ele]['color'][:3]] for ele in unique_atoms }

# Create a VTK renderer, render window, and interactor
renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

# Function to create a sphere source for each atom
def create_atom_sphere(position, radius=0.5):
    sphere_source = vtk.vtkSphereSource()
    sphere_source.SetCenter(position)
    sphere_source.SetRadius(radius)
    sphere_source.Update()
    return sphere_source

# Function to create a mapper and actor from the sphere source
def create_atom_actor(sphere_source, color=(1, 1, 1)):
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphere_source.GetOutputPort())
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)  # Set the color of the atom
    return actor


# Function to draw the unit cell using lines
def draw_unit_cell(lattice):
    # Define the corners of the unit cell
    corners = np.array([
        [0, 0, 0],
        lattice[0],
        lattice[1],
        lattice[2],
        lattice[0] + lattice[1],
        lattice[1] + lattice[2],
        lattice[0] + lattice[2],
        lattice[0] + lattice[1] + lattice[2]
    ])
    
    # Define the pairs of corners that form the edges of the unit cell
    edges = [
        (0, 1), (0, 2), (0, 3),
        (1, 4), (1, 6), (2, 4),
        (2, 5), (3, 5), (3, 6),
        (4, 7), (5, 7), (6, 7)
    ]
    
    # Create a polyline for the unit cell
    points = vtk.vtkPoints()
    for corner in corners:
        points.InsertNextPoint(corner)
    
    lines = vtk.vtkCellArray()
    for edge in edges:
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, edge[0])
        line.GetPointIds().SetId(1, edge[1])
        lines.InsertNextCell(line)
    
    poly_data = vtk.vtkPolyData()
    poly_data.SetPoints(points)
    poly_data.SetLines(lines)
    
    # Create a mapper and actor for the unit cell
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(poly_data)
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0, 0, 0)  
    actor.GetProperty().SetLineWidth(2)  
    
    return actor

# Function to add a legend to the corner
def add_legend(renderer):
    color=atoms_colors[element]
    
    # Create a vtkLegendBoxActor
    legend = vtk.vtkLegendBoxActor()
    legend.SetNumberOfEntries(len(atoms_colors))
    
    for i, (atom_type, color) in enumerate(atoms_colors.items()):
        # Create a small sphere for the legend
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(0.2)
        sphere.Update()
        
        # Create a mapper and actor for the legend sphere
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
        
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)
        
        # Add the actor to the legend
        legend.SetEntry(i, actor.GetMapper().GetInput(), atom_type, color)
    
    # Set legend properties
    legend.SetPosition(0.75, 0.05)  # Bottom right corner (can adjust as needed)
    legend.SetWidth(0.2)
    legend.SetHeight(0.15)
    legend.UseBackgroundOn()
    legend.SetBackgroundColor(1,1,1)  
    legend.GetEntryTextProperty().SetColor(1, 1, 1) 
    
    # Add the legend to the renderer
    renderer.AddActor(legend)

# Get atom positions and chemical symbols
positions = atoms.get_positions()
symbols = atoms.get_chemical_symbols()

# Add atoms as spheres
for i, atom in enumerate(atoms):
    position = positions[i]
    element = symbols[i]
    color=atoms_colors[element]
    radius=atomic_table[element]['atomic_radius']
    sphere = create_atom_sphere(position)
    actor = create_atom_actor(sphere, color)
    renderer.AddActor(actor)


# Draw the unit cell
unit_cell_actor = draw_unit_cell(lattice)
renderer.AddActor(unit_cell_actor)

# Add the legend to the corner
add_legend(renderer)

# Set the background color and initialize the renderer
renderer.SetBackground(1, 1, 1)  
renderWindow.SetSize(800, 800)

# Start the visualization
renderWindow.Render()
renderWindowInteractor.Start()
