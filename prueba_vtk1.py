import numpy as np
import vtk
import os
import random

print(vtk.__version__)

import bundleTools3 as bt3

def read_bundle_severalbundles(file_path):
    with open(file_path, 'rb') as f:
        bundle_names = []
        fibers = []
        
        while True:
            # Leer el nombre del bundle
            name_length = int.from_bytes(f.read(4), 'little')
            if not name_length:
                break
            name = f.read(name_length).decode('utf-8')
            bundle_names.append(name)
            
            # Leer el número de puntos en el bundle
            p = int.from_bytes(f.read(4), 'little')
            print(f"Bundle name: {name}, number of points: {p}")  # Mensaje de depuración
            
            if p < 0:
                raise ValueError(f"Invalid number of points: {p}")
            
            if p == 0:
                continue
            
            # Leer los vértices
            vertex = np.frombuffer(f.read(p * 3 * 4), 'f').reshape(-1, 3)
            fibers.append(vertex)
        
        return fibers, bundle_names

# Función para leer el archivo .hie y devolver un diccionario con los colores
def read_hie_file(hie_file):
    colors = {}
    with open(hie_file, 'r') as file:
        lines = file.readlines()
        name = None
        for line in lines:
            line = line.strip()
            if line.startswith("name"):
                name = line.split(" ")[1]
            elif line.startswith("color") and name is not None:
                color = list(map(float, line.split(" ")[1:]))
                colors[name] = color
                name = None
    return colors

# Función para convertir fibras a vtkPolyData con colores
def fibers_to_vtk_polydata(fibers, colors, bundle_names):
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    color_array = vtk.vtkUnsignedCharArray()
    color_array.SetNumberOfComponents(3)
    color_array.SetName("Colors")

    point_id = 0
    for bundle_idx, bundle in enumerate(fibers):
        color = colors.get(bundle_names[bundle_idx], [255, 255, 255])  # Blanco por defecto
        for fiber in bundle:
            line = vtk.vtkPolyLine()
            line.GetPointIds().SetNumberOfIds(len(fiber))
            for i, point in enumerate(fiber):
                points.InsertNextPoint(point)
                line.GetPointIds().SetId(i, point_id)
                point_id += 1
            lines.InsertNextCell(line)
            for _ in fiber:
                color_array.InsertNextTuple3(color[0], color[1], color[2])

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    polydata.GetCellData().SetScalars(color_array)
    return polydata

# Función para visualizar vtkPolyData
def visualize_polydata(polydata):
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0.0, 0.0, 0.0)  # Fondo negro

    render_window = vtk.vtkRenderWindow()
    render_window.SetAlphaBitPlanes(1)  # Para soporte de transparencia
    render_window.AddRenderer(renderer)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    render_window.Render()
    interactor.Start()

path = os.getcwd()

bundle_file = os.path.join(path, "visualization_test", "100206", "centroids.bundles")
hie_file = os.path.join(path, "visualization_test", "QB_21.hie")

# Leer fibras y nombres de los bundles
fibers, bundle_names = read_bundle_severalbundles(bundle_file)

print("\n\n",len(fibers))

# Leer colores desde el archivo .hie
colors = read_hie_file(hie_file)

# Convertir las fibras a vtkPolyData con colores
polydata = fibers_to_vtk_polydata(fibers, colors, bundle_names)

# Visualizar las fibras
visualize_polydata(polydata)
