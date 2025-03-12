import matplotlib.pyplot as plt
from Bio.PDB import PDBIO, Structure, Model, Chain, Atom, Residue
from Bio.PDB.PDBIO import Select
from python.functions_march import triangle_per_tetrahedra, diatomic_molecule_density, AtomDescriptor, polyatomic_molecule_density, PDBReader
from skimage.measure import marching_cubes
from joblib import Parallel, delayed
from tqdm import tqdm
import os

def process_block(k):
        local_triangles = []
        for j in range(1, y-1):
            for i in range(1, x-1):
                block = (i, j, k)
                for tetrahedra in cube_5_tetrahedra:
                    triangle = triangle_per_tetrahedra(tetrahedra, block, grid, threshold, "midpoint")
                    if triangle:
                        local_triangles.append(triangle)
        return local_triangles

def get_point_idx(point, points):
    for idx, p in enumerate(points):
        if p[0] == point[0] and p[1] == point[1] and p[2] == point[2]:
            return idx
    return -1

def check_incl(point, points):
    for idx, p in enumerate(points):
        if p[0] == point[0] and p[1] == point[1] and p[2] == point[2]:
            return False, idx
    return True, None

def postprocess(stream, idxs):
    added = []
    for tri in idxs:
        if len(tri) == 6:
            for i in range(3):
                if i+1 == 3:
                    added.append((tri[i], tri[0]))
                else:
                    added.append((tri[i], tri[i+1]))
            for i in range(3,6):
                if i+1 == 6:
                    added.append((tri[i], tri[3]))
                else:
                    added.append((tri[i], tri[i+1]))
        else:
            for i in range(3):
                if i+1 == 3:
                    added.append((tri[i], tri[0]))
                else:
                    added.append((tri[i], tri[i+1]))

    usable_list = []
    for (a1, a2) in added:
        if a1 > a2:
            usable_list.append((a2, a1))
        else:
            usable_list.append((a1, a2))

    set_added = set(usable_list)
    new = list(set_added)
    new.sort()

    stream.write("\n")
    for a1, a2 in new:
        stream.write(f"CONECT{a1+1:5}{a2+1:5}\n")

cube_5_tetrahedra = [
    (4, 6, 7, 8),
    (1, 5, 6, 7),
    (1, 3, 4, 7),
    (1, 2, 4, 6),
    (1, 4, 6, 7),
]

triangles = []
threshold = 0.03

if not os.path.exists("triangles.txt"):
    reader = PDBReader()
    # change with actual path
    reader.read_file('data/alanina.pdb')

    atoms = reader.atoms

    x_min = min(atom.position[0] for atom in atoms)
    x_max = max(atom.position[0] for atom in atoms)

    y_min = min(atom.position[1] for atom in atoms)
    y_max = max(atom.position[1] for atom in atoms)

    z_min = min(atom.position[2] for atom in atoms)
    z_max = max(atom.position[2] for atom in atoms)

    delta_x = 0.1
    delta_y = 0.1
    delta_z = 0.1

    print(f"x_min: {x_min}, x_max: {x_max}")
    print(f"y_min: {y_min}, y_max: {y_max}")
    print(f"z_min: {z_min}, z_max: {z_max}")

    x_range = [x_min - 5, x_max + 5]
    y_range = [y_min - 5, y_max + 5]
    z_range = [z_min - 5, z_max + 5]

    poly = polyatomic_molecule_density(atoms, x_range, y_range, z_range, 'none', delta_x, delta_y, delta_z)


    # poly = polyatomic_molecule_density(reader.atoms, x_range, y_range, z_range, 'vdw', delta_x, delta_y, delta_z)
    # poly = polyatomic_molecule_density(atom_list, x_range, y_range, z_range, 'none', delta_x, delta_y, delta_z)
    grid = poly

    # for k in range(1, z-1):
    #     for j in range(1, y-1):
    #         for i in range(1, x-1):
    #             block = (i,j,k)
    #             for tetrahedra in cube_5_tetrahedra:
    #                 triangle = triangle_per_tetrahedra(tetrahedra, block, grid, threshold, "quad")
    #                 if triangle:
    #                     triangles.append(triangle)

    x = poly.shape[0]
    y = poly.shape[1]
    z = poly.shape[2]
  
    results = Parallel(n_jobs=-1)(
        delayed(process_block)(k) 
        for k in tqdm(range(1, z-1), desc="Processing z-axis") 
    )

    for result in results:
        triangles.extend(result)

    with open("triangles.txt", "w") as f:
        for triangle in triangles:
            for point in triangle:
                f.write(f"{point[0]} {point[1]} {point[2]}\n")
            f.write("\n")

else:
    with open("triangles.txt", "r") as f:
        triangle = []
        for line in f:
            if line.strip():
                point = tuple(map(float, line.strip().split()))
                triangle.append(point)
            else:
                if triangle:
                    triangles.append(triangle)
                    triangle = []

print("# of triangles:", len(triangles))

# unique_list_of_points = []
# triangle_list = []

# for triangle in tqdm(triangles):
#     new_triangle_indexes = []
#     for point in triangle:
#         idx = get_point_idx(point, unique_list_of_points)
#         if idx != -1:
#             new_triangle_indexes.append(idx)
#         else:
#             unique_list_of_points.append(point)
#             new_triangle_indexes.append(len(unique_list_of_points)-1)
#     triangle_list.append(new_triangle_indexes)

# tri_points = {}

# idx = 0
# for triangle in tqdm(triangles, desc = "DIOCANE"): # [[25.0, 21.5, 18.5], [24.5, 21.5, 19.0], [25.0, 21.5, 19.0]]
#     triangle_idx = []
#     for number, point in enumerate(triangle): # [25.0, 21.5, 18.5]
#         point_idx = 0
#         if point not in tri_points.keys():
#             tri_points[point] = idx
#             point_idx = idx
#             idx+=1
#         else:
#             point_idx = tri_points[point]
        
#         triangle_idx.append(point_idx)
#         # print(tri_points)
#         if(idx==6):
#             print(triangle_list_idxs)
#     triangle_list_idxs.append(triangle_idx)

# unique_list_of_points = []
triangle_list = []

unique_dict_of_points = {}

# for triangle in tqdm(triangles):

# triangles = [[(1,2,3),(1,4,3),(2,2,3)],[(2,2,3),(4,0,2),(5,7,0)]]
for triangle in tqdm(triangles[0:19999]):
    new_triangle_indexes = []
    for point in triangle:
        # idx = get_point_idx(point, unique_list_of_points)

        print(point in unique_dict_of_points.keys())
        if(point in unique_dict_of_points.keys()):
            idx = unique_dict_of_points[point]
        else:
            idx = len(unique_dict_of_points)
            unique_dict_of_points[point] = idx

        new_triangle_indexes.append(idx)
        # if idx != -1:
        #     new_triangle_indexes.append(idx)
        # else:
        #     unique_list_of_points.append(point)
        #     new_triangle_indexes.append(len(unique_list_of_points)-1)
    triangle_list.append(new_triangle_indexes)


# MOLECULE
io = PDBIO()
structure = Structure.Structure("prova")
model = Model.Model(0)
chain = Chain.Chain("A")
structure.add(model)
model.add(chain)

# Add atoms
for i, coords in tqdm(enumerate(unique_dict_of_points), desc="atoms"):
    # print(i, coords)
    res = Residue.Residue((" ", i, " "), "PSE", " ")
    atom = Atom.Atom(str(i), (coords[0], coords[1], coords[2]), 1.0, 1.0, " ", str(i), i, "C")
    res.add(atom)
    chain.add(res)

io.set_structure(structure)
io.save("prova" + ".pdb")

with open("prova" + ".pdb", "r") as file:
    lines = file.readlines()

with open("prova" + ".pdb", "w") as file:
    for line in lines:
        if not (line.startswith("TER") or line.startswith("END")):
            file.write(line)

with open("prova" + ".pdb", "a") as stream:
    postprocess(stream, triangle_list)
    stream.write("END\n")

exit()

# Comment this if you want to toggle the marching_cube computation as reference
# verts, faces, normals, values = marching_cubes(grid, level=threshold)
# ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], cmap='Spectral', lw=1, alpha=0.5)

ax.view_init(elev=12., azim=295)
# ax.set_box_aspect([1,1,1])  # Aspect ratio is 1:1:1
ax.set_title("Marching Tetra Surface")

plt.show()