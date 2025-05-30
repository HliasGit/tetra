import numpy as np
import sys
import math
import matplotlib.pyplot as plt

sys.path.append('/home/elia/tesi/code/electron-density/scripts/shielding')

from skimage.measure import marching_cubes
from numpy.linalg import det
from diatomic_molecule import diatomic_molecule_density, AtomDescriptor
from pdb_reader import PDBReader
from general_molecule import polyatomic_molecule_density
from gaussian_cube_io import GaussianCubeIO
from joblib import Parallel, delayed

def midpoint(point1, point2, var): # Change val1 and val2 to val0 val1
    return (point1[var] + point2[var])/2

def linearInterpol(point1, point2, var, val0, val2, threshold): # Change val1 and val2 to val0 val1
    return point1[var] + ((point2[var]-point1[var])/(val2-val0))*(threshold-val0)

def quadraticInterpol(val0, val1, threshold, coord1, coord2, grid):
    p0 = np.array(coord1)
    p1 = np.array(coord2)
    
    p_1 = 2*p0-p1
    p2 = 2*p1-p0

    # print("p_1:", p_1)
    # print("p2:", p2)

    val_1 = grid[p_1[0], p_1[1], p_1[2]]
    val2 = grid[p2[0], p2[1], p2[2]]

    # print("val_1:", val_1)
    # print("val_2:", val2)

    S = (val0 - threshold)/(val0-val1)
    A = ((2-S)*(val_1-threshold) + (S+1)*(val2-threshold))/6
    B = A/(val0-val1)
    t = 2*S/(1+B+math.sqrt(((1+B)**2 - 4*S*B)))
    return (p0 + t*(p1-p0))

def count_points_per_tetrahedra(tetrahedra_global_coordinate, grid, threshold):
    Nm = []
    Nz = []
    Np = []
    for coord in tetrahedra_global_coordinate:
        val = (grid[coord[0], coord[1], coord[2]] - threshold)
        if val > 0:
            Np.append(coord)
        if val == 0:
            Nz.append(coord)
        if val < 0:
            Nm.append(coord)
    return Nm, Nz, Np

def check_orientation(stack, stacked_val):
    matrix = np.ones((4,4))
    for tet in range(len(stack)):
        for var in range(3):
            matrix[var+1][tet] = stack[tet][var]
    determ = det(matrix)
    if determ < 0:
        print("cose")
        count += 1
        matrix[:, [0, 1]] = matrix[:, [1, 0]]

        determ = det(matrix)

        if determ < 0:
            print("cose")

def make_triangle(stack, stacked_val, stacked_coord, pairs, int_method, grid, threshold):
    final_point = []  

    # matrix = np.ones((4,4))
    # for tet in range(len(stack)):
    #     for var in range(3):
    #         matrix[var+1][tet] = stack[tet][var]

    # determ = det(matrix)
    # if determ < 0:
    #     # stack[0], stack[1] = stack[1], stack[0]
    #     # stacked_val[0], stacked_val[1] = stacked_val[1], stacked_val[0]
    #     pairs[0][0], pairs[0][1] = pairs[0][1], pairs[0][0]
    #     pairs[1][0], pairs[1][1] = pairs[1][1], pairs[1][0]
    #     count = count +1

    for pair in pairs:
        component = []
        point1 = stack[pair[0]-1]
        point2 = stack[pair[1]-1]
        val1 = stacked_val[pair[0]-1]
        val2 = stacked_val[pair[1]-1]
        coord1 = stacked_coord[pair[0]-1]
        coord2 = stacked_coord[pair[1]-1]
        if int_method == "quad":
            component = quadraticInterpol(val1, val2, threshold, coord1, coord2, grid)
        else:
        # print("point1:", point1, "point2:", point2)
            for var in range(3):
                if int_method == "midpoint":
                    interp = midpoint(point1, point2, var)
                elif int_method == "linear":
                    interp = linearInterpol(point1, point2, var, val1, val2, threshold)
                else:
                    print("Scegline uno giusto docà")
                component.append(interp)
        final_point.append(component)
    return final_point

def get_action_value(Nm, Nz, Np):
    Nm_val = len(Nm)
    Nz_val = len(Nz)
    Np_val = len(Np)

    # if Nm_val == 0:
    #     return 0
    # if Nm_val + Nz_val == 4:
    #     return 0
    # if Nm_val == 1:
    #     if Np_val == 3:
    #         return 1
    #     if Np_val == 2:
    #         return 2
    #     if Np_val == 1:
    #         return 3
    #     if Np_val == 0:
    #         return 4
    # if Nm_val == 2:
    #     if Np_val == 2:
    #         return 7
    #     if Np_val == 1:
    #         return 5
    # if Nm_val == 3:
    #     return 6
    
    match (Nm_val, Nz_val, Np_val):
        case (0, _, _):
            return 0
        case (1, 0, 3):
            return 1
        case (1, 1, 2):
            return 2
        case (1, 2, 1):
            return 3
        case (1, 3, 0):
            return 4
        case (2, 0, 2):
            return 7
        case (2, 1, 1):
            return 5
        case (2, 2, 0):
            return 0
        case (3, 0, 1):
            return 6
        case (3, 1, 0):
            return 0
        case (4, 0, 0):
            return 0
    
def from_point_globalcube_to_local_global_point(point, block_coordinate):
    """Get the global coordinate of a vertex by specifiyng the coordinate of the
    block that contain the vertex and the vertex index you want to know

    Args:
        point (int): point index
        block_coordinate (tuple): coordinates of the block

    Raises:
        ValueError: Index must be >0, <9

    Returns:
        list: local point coordinates
        list: global point coordinates
    """
    if point < 1 or point > 8:
        raise ValueError("Number must be > 0 and < 9")

    conversion_table = [ # Note that this pattern is the 0-7 in binary from right to left
        ("i'", "j'", "k'"),
        ("i''", "j'", "k'"),
        ("i'", "j''", "k'"),
        ("i''", "j''", "k'"),
        ("i'", "j'", "k''"),
        ("i''", "j'", "k''"),
        ("i'", "j''", "k''"),
        ("i''", "j''", "k''"),
    ]

    # I want to find the 7th point of the block_coordinate = (0,0,0) cube

    one_apex = np.zeros((3), dtype = int)
    two_apex = np.zeros((3), dtype = int)
    global_coordinate = np.zeros((3), dtype = int)
    local_coordinate = np.zeros((3), dtype = int)

    for var in range(3):
        if block_coordinate[var] % 2 == 0:
            one_apex[var] = block_coordinate[var]
        else:
            one_apex[var] = block_coordinate[var]+1
        
        two_apex[var] = 2*block_coordinate[var]+1-one_apex[var]

        if conversion_table[point-1][var] == "i'" or conversion_table[point-1][var] == "j'" or conversion_table[point-1][var] == "k'":
            global_coordinate[var] = one_apex[var]
        else:
            global_coordinate[var] = two_apex[var]
        local_coordinate[var] = global_coordinate[var]-block_coordinate[var]

    return local_coordinate, global_coordinate
    
def add_point_to_stack(tetrahedra_global_coordinate, Nm, Nz, Np, Val_m, Val_z, Val_p, Coo_m, Coo_z, Coo_p, grid, threshold):
    coord = tetrahedra_global_coordinate
    # print("coord:", coord)
    val = grid[coord[0], coord[1], coord[2]]
    val_th = val-threshold
    tol = 0
    if val_th > tol:
        Np.append(coord)
        Val_p.append(val)
        Coo_p.append(coord)
    if val_th <= tol and val_th >= -tol:
        Nz.append(coord)
        Val_z.append(val)
        Coo_z.append(coord)
    if val_th < -tol:
        Nm.append(coord)
        Val_m.append(val)
        Coo_m.append(coord)
    return Nm, Nz, Np, Val_m, Val_z, Val_p, Coo_m, Coo_z, Coo_p
    

def get_grid_pts_from_actval(act_val):
    if act_val == 0:
        return
    switcher = {
        1: [[1, 2], [1, 3], [1, 4]],
        2: [[2, 2], [1, 3], [1, 4]],
        3: [[2, 2], [3, 3], [1, 4]],
        4: [[2, 2], [3, 3], [4, 4]],
        5: [[1, 4], [2, 4], [3, 3]],
        6: [[1, 4], [2, 4], [3, 4]],
        7: [[1, 4], [2, 4], [1, 3], [2, 4], [2, 3], [1, 3]]
    }
    return switcher.get(act_val, "Invalid action value")

def triangle_per_tetrahedra(tetrahedra, block, grid, threshold, int_method):
    coord = []
    stack = []
    stacked_val = []
    stacked_coord = []
    Nm = []
    Nz = []
    Np = []
    Val_m = []
    Val_z = []
    Val_p = []
    Coo_m = []
    Coo_z = []
    Coo_p = []
    for point in tetrahedra:
        _, global_coordinate = from_point_globalcube_to_local_global_point(point, block)
        Nm, Nz, Np, Val_m, Val_z, Val_p, Coo_m, Coo_z, Coo_p = add_point_to_stack(global_coordinate, Nm, Nz, Np, Val_m, Val_z, Val_p, Coo_m, Coo_z, Coo_p, grid, threshold)
        coord.append(global_coordinate)

    for i in Nm:
        stack.append(i)
    for i in Nz:
        stack.append(i)
    for i in Np:
        stack.append(i)

    for i in Val_m:
        stacked_val.append(i)
    for i in Val_z:
        stacked_val.append(i)
    for i in Val_p:
        stacked_val.append(i)

    for i in Coo_m:
        stacked_coord.append(i)
    for i in Coo_z:
        stacked_coord.append(i)
    for i in Coo_p:
        stacked_coord.append(i)

    # print(stack)
    # print(stacked_val)

    points = count_points_per_tetrahedra(coord, grid, threshold)
    # print(points)

    act_val = get_action_value(points[0], points[1], points[2])
    # print("Act_val for the tetrahedra number", act_val)

    pairs = get_grid_pts_from_actval(act_val)
    if pairs == None:
        return
    
    triangle = make_triangle(stack, stacked_val, stacked_coord, pairs, int_method, grid, threshold)

    return triangle