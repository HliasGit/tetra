import numpy as np
import matplotlib.pyplot as plt

def count_points_per_tetrahedra(tetrahedra_global_coordinate):
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


def make_triangle(stack, pairs):
    final_point = []    
    for pair in pairs:
        component = []
        point1 = stack[pair[0]-1]
        point2 = stack[pair[1]-1]
        for var in range(3):
            mid = (point1[var] + point2[var])/2
            component.append(mid)
        final_point.append(component)
    return final_point

def get_action_value(Nm, Nz, Np):
    Nm_val = len(Nm)
    Nz_val = len(Nz)
    Np_val = len(Np)
    if Nm_val == 0:
        return 0
    if Nm_val + Nz_val == 4:
        return 0
    if Nm_val == 1:
        if Np_val == 3:
            return 1
        if Np_val == 2:
            return 2
        if Np_val == 1:
            return 3
        if Np_val == 0:
            return 4
    if Nm_val == 2:
        if Np_val == 2:
            return 7
        if Np_val == 1:
            return 5
    if Nm_val == 3:
        return 6
    
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
    
def add_point_to_stack(tetrahedra_global_coordinate, Nm, Nz, Np):
    coord = tetrahedra_global_coordinate
    val = (grid[coord[0], coord[1], coord[2]] - threshold)
    if val > 0:
        Np.append(coord)
    if val == 0:
        Nz.append(coord)
    if val < 0:
        Nm.append(coord)
    return Nm, Nz, Np
    

def get_grid_pts_from_actval(act_val):
    if act_val == 0:
        return
    switcher = {
        1: [(1,2), (1,3),(1,4)],
        2: [(2,2), (1,3),(1,4)],
        3: [(2,2), (3,3),(1,4)],
        4: [(2,2), (3,3),(4,4)],
        5: [(1,4), (2,4),(3,3)],
        6: [(1,4), (2,4),(3,4)],
        7: [(1,4), (2,4),(1,3), (2,4), (2,3),(1,3)]
    }
    return switcher.get(act_val, "Invalid action value")

cube_5_tetrahedra = [
    (4, 6, 7, 8),
    (1, 5, 6, 7),
    (1, 3, 4, 7),
    (1, 2, 4, 6),
    (1, 4, 6, 7),
]

n = 4
c = n+0.5
threshold = 3.5

grid = np.ones((n+1, n+1, n+1), dtype=float)

for k in range (0, n+1):
    for j in range (0, n+1):
        for i in range (0, n+1):
            if k == n-1:
                grid[i,j,k] = n
            elif k == n:
                grid[i,j,k] = n+1

grid[1,1,1] = n+1

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
spacing = 0.1

# for idx, tetrahedra in enumerate(cube_5_tetrahedra):

triangles = []

x = 3
y = 3
z = 3

for k in range(z):
    for j in range(y):
        for i in range(x):
            print(i,j,k)
            block = (i,j,k)
            for tetrahedra in cube_5_tetrahedra:
                coord = []
                stack = []
                Nm = []
                Nz = []
                Np = []
                for point in tetrahedra:
                    _, global_coordinate = from_point_globalcube_to_local_global_point(point, block)
                    Nm, Nz, Np = add_point_to_stack(global_coordinate, Nm, Nz, Np)
                    coord.append(global_coordinate)

                for i in Nm:
                    stack.append(i)
                for i in Nz:
                    stack.append(i)
                for i in Np:
                    stack.append(i)

                # print(stack)

                points = count_points_per_tetrahedra(coord)
                # print(points)

                act_val = get_action_value(points[0], points[1], points[2])
                # print("Act_val for the tetrahedra number", idx, ", ", act_val)

                pairs = get_grid_pts_from_actval(act_val)
                # print(pairs)
                if pairs == None:
                    continue

                triangle = make_triangle(stack, pairs)
                # print(triangle)
                triangles.append(triangle)


# Plot the grid points with color based on the value of the grid
for i in range(x+1):
    for j in range(y+1):
        for k in range(z+1):
            ax.scatter(i, j, k, c=plt.cm.viridis(grid[i, j, k] / grid.max()), marker='o')

ax.view_init(elev=12, azim=270)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_title("Grid plot with function value in colours")

for triangle in triangles:
    for i in range(len(triangle)):
        point1 = triangle[i]
        point2 = triangle[(i + 1) % len(triangle)]
        ax.plot([point1[0], point2[0]], [point1[1], point2[1]], [point1[2], point2[2]], color="black")


ax.view_init(elev=12., azim=295)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_title("Grid plot with function value in colours")

plt.show()