#include "triangles.h"

/**
 * @brief Returns the pointer to the triangle generated from the interpolation
 * 
 * @param stack Pointer to the stack beginning
 * @param pairs Pointer to the pairs in order to access the stack
 * @param two_triangles Boolean to handle the case of building two triangles
 * @param threshold Threshold of the isosurface
 * @param func_ptr Function pointer to invoke dynamically the function
 */
Triangle *make_triangle(StackNode *stack, int *pairs, bool two_triangles, dim_t threshold, void (*func_ptr)(TriangleVertex*, CubeVertex*, CubeVertex*, dim_t*, dim_t*, dim_t)){

    Triangle *triangle = (Triangle *)malloc(sizeof(Triangle));

    triangle->v1 = (TriangleVertex *)malloc(sizeof(TriangleVertex));
    triangle->v2 = (TriangleVertex *)malloc(sizeof(TriangleVertex));
    triangle->v3 = (TriangleVertex *)malloc(sizeof(TriangleVertex));

    for(int idx=0; idx<6; idx+=2){
        CubeVertex *point1 = get_coordinate_by_idx(stack, pairs[idx]-1);
        CubeVertex *point2 = get_coordinate_by_idx(stack, pairs[idx+1]-1);
        dim_t *val1 = get_value_by_idx(stack, pairs[idx]-1);
        dim_t *val2 = get_value_by_idx(stack, pairs[idx+1]-1);

        if(idx==0){
            func_ptr(triangle->v1, point1, point2, val1, val2, threshold);
        }        
        if(idx==2){
            func_ptr(triangle->v2, point1, point2, val1, val2, threshold);
        }
        if(idx==4){
            func_ptr(triangle->v3, point1, point2, val1, val2, threshold);
        }
    }

    if (two_triangles)
    {
        for(int idx=6; idx<12; idx+=2){
            CubeVertex *point1 = get_coordinate_by_idx(stack, pairs[idx]-1);
            CubeVertex *point2 = get_coordinate_by_idx(stack, pairs[idx+1]-1);
            dim_t *val1 = get_value_by_idx(stack, pairs[idx]-1);
            dim_t *val2 = get_value_by_idx(stack, pairs[idx+1]-1);
    
            if(idx==6){
                func_ptr(triangle->v1, point1, point2, val1, val2, threshold);
            }
            if(idx==8){
                func_ptr(triangle->v2, point1, point2, val1, val2, threshold);
            }
            if(idx==10){
                func_ptr(triangle->v3, point1, point2, val1, val2, threshold);
            }
        }
    }
    return triangle;
}

/**
 * @brief Compute the midpoint interpolation
 * 
 * @param vtx Pointer to the vertex to be computed, to store it
 * @param point1 Pointer to the first CubeVertex
 * @param point1 Pointer to the second CubeVertex
 * @param val1 Pointer to the grid value of the first point. Needed just for the function pointer
 * @param val2 Pointer to the grid value of the second point.Needed just for the function pointer
 * @param threshold Value of the threshold value. Needed just for the function pointer
 */
void midpoint_interpol(TriangleVertex *vtx, CubeVertex *point1, CubeVertex *point2, dim_t *val1, dim_t *val2, dim_t threshold){
    vtx->x = (coord_t)(point2->x+point1->x)/2;
    vtx->y = (coord_t)(point2->y+point1->y)/2;
    vtx->z = (coord_t)(point2->z+point1->z)/2;
}

/**
 * @brief Compute the linear interpolation
 * 
 * @param vtx Pointer to the vertex to be computed, to store it
 * @param point1 Pointer to the first CubeVertex
 * @param point1 Pointer to the second CubeVertex
 * @param val1 Pointer to the grid value of the first point
 * @param val2 Pointer to the grid value of the second point
 * @param threshold Value of the threshold value 
 */
void linear_interpol(TriangleVertex *vtx,CubeVertex *point1, CubeVertex *point2, dim_t *val1, dim_t *val2, dim_t threshold){
    
    if ((*val2) - (*val1) == 0) {
        fprintf(stderr, "Error: Division by zero in linear interpolation\n");
        exit(-1);
    }

    vtx->x = ((coord_t)(point1->x) + ((((coord_t)point2->x - (coord_t)point1->x) / ((*val2) - (*val1))) * (threshold - (*val1))));
    vtx->y = ((coord_t)(point1->y) + ((((coord_t)point2->y - (coord_t)point1->y) / ((*val2) - (*val1))) * (threshold - (*val1))));
    vtx->z = ((coord_t)(point1->z) + ((((coord_t)point2->z - (coord_t)point1->z) / ((*val2) - (*val1))) * (threshold - (*val1))));
}