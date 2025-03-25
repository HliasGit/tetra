#include "struct.h"

/**
 * @brief Push a new node (vertex) into the stack
 * 
 * @param start Pointer to pointer to the beginning of the stack
 * @param value Value to be added to the stack
 * @param vert Vertex (Cube vertex) to be added to the stack
 */
void push_into_stack(StackNode **start, dim_t value, CubeVertex vert){
    StackNode *new = (StackNode*) malloc(sizeof(StackNode));
    new->owned_value = value;
    new->coordinate = vert;
    new->next = NULL;

    if(*start == NULL || (*start)->owned_value >= value) {
        new->next = *start;
        *start = new;
    } else {
        StackNode *current = *start;
        while(current->next != NULL && current->next->owned_value < value) {
            current = current->next;
        }
        new->next = current->next;
        current->next = new;
    }
}

/**
 * @brief Free the stack that stores the tetrahedra
 * 
 * @param start Pointer to pointer to the beginning of the stack
 */
void free_stack(StackNode **start){
    StackNode *current = *start;
    StackNode *next;

    while(current != NULL){
        next = current->next;
        free(current);
        current = next;
    }

    *start = NULL;
}

/**
 * @brief Access the stack to get the coordinates
 * 
 * @param start Pointer to the stack beginning
 * @param idx Index needed to be accessed
 */
CubeVertex *get_coordinate_by_idx(StackNode *start, int idx){
    int i=0;
    StackNode *ptr = start;
    while(i<idx){
        i++;
        if(ptr == NULL){
            fprintf(stderr, "Stack is smaller than idx\n");
            exit(-1);
        }
        verbose_print("Iter for the coordinates\n");
        ptr = ptr->next;
    }

    verbose_print("Found\n");
    return &ptr->coordinate;
}

/**
 * @brief Access the stack to get the value
 * 
 * @param start Pointer to the stack beginning
 * @param idx Index needed to be accessed
 */
dim_t *get_value_by_idx(StackNode *start, int idx){
    int i=0;
    StackNode *ptr = start;
    while(i<idx){
        i++;
        if(ptr == NULL){
            fprintf(stderr, "Stack is smaller than idx\n");
            exit(-1);
        }
        verbose_print("Iter for the coordinates\n");
        ptr = ptr->next;
    }

    verbose_print("Found\n");
    return &ptr->owned_value;
}

int push_vertex(VertexNode **start, TriangleVertex *m_vertex, int *vertex_counter){
    VertexNode *new = (VertexNode*) malloc(sizeof(VertexNode));
    new->vertex = m_vertex;
    new->next = NULL;

    if(*start == NULL || coordinate_less_than(m_vertex, (*start)->vertex)){
        (*vertex_counter)++;
        new->idx = *vertex_counter;
        new->next = *start;
        *start = new;
        return new->idx;
    } else {
        VertexNode *current = *start;
        while(current->next != NULL && coordinate_less_than(current->next->vertex, m_vertex)){
            current = current->next;
        }
        if(current->next != NULL){
            current = current->next;
        }
        if(!coordinate_equals(current->vertex, m_vertex)){
            (*vertex_counter)++;
            new->idx = *vertex_counter;
            new->next = current->next;
            current->next = new;
            return new->idx;
        } else {
            return current->idx;
        }
    }
}

void print_vertex_list(VertexNode *start){
    int count = 0;
    while(start != NULL){
        printf("    Vertex %d:\n", start->idx);
        printf("        x: %f\n", start->vertex->x);
        printf("        y: %f\n", start->vertex->y);
        printf("        z: %f\n", start->vertex->z);
        
        start = start->next;
    }
}

void push_triangle(Polyhedra **p, Triangle *triangle, int *vertex_counter){
    TriangleNode *new = (TriangleNode*) malloc(sizeof(TriangleNode));
    new->vert1 = push_vertex(&(*p)->vertices, triangle->v1, &(*vertex_counter));
    new->vert2 = push_vertex(&(*p)->vertices, triangle->v2, &(*vertex_counter));
    new->vert3 = push_vertex(&(*p)->vertices, triangle->v3, &(*vertex_counter));

    // printf("Added new triangle:\n");
    // printf("    Vertices: %d, %d, %d\n", new->vert1, new->vert2, new->vert3);

    new->next = (*p)->triangles;
    (*p)->triangles = new;
}

void print_triangle_list(TriangleNode *start){
    int count = 0;
    while(start != NULL){
        count++;
        printf("Triangle %d\n", count);
        printf("    Vertex:\n");
        printf("        x: %d\n", start->vert1);
        printf("        y: %d\n", start->vert2);
        printf("        z: %d\n", start->vert3);
        
        start = start->next;
    }
}