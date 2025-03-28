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

void push_triangle(Polyhedra **p, Triangle *triangle, size_t *vertex_counter){
    TriangleNode *new = (TriangleNode*) malloc(sizeof(TriangleNode));

    // new->vert1 = push_vertex(&(*p)->vertices, triangle->v1, &(*vertex_counter));
    // new->vert2 = push_vertex(&(*p)->vertices, triangle->v2, &(*vertex_counter));
    // new->vert3 = push_vertex(&(*p)->vertices, triangle->v3, &(*vertex_counter));

    double coord1[] = {triangle->v1->x, triangle->v1->y, triangle->v1->z};
    double coord2[] = {triangle->v2->x, triangle->v2->y, triangle->v2->z};
    double coord3[] = {triangle->v3->x, triangle->v3->y, triangle->v3->z};

    new->vert1 = add(&(*p)->root, coord1, &(*vertex_counter));
    new->vert2 = add(&(*p)->root, coord2, &(*vertex_counter));
    new->vert3 = add(&(*p)->root, coord3, &(*vertex_counter));

    verbose_print("Added new triangle:\n");
    verbose_print("    Vertices: %ld, %ld, %ld\n", new->vert1, new->vert2, new->vert3);

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

int add(TriangleCoordNode **root, double *full_coordinate, size_t *idx) {
    if (*root == NULL) {
        // printf("Generating root block\n");
        TriangleCoordNode *first = malloc(sizeof(TriangleCoordNode));
        TriangleCoordNode *second = malloc(sizeof(TriangleCoordNode));
        TriangleCoordNode *third = malloc(sizeof(TriangleCoordNode));
        
        if (!first || !second || !third) {
            printf("Memory allocation failed\n");
            return false;
        }

        first->coordinate = full_coordinate[0];
        first->level = 1;
        first->next_list = NULL;
        first->next_level = second;
        
        second->coordinate = full_coordinate[1];
        second->level = 2;
        second->next_list = NULL;
        second->next_level = third;
        
        third->coordinate = full_coordinate[2];
        third->level = 3;
        third->next_list = NULL;
        third->next_level = NULL;
        third->index = (*idx)++;

        // printf("Block root coordinates: %d -> %d -> %d\n", first->coordinate, second->coordinate, third->coordinate);
        *root = first;
        return third->index;
    }
    
    TriangleCoordNode *curr = *root;
    TriangleCoordNode *prev = NULL;
    bool create = true;
    
    while (curr) {
        if (curr->coordinate == full_coordinate[0]) {
            create = false;
            break;
        }
        prev = curr;
        curr = curr->next_list;
    }
    
    if (create) {
        // printf("Generating complete block\n");
        TriangleCoordNode *first = malloc(sizeof(TriangleCoordNode));
        TriangleCoordNode *second = malloc(sizeof(TriangleCoordNode));
        TriangleCoordNode *third = malloc(sizeof(TriangleCoordNode));
        
        if (!first || !second || !third) {
            printf("Memory allocation failed\n");
            exit(-1);
        }

        first->coordinate = full_coordinate[0];
        first->level = 1;
        first->next_list = NULL;
        first->next_level = second;
        
        second->coordinate = full_coordinate[1];
        second->level = 2;
        second->next_list = NULL;
        second->next_level = third;
        
        third->coordinate = full_coordinate[2];
        third->level = 3;
        third->next_list = NULL;
        third->next_level = NULL;

        // printf("Block coordinates: %d -> %d -> %d\n", first->coordinate, second->coordinate, third->coordinate);
        if (prev) prev->next_list = first;

        third->index = (*idx)++;

        return third->index;
    }
    
    create = true;
    prev = NULL;
    curr = curr->next_level;
    
    while (curr) {
        if (curr->coordinate == full_coordinate[1]) {
            create = false;
            break;
        }
        prev = curr;
        curr = curr->next_list;
    }
    
    if (create) {
        // printf("Generating two-block\n");
        TriangleCoordNode *second = malloc(sizeof(TriangleCoordNode));
        TriangleCoordNode *third = malloc(sizeof(TriangleCoordNode));
        
        if (!second || !third) {
            printf("Memory allocation failed\n");
            exit(-1);
        }

        second->coordinate = full_coordinate[1];
        second->level = 2;
        second->next_list = NULL;
        second->next_level = third;
        
        third->coordinate = full_coordinate[2];
        third->level = 3;
        third->next_list = NULL;
        third->next_level = NULL;

        // printf("Block coordinates: %d -> %d\n", second->coordinate, third->coordinate);
        if (prev) prev->next_list = second;
        
        third->index = (*idx)++;

        return third->index;
    }
    
    create = true;
    prev = NULL;
    curr = curr->next_level;
    
    while (curr) {
        if (curr->coordinate == full_coordinate[2]) {
            create = false;
            return curr->index;
        }
        prev = curr;
        curr = curr->next_list;
    }
    
    if (create) {
        // printf("Generating one-block\n");
        TriangleCoordNode *third = malloc(sizeof(TriangleCoordNode));
        
        if (!third) {
            printf("Memory allocation failed\n");
            exit(-1);
        }

        third->coordinate = full_coordinate[2];
        third->level = 3;
        third->next_level = NULL;
        third->next_list = NULL;

        // printf("Block coordinate: %d\n", third->coordinate);
        if (prev) prev->next_list = third;
        
        third->index = (*idx)++;

        return third->index;
    }
}

void print_vertices(TriangleCoordNode *TriangleCoordNode, double *first, double *second, FILE *fptr) {
    if (TriangleCoordNode == NULL) return;
    
    while (TriangleCoordNode) {
        if(TriangleCoordNode->level == 1){
            *first = TriangleCoordNode->coordinate;
        }
        if(TriangleCoordNode->level == 2){
            *second = TriangleCoordNode->coordinate;
        }
        if(TriangleCoordNode->level == 3){
            // printf("Level: %d, Coordinate: %f, %f, %f\n", TriangleCoordNode->level, *first, *second, TriangleCoordNode->coordinate);
            // printf("%5ld", TriangleCoordNode->index);
            // exit(-1);
            fprintf( fptr, "ATOM %6ld 0    PSE A   0      %6.3f  %6.3f  %6.3f  1.00  1.00           C\n", 
                    TriangleCoordNode->index, *first, *second, TriangleCoordNode->coordinate);
        }
        print_vertices(TriangleCoordNode->next_level, first, second, fptr);
        TriangleCoordNode = TriangleCoordNode->next_list;
    }
}

void free_tree(TriangleCoordNode *TriangleCoordNode) {
    if (TriangleCoordNode == NULL) return;
    free_tree(TriangleCoordNode->next_level);
    free_tree(TriangleCoordNode->next_list);
    free(TriangleCoordNode);
}

void free_list(TriangleNode *start){
    TriangleNode *curr;
    while(start != NULL){
        curr=start;
        start = start->next;
        free(curr);
    }
}