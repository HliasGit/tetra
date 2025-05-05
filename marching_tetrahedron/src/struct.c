#include "struct.h"

/**
 * @brief Push a new node (vertex) into the stack
 * 
 * @param start Pointer to pointer to the beginning of the stack
 * @param value Value to be added to the stack
 * @param vert Vertex (Cube vertex) to be added to the stack
 */
void push_into_stack(StackNode **start, dim_t value, CubeVertex vert, int point, int exp_pos, int *permutations){
    StackNode *new = (StackNode*) malloc(sizeof(StackNode));
    new->owned_value = value;
    new->coordinate = vert;
    new->point = point;
    new->next = NULL;

    if(*start == NULL || (*start)->owned_value > value) {
        new->next = *start;
        *start = new;
    } else {
        StackNode *current = *start;
        int real_pos = 1;
        while(current->next != NULL && current->next->owned_value <= value) {
            real_pos++;
            current = current->next;
        }
        // printf("Real Pos: %d\n", real_pos);
        // printf("Exp Pos: %d\n", exp_pos);
        if(real_pos!=exp_pos){
            // printf("INSIDE\n");
            (*permutations)++;
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

/**
 * @brief Push a triangle to the polyhedra data structure
 * 
 * Create the triangle node from the triangle. Push it into the triangle list
 * 
 * @param p Pointer to Pointer to the polyhedra
 * @param triangle Pointer to the triangle to be pushed
 * @param vertex_counter Pointer to the global vertex counter variable
 */
void push_triangle(Polyhedra **p, Triangle *triangle, size_t *vertex_counter){
    TriangleNode *new = (TriangleNode*) malloc(sizeof(TriangleNode));

    double coord1[] = {triangle->v1->x, triangle->v1->y, triangle->v1->z};
    double coord2[] = {triangle->v2->x, triangle->v2->y, triangle->v2->z};
    double coord3[] = {triangle->v3->x, triangle->v3->y, triangle->v3->z};

    new->vert1 = add(&(*p)->root_vertices, coord1, &(*vertex_counter));
    new->vert2 = add(&(*p)->root_vertices, coord2, &(*vertex_counter));
    new->vert3 = add(&(*p)->root_vertices, coord3, &(*vertex_counter));

    verbose_print("Added new triangle:\n");
    verbose_print("    Vertices: %ld, %ld, %ld\n", new->vert1, new->vert2, new->vert3);

    new->next = (*p)->triangles;
    (*p)->triangles = new;
}

/**
 * @brief Create the coordinate in a unique way and link it to the suffix tree that keeps the vertices. Returns
 *              the address of the created object
 * 
 * @param root Pointer to Pointer to the beginning of the vertices data structure
 * @param full_coordinate Pointer to the array containing the 3 coordinates of the triangle
 * @param idx Pointer to the variable of the unique idxs. Needed for the global uniqueness
 */
TriangleCoordNode* add(TriangleCoordNode **root, double *full_coordinate, size_t *idx) {
    if (*root == NULL) {
        // printf("Generating root block\n");
        TriangleCoordNode *first = malloc(sizeof(TriangleCoordNode));
        TriangleCoordNode *second = malloc(sizeof(TriangleCoordNode));
        TriangleCoordNode *third = malloc(sizeof(TriangleCoordNode));
        
        if (!first || !second || !third) {
            printf("Memory allocation failed\n");
            return false;
        }

        first->coordinate1 = full_coordinate[0];
        first->level = 1;
        first->next_list = NULL;
        first->next_level = second;
        
        second->coordinate2 = full_coordinate[1];
        second->level = 2;
        second->next_list = NULL;
        second->next_level = third;
        
        third->coordinate1 = full_coordinate[0];
        third->coordinate2 = full_coordinate[1];
        third->coordinate3 = full_coordinate[2];
        third->level = 3;
        third->next_list = NULL;
        third->next_level = NULL;
        third->index = (*idx)++;

        *root = first;
        return third;
    }
    
    TriangleCoordNode *curr = *root;
    TriangleCoordNode *prev = NULL;
    bool create = true;
    
    while (curr) {
        if (curr->coordinate1 == full_coordinate[0]) {
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

        first->coordinate1 = full_coordinate[0];
        first->level = 1;
        first->next_list = NULL;
        first->next_level = second;
        
        second->coordinate2 = full_coordinate[1];
        second->level = 2;
        second->next_list = NULL;
        second->next_level = third;
        
        third->coordinate1 = full_coordinate[0];
        third->coordinate2 = full_coordinate[1];
        third->coordinate3 = full_coordinate[2];
        third->level = 3;
        third->next_list = NULL;
        third->next_level = NULL;

        if (prev) prev->next_list = first;

        third->index = (*idx)++;

        return third;
    }
    
    create = true;
    prev = NULL;
    curr = curr->next_level;
    
    while (curr) {
        if (curr->coordinate2 == full_coordinate[1]) {
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

        second->coordinate2 = full_coordinate[1];
        second->level = 2;
        second->next_list = NULL;
        second->next_level = third;
        
        third->coordinate1 = full_coordinate[0];
        third->coordinate2 = full_coordinate[1];
        third->coordinate3 = full_coordinate[2];
        third->level = 3;
        third->next_list = NULL;
        third->next_level = NULL;

        if (prev) prev->next_list = second;
        
        third->index = (*idx)++;

        return third;
    }
    
    create = true;
    prev = NULL;
    curr = curr->next_level;
    
    while (curr) {
        if (curr->coordinate3 == full_coordinate[2]) {
            create = false;
            return curr;
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

        third->coordinate1 = full_coordinate[0];
        third->coordinate2 = full_coordinate[1];
        third->coordinate3 = full_coordinate[2];
        third->level = 3;
        third->next_level = NULL;
        third->next_list = NULL;

        if (prev) prev->next_list = third;
        
        third->index = (*idx)++;

        return third;
    }
}

/**
 * @brief free the suffix tree containing the unique vertices
 * 
 * @param TriangleCoordNode Pointer to the beginning of the tree
 */
void free_tree(TriangleCoordNode *TriangleCoordNode) {
    if (TriangleCoordNode == NULL) return;
    free_tree(TriangleCoordNode->next_level);
    free_tree(TriangleCoordNode->next_list);
    free(TriangleCoordNode);
}

/**
 * @brief free the list containing the triangles
 * 
 * @param start Pointer to the beginning of the triangles list
 */
void free_list(TriangleNode *start){
    TriangleNode *curr;
    while(start != NULL){
        curr=start;
        start = start->next;
        free(curr);
    }
}

/**
 * @brief reverse the TriangleNode list
 * 
 * Needed because that list is populated pushing in head the new data, with increasing idx.
 * I need it with the head pointing the id 0.
 * 
 * @param head Pointer to Pointer to the head of the list 
 */
void reverse_list(TriangleNode **head){
    if(head == NULL || *head == NULL){
        fprintf(stderr, "empty triangles\n");
        exit(-1);
    }

    if ((*head)->next == NULL){
        fprintf(stderr, "only one triangle\n");
        exit(-1);
    }

    TriangleNode *prev = NULL;
    TriangleNode *curr = *head;
    TriangleNode *next = NULL;

    while (curr != NULL) {
        next = curr->next;
        curr->next = prev;
        prev = curr;
        curr = next;
    }

    *head = prev;
}