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
