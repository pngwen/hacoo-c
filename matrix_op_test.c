/* Test the matrix operations */
#include <stdio.h>
#include "matrix.h"

int main()
{
    double a[] = { 1, -3, 7, -1, 4, -7, -1, 3, -6 };
    double b[] = { -3, 3, -7, 1, 1, 0, 1, 0, 1 };
    double c[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    double d[] = { 2, 0, 0, 0, 2, 0, 0, 0, 2 };
    matrix_t *m1 = array_to_matrix(a, 3, 3);
    matrix_t *m2 = array_to_matrix(b, 3, 3);
    matrix_t *m3 = array_to_matrix(c, 3, 3);
    matrix_t *m4 = array_to_matrix(d, 3, 3);
    matrix_t *inv= new_matrix(3, 3);
    matrix_t *res= new_matrix(3, 3);

    invert_matrix(inv, m1);
    sub_matrix(res, m2, inv);
    printf("Matrix 1:\n");
    print_matrix(m1);
    printf("\nInverse: \n");
    print_matrix(inv);
    printf("\nInverse Error: \n");
    print_matrix(res);


    printf("\nMultiplication Test\n");
    mul_matrix(res, m3, m4);
    print_matrix(res);

    printf("\nTranspose mul test\n");
    mul_transpose_matrix(res, m3, m4);
    print_matrix(res);

    printf("\nAdd Matrix Test\n");
    add_matrix(res, m3, m4);
    print_matrix(res);
}
