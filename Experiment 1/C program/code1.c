#include "stdio.h"
int main(){
    int rows;
    printf("Enter the number: ");
    scanf("%d", &rows);
    for (int i = 1; i <= rows; i++)
    {
        for (int j = 0; j < 2 * i - 1; j++){
            printf("*");
        }
        printf("\n");
    }
    return 0;
    
}
        