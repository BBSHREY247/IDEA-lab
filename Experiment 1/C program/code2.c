#include <stdio.h>
#include <stdlib.h>
int main(int argc, char const *argv[])
{
    char name[50];
    printf("Enter Your Name: ");
    fgets(name, sizeof(name), stdin);
    printf("Hello %s", name);
    return 0;
}
