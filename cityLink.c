/** @file cityLink.c
 *  @brief The program find the transitional closure of a graph,
 *  and a path between two cities.
 * 
 *  The program reads a graph from an input file, processes the data,
 *  and performs operations based on user-defined arguments. Based on that,
 *  creates a neighbor table from the input file, generates the R* table representing
 *  edges of the graph, computes the transitive closure of the graph using the R* table, 
 *  and find paths between specified source and destination cities.
 * 
 *  @author Elena Eleftheriou
 *  @bug No know bugs.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <unistd.h>
#include <stdbool.h>

/**
 * @brief Reads the input file.
 * 
 *
 * This function reads the input file containing the graph information,
 * processes the data, and performs operations based on the input arguments provided.
 * More specifically, it creates the neighbor, and prints it using printNeighbor, and
 * then with prnt checks if the user wants to create R and find the transitional closure or the path
 * between two cities.
 *
 * @param filename input file.
 * @param prnt Array indicating which information to print.
 * @param sourceCity Source city
 * @param destCity Destination city
 * @return void
 */
void readFile(char *, bool *, int, int);

/**
 * @brief Prints the R* table to a file.
 *
 * This function prints the R* table to the out-filename, where the filename is the 
 * given input filename.
 *
 * @param filename Name of the output file.
 * @param R R* table.
 * @param size Size of the R* table.
 * @return void
 */
void printFile(char *, int **, int);

/**
 * @brief Reads command line arguments and processes them.
 *
 * This function parses the command line arguments and performs
 * necessary operations based on the arguments provided, using getopt.
 * The valid arguments are i, o, p, r(with i being necessary) and are used as descriped.
 * If the user gave some other character or not i, or something else goes wrong, 
 * the program stops running.
 * Also, it has a dynamic array, prnt, that is used for us to know what the program
 * has to print and where.
 * Lastly, it calls readFile() to continue the program. * 
 *
 * @param argumentc Number of command line arguments.
 * @param argumentv Array of command line argument strings.
 * @return void
 */
void readArguments(int, char **);

/**
 * @brief Prints the neighbor table.
 *
 * This function prints the neighbor table.
 *
 * @param G G
 * @param N size of G
 * @return void
 */
void printNeighbor(int **, int);

/**
 * @brief Prints the R* table.
 *
 * This function prints the R* table.
 *
 * @param R R* table.
 * @param sizeR Size of the R* table.
 * @return void
 */
void printR(int **, int);

/**
 * @brief Creates the R* table based on the given graph.
 *
 * This function creates the R* table, which represents the edges of the graph.
 * It gets from the neighbor table the values that are 1 and puts them in R*, while
 * creating its columns. * 
 *
 * @param R R* table.
 * @param G G
 * @param N size of G
 * @return Size of the R* table.
 */
int createR(int ***, int **, int *);

/**
 * @brief Computes the transitive closure of the graph.
 *
 * This function computes the transitive closure of the graph represented by the R* table.
 * More specifically, it uses the given pseudocode to find the connections and then,
 * based on prnt it prints in a new file or in the terminal the new R, or find the path between
 * two cities.
 *
 * @param R R* table.
 * @param N size of R*
 * @param prnt Array indicating which information to print.
 * @param sourceCity Source city
 * @param destCity Destination city
 * @param filename filename
 * @return void
 */
void transitionalClosure(int **, int, bool *, int, int, char *);

/**
 * @brief Sets a new entry in the R* table.
 *
 * This function adds a new edge to the R* table.
 * More specifically, it creates a new dynamic array with the length+1 of the given one
 * so that we can add another row to it.
 * Also, if the flag is 1 then it adds the given u and w to the new array and then returns it.
 *
 * @param R R* table.
 * @param N size of R*
 * @param flag Flag indicating whether to add a new edge.
 * @param u Source city
 * @param w Destination city
 * @return Pointer to the new R* table.
 */
int ** setR(int **, int, int, int, int);

/**
 * @brief Checks if the values are already in R*
 * 
 * Checks if u and w is already in the array R*.
 * 
 * @param newR newR*table
 * @param u u
 * @param w w
 * @param N size of newR*
 * @return true if there are not the same values
*/
bool checkW(int **, int, int, int);

/**
 * @brief Finds the path between two cities.
 *
 * This function finds the path between two cities using the R* table.
 * The function is reccursive and uses findA() to find the index of the array
 * that has the destination city in. Then it calls itself with j as the source city of the destination city
 * that findA returned.
 *
 * @param R R* table.
 * @param size Size of the R* table.
 * @param i Source city.
 * @param j Destination city.
 * @param visited Array to keep track of visited cities.
 * @return void
 */
bool findLink(int **, int, int , int, bool *, int);

/**
 * @brief Checks if a link exists between two cities.
 *
 * This function checks if a direct link exists between two cities in the newR* table
 * that was made with transitionalClosure()
 *
 * @param R newR* table.
 * @param size Size of the R* table.
 * @param i Source city.
 * @param j Destination city.
 * @return True if a link exists.
 */
bool checkLink(int **, int, int, int);

/** @brief Program entrypoint.
 *  
 *  It calls the function readArguments() with the arguments argc and argv
 *  to start reading what the user wants to do with this program from the terminal.
 *
 * @return 0
 */
int main(int argc, char *argv[]){
    readArguments(argc, argv);
    return 0;
}


void readArguments(int argumentc, char **argumentv){
    char *filename = NULL;
    bool *prnt = (bool *)malloc(4*sizeof(bool));
    if (prnt == NULL) { // Check if memory allocation was successful
        perror( "Memory allocation failed.\n");
        exit(-1);
    }
    int sourceCity = 0, destCity = 0;
    char *filenameOut = NULL; 
    int index;
    int c;
    bool flag= false;
    opterr = 0;
    if(argumentc == 1){
        printf("No command line arguments given!\nUsage: <executable> -i <inputfile> [-r <source >,<destination> -p -o]\n");
        exit(-1);
    }   
    while ((c = getopt (argumentc, argumentv, "i:por:")) != -1) {
        switch (c){
            case 'r':
                if (sscanf(optarg, "%d,%d", &sourceCity, &destCity) != 2) {
                    printf("Invalid -r format\nUsage: <executable> -i <inputfile> [-r <source >,<destination> -p -o]\n");
                    exit(-1);
                }
                prnt[0] = 1;
                prnt[1] = 1; //print call for r
                
            break;
            case 'p':
                prnt[0] = 1;
                prnt[2] = 1; //print call for p
            break;
            case 'o':
                prnt[0] = 1;
                prnt[3] = 1; //print call for o
            break;
            case 'i':
                if(optarg!=NULL && strstr(optarg, ".txt") != NULL){
                    filename = (char *)malloc(strlen(optarg)*sizeof(char));
                    if (filename == NULL) { // Check if memory allocation was successful
                        perror( "Memory allocation failed.\n");
                        exit(-1);
                    }
                    strcpy(filename, optarg);
                    flag = 1;
                }
                else{
                    printf("No input file given!\nUsage: <executable> -i <inputfile> [-r <source >,<destination> -p -o]\n");
                    exit(-1);
                }
            break;
            case '?':
                if(optopt == 'i')
                    printf("./cityLink: option requires an argument -- 'i'\nUsage: <executable> -i <inputfile> [-r <source >,<destination> -p -o]\n");
                else if (isprint (optopt))
                    fprintf (stderr, "./cityLink: invalid option -- '%c'\nUsage: <executable> -i <inputfile> [-r <source >,<destination> -p -o]\n", optopt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                exit(-1);
        }
    }
    if(!flag){
        printf("No input file given!\nUsage: <executable> -i <inputfile> [-r <source >,<destination> -p -o]\n");
        exit(-1);
    }
    readFile(filename, prnt, sourceCity, destCity);
    free(filename);
    free(prnt);
}

void readFile(char *filename, bool *prnt, int sourceCity, int destCity){
    FILE *fp;
    fp = fopen(filename,"r");
    if(fp == NULL){
        printf("Input file cannot be read!\n");
        exit(-1);
    }
    int N;
    int **G;
    fscanf(fp,"%d",&N);
    fscanf(fp,"\n");
    G = (int **)malloc(N*sizeof(int *));
    if (G == NULL) { // Check if memory allocation was successful
        perror( "Memory allocation failed.\n");
        exit(-1);
    }
    for(int i = 0; i<N; i++){
        *(G+i) = (int *)malloc(N*sizeof(int));
        if (*(G+i) == NULL) { // Check if memory allocation was successful
        perror( "Memory allocation failed.\n");
        free(G);
        exit(-1);
        }
        for(int j = 0; j<N; j++){
            int temp = fgetc(fp);
            if(temp == EOF){
                perror("File not complete");
                exit(-1);
            }
            if((char)temp == ' '){
                temp = fgetc(fp);
            }
            if (temp == '\n' || temp == '\r') {
                // Handle newline characters (UNIX or Windows style)
                // Skip them and continue to the next iteration
                j--; // Decrement j to read the same array index in the next iteration
                continue;
            }
            *(*(G + i) + j) = temp - '0';

        }
        int newline = fgetc(fp);
        if (newline != '\n' && newline != '\r' && newline != EOF) {
            perror("File not complete");
            exit(-1);
        }
        if (newline == EOF && i < N - 1) {
            perror("File not complete");
            exit(-1);
        }
    }    
    fclose(fp);  
    printNeighbor(G, N);
    printf("\n");
    if(prnt[0]){
        int **R ;
        int sizeR = createR(&R, G, &N);
        transitionalClosure(R, sizeR, prnt, sourceCity, destCity, filename);
        for(int i = 0; i<sizeR; i++)
            free(*(R+i));
        free(R);
    }
    for(int i = 0; i<N; i++)
        free(*(G+i));
    free(G);
}

int createR(int ***R, int **G, int *N){
   *R  = (int **)malloc(((*N)*(*N) -(*N))*sizeof(int *));
   if (*R == NULL) { // Check if memory allocation was successful
        perror( "Memory allocation failed.\n");
        exit(-1);
    }
   int sizeR = 0;
    for(int i=0; i<*N; i++){
        for(int j=0; j<*N; j++){
            if(G[i][j] == 1){
                *(*R+sizeR) = (int *)malloc(2*sizeof(int));
                if (*(*R+sizeR) == NULL) { // Check if memory allocation was successful
                    perror( "Memory allocation failed.\n");
                    free(R);
                    exit(-1);
                }
                (*R)[sizeR][0] = i;
                (*R)[sizeR][1] = j;
                sizeR++;
            }
        }
    }
    return sizeR;
}

void printR(int **R, int sizeR){
    printf("R* table\n");
    for(int i=0;i<sizeR;i++){
        printf("%d -> %d\n", R[i][0], R[i][1]);
    }
}

void printFile(char *filename, int **R, int size){    
    char *str1 = "out-";
    char *file = malloc(strlen(str1) + strlen(filename) + 1);
    if (file == NULL) { // Check if memory allocation was successful
        perror( "Memory allocation failed.\n");
        exit(-1);
    }
    strcpy(file, str1);
    strcat(file, filename);
    FILE *fp;
    fp = fopen(file,"w");
    if(fp == NULL){
        perror("Problem in file open");
        exit(-1);
    }
    fprintf(fp, "R* table\n");
    for(int i=0;i<size;i++){
        fprintf(fp, "%d -> %d\n", R[i][0], R[i][1]);
    }
    printf("Saving %s...\n", file);
    fclose(fp);
    free(file);
}

void printNeighbor(int **G, int N){
    printf("Neighbor table\n");
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++)
            printf("%d ", *(*(G + i)+j));
        printf("\n");
    }
}

int ** setR(int **R, int N, int flag, int u, int w){
    int **newR = (int **)malloc((N+flag)*sizeof(int *));
    if (newR == NULL) { // Check if memory allocation was successful
        perror( "Memory allocation failed.\n");
        exit(-1);
    }
    for(int i=0; i<N; i++){
        newR[i] = (int *)malloc(2*sizeof(int));
        if (newR[i] == NULL) { // Check if memory allocation was successful
            perror( "Memory allocation failed.\n");
            free(newR);
            exit(-1);
        }
        for(int j=0 ;j<2; j++){
            newR[i][j] = R[i][j];
        }
    }
    if(flag == 1){
        newR[N] = (int *)malloc(2*sizeof(int));
        if (newR[N] == NULL) { // Check if memory allocation was successful
            perror( "Memory allocation failed.\n");
            free(newR);
            exit(-1);
        }
        newR[N][0] = u;
        newR[N][1] = w;
    }
    return newR;
}

void transitionalClosure(int **R, int N, bool *prnt, int sourceCity, int destCity, char *filename){
    int **newR = setR(R, N, 0, 0, 0);
    bool flag = false;
    int sizeNewR = N;
    do{       
        for(int i=0; i<sizeNewR; i++){
            for(int j=0; j<N; j++){
                if(newR[i][1] == R[j][0]){
                    flag = true;
                    int u = newR[i][0] ;
                    int w = R[j][1];
                    if(u!=w && checkW(newR, u, w, sizeNewR)){
                        int **twinR = setR(newR, sizeNewR, 1, u, w);
                        newR = twinR;
                        sizeNewR++;                        
                    }
                    else{
                        flag = false;
                    }
                }
                else{
                    flag = false;
                }
            }
        }
    }while(flag);
    if(prnt[2])
        printR(newR, sizeNewR);
    if(prnt[1])
        if(checkLink(newR,sizeNewR, sourceCity, destCity) ){
            bool *visited = (bool *)malloc(N*sizeof(bool));
            findLink(R, N,sourceCity, destCity, visited, destCity);
            printf("\n");
            free(visited);
        }
    if (prnt[3])
        printFile(filename, newR, sizeNewR);   
    for(int i = 0; i<sizeNewR; i++)
        free(*(newR+i));
    free(newR);
}

bool checkW(int **newR, int u, int w, int N){
    for(int i=0; i<N; i++){
        if(newR[i][0] == u && newR[i][1] == w)
            return false;
    }
    return true;
}

bool checkLink(int **R, int size, int i, int j){
    for(int a=0; a<size; a++){
        if(R[a][0] == i && R[a][1] ==j){
            printf("Yes Path Exists!\n");
            return 1;
        }
    }
    printf("No Path Exists!\n");
    return 0;
}

bool findLink(int **R, int size, int i, int j, bool *visited, int dest) {
    if(i == j){
        printf("%d", j);
        return true;
    }
    visited[dest] = true;
    for (int a = size-1; a >=0; a--) {
        if (R[a][1] == j && !visited[R[a][0]]) {
            if(findLink(R, size, i, R[a][0], visited, j) ){
                printf(" => %d", j);
                return true;
            }
        }        
    }
    return false;
}