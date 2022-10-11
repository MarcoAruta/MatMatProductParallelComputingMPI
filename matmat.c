#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#define DIM 2

void BMR(int, int, int, double**, double**, double**, int*, MPI_Comm *, MPI_Comm *, MPI_Comm *);
void localProduct(double**, double**, double**, int, int);
void createResult(double**, double**, int, int, int, int);
void createMat(double***, int, int, bool);
void createGrid(MPI_Comm*, MPI_Comm*, MPI_Comm*, int, int, int, int, int*);
void splitMatrix(int, int, double**, double**, int, int);
void copyMatrix(double**, double**, int, int);


int main(int argc, char* argv[]){
    srand(time(NULL));
    int menum, nproc, check, dimMat, dimGrid, dimSubatrix, i, j, k, *coordinate;
	double **matrixA, **MatrixB, **result, **partialResult, **submatrixA, **submatrixB, startTime, stopTime;
    int submatrixParameters[4];
	MPI_Comm grid, gridr, gridc;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &menum);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    /*
        Check dei parametri in input da parte del processore
    */
    if(menum == 0){
		if(argc != 2){
			printf("Errore: numeri di parametri non valido.\n");
			check = -1;
		}
		else{
			dimMat = atoi(argv[1]);
            dimGrid = sqrt(nproc);

			//check integrità parametri e validità matrice
			if(nproc != dimGrid*dimGrid || dimMat % dimGrid != 0){ 
				printf("Errore: impossibile accettare parametri in input.\n");
				check = -1;
			}
			else{
				check = 1;
                
				createMat(&matrixA, dimMat, dimMat, true);	
                createMat(&MatrixB, dimMat, dimMat, true);
                createMat(&result, dimMat, dimMat, false);
                
                /*printf("matrixA:\n");
                for(i = 0; i < dimMat; i++){
                    for(j = 0; j < dimMat; j++)
                        printf("%.2lf ", matrixA[i][j]);
                    printf("\n");
                }
                printf("\n");printf("\n");

                printf("MatrixB:\n");
                for(i = 0; i < dimMat; i++){
                    for(j = 0; j < dimMat; j++)
                        printf("%.2lf ", MatrixB[i][j]);
                    printf("\n");
                }
                printf("\n");printf("\n");  */			}
        }
	}
    MPI_Bcast(&check, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(check == 1){
        MPI_Bcast(&dimMat, 1, MPI_INT, 0, MPI_COMM_WORLD);        
        MPI_Bcast(&dimGrid, 1, MPI_INT, 0, MPI_COMM_WORLD);

        dimSubatrix = dimMat/dimGrid;

        createMat(&partialResult, dimSubatrix, dimSubatrix, false);        
        submatrixA = (double**) calloc(dimSubatrix, sizeof(double*));
        for(i = 0; i < dimSubatrix; i++) 
			submatrixA[i] = (double*) calloc(dimSubatrix, sizeof(double));
        submatrixB = (double**) calloc(dimSubatrix, sizeof(double*));
        for(i = 0; i < dimSubatrix; i++)	
			submatrixB[i] = (double*) calloc(dimSubatrix, sizeof(double));

        splitMatrix(menum, nproc, matrixA, submatrixA, dimSubatrix, dimGrid);
        splitMatrix(menum, nproc, MatrixB, submatrixB, dimSubatrix, dimGrid);

        coordinate = (int*) calloc(2, sizeof(int));
        createGrid(&grid, &gridr, &gridc, menum, nproc, dimGrid, dimGrid, coordinate);

        startTime = MPI_Wtime();
        BMR(menum, dimSubatrix, dimGrid, partialResult, submatrixA, submatrixB, coordinate, &grid, &gridr, &gridc);
        createResult(partialResult, result, menum, nproc, dimMat, dimSubatrix);
        stopTime = MPI_Wtime();

        if(menum == 0){
             /*printf("result:\n");
            for(i = 0; i < dimMat; i++){
                for(j = 0; j < dimMat; j++)
                    printf("%.2lf ", result[i][j]);
            printf("\n");
            } */
            printf("\nTime: %lf seconds\n", stopTime-startTime);
        }
    }

    MPI_Finalize();
    return 0;
}

//Funzione creazione matrice rows*columns
void createMat(double*** matrix, int rows, int columns, bool fill){
    int i, j;

    //allocazione matrice
    *matrix = (double**) calloc(rows, sizeof(double*));
    for(i = 0; i < rows; i++)
        (*matrix)[i] = (double*) calloc(columns, sizeof(double));

    //popolamento matrice
    if(fill){
        for(i = 0; i < rows; i++){
            for(j = 0; j < columns; j++)
                (*matrix)[i][j] = (double)rand() / RAND_MAX * 150;
        }
    }

    return;
}

//Funzione creazione griglia
void createGrid(MPI_Comm *grid , MPI_Comm *gridr, MPI_Comm *gridc , int menum , int nproc , int row, int colGrid, int *coordinate){
	int *ndim, reorder, *period, vc[2];
	
	ndim = (int*) calloc(DIM, sizeof(int));
	ndim[0] = row;
	ndim[1] = colGrid;
	period = (int*) calloc(DIM, sizeof(int));
	period[0] = period[1] = 1;
	reorder = 0;
	
	MPI_Cart_create(MPI_COMM_WORLD, DIM, ndim, period, reorder, grid);
	MPI_Cart_coords(*grid, menum, DIM, coordinate);
	
	// Creazione sottogrid
	vc[0] = 0;
	vc[1] = 1;
	MPI_Cart_sub(*grid, vc, gridr);
	
	vc[0] = 1;
	vc[1] = 0;
	MPI_Cart_sub(*grid, vc, gridc); 
	
	return;
}

//Funzione di inoltro ai processori di sottomatrici da parte del processore 0 
void splitMatrix(int menum, int nproc, double** matrixA, double** submatrixA, int dimSubatrix, int dimGrid){
    int i, j, k, l, m, count, startingLine, startingColumn, tag;
    double submatrix[dimSubatrix][dimSubatrix];
    MPI_Status status;
    
    if(menum == 0){
        for(i = 0; i < dimSubatrix; i++){
            for(j = 0; j < dimSubatrix; j++)
                submatrixA[i][j] = matrixA[i][j];
        }

        for(i = 1; i < nproc; i++){
            //popolamento sottomatrici
            l = 0, m = 0, count = 0;
            startingLine = dimSubatrix * (i / dimGrid);
            startingColumn = dimSubatrix * (i % dimGrid);
            k = startingLine; j = startingColumn;

            while(count < dimSubatrix*dimSubatrix){
                submatrix[l][m] = matrixA[k][j];
                count++;
                j++; m++;
                if(m == dimSubatrix){ 
					j = startingColumn; 
					k++; 
					l++; 
					m = 0;
				}
            }
            //Inoltro sottomatrice
            for(j = 0; j < dimSubatrix; j++){
                tag = 25 + i;
                MPI_Send(submatrix[j], dimSubatrix, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
            }
        }
    }
    else{
        //Ricezione sottomatrice da parte dei processori
        for(j = 0; j < dimSubatrix; j++){
            tag = 25 + menum;
            MPI_Recv(submatrixA[j], dimSubatrix, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        }
    } 

    return;
}

void BMR(int menum, int dimSubatrix, int dimGrid, double** partialResult, double** submatrixA, double** submatrixB, int* coordinate, MPI_Comm *grid, MPI_Comm *gridr, MPI_Comm *gridc){
    int step, tag, j, i, sender, menumRow, menumCol, senderCol, receiverCol;
    double **bufferA;
    MPI_Status status;
    createMat(&bufferA, dimSubatrix, dimSubatrix, false);

    MPI_Comm_rank(*gridc, &menumCol);
    MPI_Comm_rank(*gridr, &menumRow);

    /*
        Ad ogni step:
         1. il processore sender invia la sottomatrice A agli altri processori sulla riga;
         2. ogni processore moltiplica la sottomatrice A che ha ricevuto/copiato nel caso del step 0 (bufferA) per la sottomatrice B che possiede;
         3. ogni processore riceve una sottomatrice B dal processore nella stessa colonna, nella riga successiva ed aggiorna la propria. 
    */
    for(step = 0; step < dimGrid; step++){
        if(coordinate[1] == (coordinate[0] + step) % dimGrid){ 
            sender = menumRow;
            copyMatrix(bufferA, submatrixA, dimSubatrix, dimSubatrix);
        }
        else
            sender = (coordinate[0] + step) % dimGrid;
        
        for(j = 0; j < dimSubatrix; j++)
            MPI_Bcast(bufferA[j], dimSubatrix, MPI_DOUBLE, sender, *gridr);
      
        localProduct(bufferA, submatrixB, partialResult, dimSubatrix, dimSubatrix);

		if (menumCol-1 < 0)
			receiverCol = (menumCol-1) + dimGrid;
		else
			receiverCol = (menumCol-1) % dimGrid;
			
        if (menumCol+1 < 0)
			senderCol = (menumCol+1) + dimGrid;
		else
			senderCol = (menumCol+1) % dimGrid;
        for(j = 0; j < dimSubatrix; j++){
            MPI_Send(submatrixB[j], dimSubatrix, MPI_DOUBLE, receiverCol, 20 + receiverCol, *gridc);
            MPI_Recv(submatrixB[j], dimSubatrix, MPI_DOUBLE, senderCol, 20 + menumCol, *gridc, &status);
        }
    }
}

void localProduct(double** m1, double** m2, double** res, int colsM1, int rowsM2){
    int i, j, k;

    for(i = 0; i < rowsM2; i++){
        for(j = 0; j < colsM1; j++){
            for(k = 0; k < colsM1; k++)
                res[i][j] += m1[i][k]*m2[k][j];
        }
    }

    return;
}


void copyMatrix(double** m1, double** m2, int rowsM2, int colsM2){
    int i, j;
    for(i = 0; i < rowsM2; i++){
        for(j = 0; j < colsM2; j++)
            m1[i][j] = m2[i][j]; //m1: risultato
    }

    return;
}

void createResult(double** partial, double** final, int menum, int nproc, int dimMat, int dimSubatrix){
    int i, j, l,p,k;
    MPI_Status status;
    double **A;

    createMat(&A, dimSubatrix,dimSubatrix, false);

    if(menum == 0){
        for(i = 0; i < dimSubatrix; i++){
            for(j = 0; j < dimSubatrix; j++)
                final[i][j] = partial[i][j];
        }

        i = 1;p=0;
        while(i < nproc){
            k = (i * dimSubatrix) % dimMat;
            for(j = 0; j < dimSubatrix; j++){
                MPI_Recv(A[j], dimSubatrix, MPI_DOUBLE, i, 20 + i, MPI_COMM_WORLD, &status);
            }
            for(l = 0; l < dimSubatrix; l++){
                for(j = 0; j < dimSubatrix; j++)
                    final[j+p][k] = A[j][l];
                
            k++;
            if(k==dimMat) p=p+dimSubatrix;
            }
            ++i;
        }
    }
    else{
        for(i = 0; i < dimSubatrix; i++)
            MPI_Send(partial[i], dimSubatrix, MPI_DOUBLE, 0, 20 + menum, MPI_COMM_WORLD);
    }

    return;
}
