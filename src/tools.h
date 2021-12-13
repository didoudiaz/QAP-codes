/*
 *  Utilities for Local Search procedures
 *
 *  Copyright (C) 2002-2015 Daniel Diaz
 *
 *  tools.h: utilities - header file
 */

#ifndef _TOOLS_H
#define _TOOLS_H

long Real_Time(void);

long User_Time(void);


void Randomize_Seed(unsigned seed);

unsigned Randomize(void);

double Random_Double(void);

double Random_Double1(void);

unsigned Random(unsigned n);

int Random_Interval(int inf, int sup);

double Random_Interval_Double(double inf, double sup);


void Random_Array_Permut(int *vec, int size);


void Random_Permut(int *vec, int size, const int *actual_value, int base_value);

void Random_Permut_Repair(int *vec, int size, const int *actual_value, int base_value);

int Random_Permut_Check(int *vec, int size, const int *actual_value, int base_value);


#ifndef No_Gcc_Warn_Unused_Result
#define No_Gcc_Warn_Unused_Result(t) do { if(t) {} } while(0)
#endif

void *Malloc_Check(size_t size, char *src_file, int src_line);

void *Calloc_Check(size_t nb, size_t size, char *src_file,
                   int src_line);

void *Realloc_Check(void *ptr, size_t size, char *src_file, int src_line);

#define Malloc(size)       Malloc_Check(size, __FILE__, __LINE__)

#define Calloc(nb, size)   Calloc_Check(nb, size, __FILE__, __LINE__)

#define Realloc(ptr, size) Realloc_Check(ptr, size, __FILE__, __LINE__)

#define Free(ptr)          free(p@tr)


#define Alloc_Vector(vec, nb) vec = Malloc(nb * sizeof(*vec))

#define Alloc_Vector0(vec, nb) vec = Calloc(nb, sizeof(*vec))

#define Free_Vector(vec) Free(vec)

#define Copy_Vector(dst, src, nb)   memcpy((dst), (src), (nb) * sizeof(*dst))


void Fatal_Error(char *format, ...);

#endif /* _TOOLS_H */
