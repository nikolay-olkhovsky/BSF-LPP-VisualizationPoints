/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Visualization of Linear Programming Problem (ViLiPP)
Module: BSF-ProblemFunctions.h (Predefined Problem Function Forwards)
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "Problem-bsfTypes.h"	// Predefined Problem Types

void PC_bsf_CopyParameter(
	PT_bsf_parameter_T parameterIn,
	PT_bsf_parameter_T* parameterOutP
);
void PC_bsf_Start(
	bool* success
);
void PC_bsf_Init(
	bool* success,
	PT_bsf_parameter_T* parameter
);
void PC_bsf_IterOutput(
	PT_bsf_reduceElem_T* reduceResult, 
	int reduceCounter, 
	PT_bsf_parameter_T parameter,
	double elapsedTime, 
	int nextJob
);
void PC_bsf_IterOutput_1(
	PT_bsf_reduceElem_T_1* reduceResult, 
	int reduceCounter, 
	PT_bsf_parameter_T parameter,
	double elapsedTime, 
	int nextJob
);
void PC_bsf_IterOutput_2(
	PT_bsf_reduceElem_T_2* reduceResult, 
	int reduceCounter, 
	PT_bsf_parameter_T parameter,
	double elapsedTime, 
	int nextJob
);
void PC_bsf_IterOutput_3(
	PT_bsf_reduceElem_T_3* reduceResult, 
	int reduceCounter, 
	PT_bsf_parameter_T parameter,
	double elapsedTime, 
	int nextJob
);
void PC_bsf_JobDispatcher(
	PT_bsf_parameter_T* parameter,
	int* job,
	bool* exit,
	double t
);
//
void PC_bsf_MapF(
	PT_bsf_mapElem_T* mapElem, 
	PT_bsf_reduceElem_T* reduceElem, 
	int* success
);
void PC_bsf_MapF_1(
	PT_bsf_mapElem_T* mapElem, 
	PT_bsf_reduceElem_T_1* reduceElem, 
	int* success
);
void PC_bsf_MapF_2(
	PT_bsf_mapElem_T* mapElem, 
	PT_bsf_reduceElem_T_2* reduceElem, 
	int* success
);
void PC_bsf_MapF_3(
	PT_bsf_mapElem_T* mapElem, 
	PT_bsf_reduceElem_T_3* reduceElem, 
	int* success
);
void PC_bsf_ParametersOutput(
	PT_bsf_parameter_T parameter
);
void PC_bsf_ProblemOutput(
	PT_bsf_reduceElem_T* reduceResult,
	int reduceCounter,
	PT_bsf_parameter_T parameter,
	double t
);
void PC_bsf_ProblemOutput_1(
	PT_bsf_reduceElem_T_1* reduceResult, 
	int reduceCounter, 
	PT_bsf_parameter_T parameter,
	double t
);
void PC_bsf_ProblemOutput_2(
	PT_bsf_reduceElem_T_2* reduceResult, 
	int reduceCounter, 
	PT_bsf_parameter_T parameter,
	double t
);
void PC_bsf_ProblemOutput_3(
	PT_bsf_reduceElem_T_3* reduceResult, 
	int reduceCounter, 
	PT_bsf_parameter_T parameter,
	double t
);
void PC_bsf_ProcessResults(
	PT_bsf_reduceElem_T* reduceResult, 
	int reduceCounter, 
	PT_bsf_parameter_T* parameter,
	int* nextJob,
	bool* exit 
);
void PC_bsf_ProcessResults_1(
	PT_bsf_reduceElem_T_1* reduceResult, 
	int reduceCounter, 
	PT_bsf_parameter_T* parameter,
	int* nextJob,
	bool* exit 
);
void PC_bsf_ProcessResults_2(
	PT_bsf_reduceElem_T_2* reduceResult, 
	int reduceCounter, 
	PT_bsf_parameter_T* parameter,
	int* nextJob,
	bool* exit 
);
void PC_bsf_ProcessResults_3(
	PT_bsf_reduceElem_T_3* reduceResult, 
	int reduceCounter, 
	PT_bsf_parameter_T* parameter,
	int* nextJob,
	bool* exit 
);
void PC_bsf_ReduceF(
	PT_bsf_reduceElem_T* x, 
	PT_bsf_reduceElem_T* y, 
	PT_bsf_reduceElem_T* z
);
void PC_bsf_ReduceF_1(
	PT_bsf_reduceElem_T_1* x, 
	PT_bsf_reduceElem_T_1* y, 
	PT_bsf_reduceElem_T_1* z
);
void PC_bsf_ReduceF_2(
	PT_bsf_reduceElem_T_2* x, 
	PT_bsf_reduceElem_T_2* y, 
	PT_bsf_reduceElem_T_2* z
);
void PC_bsf_ReduceF_3(
	PT_bsf_reduceElem_T_3* x, 
	PT_bsf_reduceElem_T_3* y, 
	PT_bsf_reduceElem_T_3* z
);
void PC_bsf_SetListSize(
	int* listSize
);
void PC_bsf_SetInitParameter(
	PT_bsf_parameter_T* parameter
);
void PC_bsf_SetMapListElem(
	PT_bsf_mapElem_T* elem, 
	int i
);