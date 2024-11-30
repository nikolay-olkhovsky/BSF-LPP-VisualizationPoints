/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Visualization of Linear Programming Problem (ViLiPP)
Module: Problem-bsfCode.cpp (Implementation of the Problem)
Prefix: PC
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
// PP_STATE_NEW_PROBLEM
// reads new problem from LPP file AND first trace point
// PP_STATE_NEW_POINT
// reads only new trace point

#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	parameter->k = 0;
}

void PC_bsf_Start(bool* success) {
	ini::IniFile config;
	int packageSize = 0;

	config.load(PP_FILE_INI);
	PP_PATH = config["general"]["PP_PATH"].as<string>();
	PP_PROBLEM_NAME = config["general"]["PP_PROBLEM_NAME"].as<string>();
	PP_MTX_PREFIX = config["general"]["PP_MTX_PREFIX"].as<string>();
	PP_MTX_POSTFIX_TR = config["general"]["PP_MTX_POSTFIX_TR"].as<string>();

	PP_N = config["general"]["PP_N"].as<int>();
	PP_LPP_FILE = config["general"]["PP_LPP_FILE"].as<string>();
	PP_ETA = config["visualization"]["PP_ETA"].as<int>();
	PP_DELTA = config["visualization"]["PP_DELTA"].as<double>();
	PP_MTX_POSTFIX_RET = config["visualization"]["PP_MTX_POSTFIX_RET"].as<string>();
	PP_MTX_POSTFIX_IM = config["visualization"]["PP_MTX_POSTFIX_IM"].as<string>();
	PP_IMAGE_OUT = config["visualization"]["PP_IMAGE_OUT"].as<bool>();
	PP_RECEPTIVE_FIELD_OUT = config["visualization"]["PP_RECEPTIVE_FIELD_OUT"].as<bool>();
	PP_CROSSFILED = config["visualization"]["PP_CROSSFILED"].as<bool>();

	PD_problemName = PP_PROBLEM_NAME;
	PD_id = 1;
	PD_currentProblem = 1;
	PD_initState = true;
	PD_time = clock();
	PD_images.clear();

#ifdef PP_DATABASE_OUTPUT
	//storage.sync_schema(true);
	try {
		if(BSF_sv_mpiRank == BSF_sv_mpiMaster)
			storage.sync_schema(true);
		PD_problemIds = storage.select(&Problem::id);
		PD_problemsNumber = PD_problemIds.size();
		PD_problems = storage.get_all<Problem>();
		PD_inequalities = storage.get_all<Inequality>();
		PD_points = storage.get_all<Precedent>();
		cout << "[" << BSF_sv_mpiRank << "] : " << "DB read success." << endl;
	}
	catch (...) {
		cout << "[" << BSF_sv_mpiRank << "] : " << "DB read failure." << endl; 
		try {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				storage.sync_schema(true);
			PD_problemIds = storage.select(&Problem::id);
			PD_problemsNumber = PD_problemIds.size();
			PD_problems = storage.get_all<Problem>();
			PD_inequalities = storage.get_all<Inequality>();
			PD_points = storage.get_all<Precedent>();
			cout << "[" << BSF_sv_mpiRank << "] : " << "DB second read success." << endl;
		}
		catch (...) {
			cout << "[" << BSF_sv_mpiRank << "] : " << "DB second read failure." << endl;
		}
	}
#else
	PD_lppFilename = PP_PATH;
	PD_lppFilename += PP_LPP_FILE;
	PD_stream_lppFile = fopen(PD_lppFilename.c_str(), "r");
	if (PD_stream_lppFile == NULL) {
		cout << '[' << BSF_sv_mpiRank << "]: Failure of opening file '" << PD_lppFilename << "'.\n";
		*success = false;
		return;
	}

	if (fscanf(PD_stream_lppFile, "%d", &PD_problemsNumber) < 1) {
		cout << '[' << BSF_sv_mpiRank << "]: Can't read package size in " << PD_lppFilename << endl;
		*success = false;
		return;
	}

	PD_traceFilename = PP_PATH;
	PD_traceFilename += PP_MTX_PREFIX;
	PD_traceFilename += PD_problemName;
	PD_traceFilename += PP_MTX_POSTFIX_TR;
	PD_stream_traceFile = fopen(PD_traceFilename.c_str(), "r");
	if (PD_stream_traceFile == NULL) {
		cout << '[' << BSF_sv_mpiRank << "]: Failure of opening file '" << PD_traceFilename << "'.\n";
		*success = false;
		return;
	}
	if (fscanf(PD_stream_traceFile, "%d", &packageSize) < 1) {
		cout << '[' << BSF_sv_mpiRank << "]: Can't read package size in " << PD_traceFilename << endl;
		*success = false;
		return;
	}
	if (packageSize != PD_problemsNumber) {
		cout << '[' << BSF_sv_mpiRank << "]: Wrong package size in " << PD_traceFilename << endl;
		*success = false;
		return;
	}

	PD_outFilename = PP_PATH;
	PD_outFilename += PP_MTX_PREFIX;
	PD_outFilename += PD_problemName;
	PD_outFilename += PP_MTX_POSTFIX_IM;
	PD_stream_outFile = fopen(PD_outFilename.c_str(), "w");
	if (PD_stream_outFile == NULL) {
		cout << '[' << BSF_sv_mpiRank << "]: Failure of opening file '" << PD_outFilename << "'.\n";
		*success = false;
		return;
	}

	//--------------- Opening retina file ------------------
	PD_retFilename = PP_PATH;
	PD_retFilename += PP_MTX_PREFIX;
	PD_retFilename += PD_problemName;
	PD_retFilename += PP_MTX_POSTFIX_RET;

	PD_stream_retFile = fopen(PD_retFilename.c_str(), "w");
	if (PD_stream_retFile == NULL) {
		cout << '[' << BSF_sv_mpiRank << "]: Failure of opening file '" << PD_retFilename << "'.\n";
		*success = false;
		return;
	}
	/*
	if (fprintf(PD_stream_retFile, "%d\n", PD_problemsNumber) < 1) {
		//
		cout
			<< "Can't write package size to " << PD_retFilename << endl;
		*success = false;
		return;
	}
	*/
#endif // PP_DATABASE_OUTPUT
}

void PC_bsf_Init(bool* success, PT_bsf_parameter_T* parameter) {
	float buf;
	int readNumber, pid;

	if (PD_initState) {
		PD_read_state = PP_STATE_NEW_PROBLEM;
		PD_initState = false;
		parameter->k = 0;
	}
	else
		PD_read_state = parameter->state;

	switch (PD_read_state) {
	case PP_STATE_NEW_PROBLEM:
		do {
#ifdef PP_DATABASE_OUTPUT
			vector<double> _buff;
			vector<Inequality> inequalities;
			//auto problem = storage.get<Problem>(PD_problemIds[PD_currentProblem]);
			auto problem = PD_problems[PD_currentProblem];
			//auto inequalities = storage.get_all<Inequality>(sqlite_orm::where(sqlite_orm::c(&Inequality::problem_id) == problem.id));
			
			copy_if(PD_inequalities.begin(), PD_inequalities.end(), back_inserter(inequalities), [&problem](const Inequality& item)->bool {return item.problem_id == problem.id; });
			//PD_traceIds = storage.select(&SurfacePoint::id, sqlite_orm::where(sqlite_orm::c(&SurfacePoint::problem_id) == problem.id));
			PD_traces.clear();
			copy_if(PD_points.begin(), PD_points.end(), back_inserter(PD_traces), [&problem](const Precedent& item)->bool {return item.problem_id == problem.id; });

			//auto currentPoint = storage.get<SurfacePoint>(PD_traceIds[0]);
			auto currentPoint = PD_traces[0];
			PD_tracesNumber = PD_traces.size();
			if (PD_tracesNumber < 1) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "[" << BSF_sv_mpiRank << "] :"
					<< "Wrong trace for problem ID: " << pid << endl
					<< "Traces number: " << PD_tracesNumber << " (must be > 2)." << endl;
				PD_currentProblem++;
				PD_currentTrace = 0;
				continue;
			}

			PD_n = problem.N;
			PD_m = inequalities.size();
			if (PD_m > PP_MAX_M || PD_n > PP_MAX_N) {
				cout << "[" << BSF_sv_mpiRank << "] :" << endl
					<< "Wrong problem parameters!" << endl
					<< "Number of inequalities: " << PD_m << " (max: " << PP_MAX_M << ")." << endl
					<< "Number of dimensions: " << PD_n << " (max: " << PP_MAX_N << ")." << endl;
				*success = false;
				return;
			}
			PD_m = 0;

			for (unsigned i = 0; i < PD_n; i++) {
				for (unsigned j = 0; j < PD_n; j++)
					if (i == j) PD_A[PD_m][j] = 1.;
					else		PD_A[PD_m][j] = 0.;
				PD_b[PD_m] = problem.high;
				PD_m++;
			}
			for (auto& inequality : inequalities) {
				_buff = charToDouble(inequality.coefficients);
				for (unsigned j = 0; j < PD_n; j++)
					PD_A[PD_m][j] = _buff[j];
				PD_b[PD_m] = inequality.b;
				PD_m++;
			}
			for (unsigned i = 0; i < PD_n; i++) {
				for (unsigned j = 0; j < PD_n; j++)
					if (i == j) PD_A[PD_m][j] = -1.;
					else		PD_A[PD_m][j] = 0.;
				PD_b[PD_m] = problem.low;
				PD_m++;
			}

			_buff = charToDouble(problem.c);
			for (unsigned i = 0; i < PD_n; i++) {
				PD_c[i] = _buff[i];
			}

			//PD_currentPointId = currentPoint.id;
			_buff = charToDouble(currentPoint.coefficients);
			for (unsigned i = 0; i < PD_n; i++) {
				PD_nextPoint[i] = _buff[i];
			}
#else
			if (fscanf(PD_stream_lppFile, "%d%d%d", &pid, &PD_m, &PD_n) == 0) {
				cout << '[' << BSF_sv_mpiRank << "]: Unexpected end of LPP file header" << endl;
				*success = false;
				return;
			}

			if (pid != PD_currentProblem + 1) {
				cout << "[" << BSF_sv_mpiRank << "] :"
					<< "Wrong problem ID in LPP file: " << pid
					<< " (expected: " << PD_currentProblem + 1 << ")" << endl;
				*success = false;
				return;
			}

			if (PD_m > PP_MAX_M || PD_n > PP_MAX_N) {
				cout << "[" << BSF_sv_mpiRank << "] :" << endl
					<< "Wrong problem parameters!" << endl
					<< "Number of inequalities: " << PD_m << " (max: " << PP_MAX_M << ")." << endl
					<< "Number of dimensions: " << PD_n << " (max: " << PP_MAX_N << ")." << endl;
				*success = false;
				return;
			}

			for (int i = 0; i < PD_m; i++) {
				for (int j = 0; j < PD_n; j++) {
					if (fscanf(PD_stream_lppFile, "%f", &buf) == 0) {
						cout << '[' << BSF_sv_mpiRank << "]: Unexpected end of LPP file (A row)" << endl;
						*success = false;
						//				system("pause");
						return;
					}
					PD_A[i][j] = buf;
				}
				if (fscanf(PD_stream_lppFile, "%f", &buf) == 0) {
					cout << '[' << BSF_sv_mpiRank << "]: Unexpected end of LPP file (b coefficient)" << endl;
					*success = false;
					//			system("pause");
					return;
				}
				PD_b[i] = buf;
			}

			for (int j = 0; j < PD_n; j++) {
				if (fscanf(PD_stream_lppFile, "%f", &buf) == 0) {
					cout << '[' << BSF_sv_mpiRank << "]: Unexpected end of LPP file (c vector)" << endl;
					*success = false;
					//			system("pause");
					return;
				}
				PD_c[j] = buf;
				//PD_z[j] = buf + 1.;
			}

			if (fscanf(PD_stream_traceFile, "%d%d%d", &pid, &PD_tracesNumber, &readNumber) == 0) {
				cout << '[' << BSF_sv_mpiRank << "]: Unexpected end of TRACE file header" << endl;
				*success = false;
				return;
			}

			if (pid != PD_currentProblem + 1) {
				cout << "[" << BSF_sv_mpiRank << "] :"
					<< "Wrong problem ID in TRACES file: " << pid
					<< " (expected: " << PD_currentProblem + 1 << ")" << endl;
				*success = false;
				return;
			}

			if (readNumber != PD_n) {
				cout << "[" << BSF_sv_mpiRank << "] :" << endl
					<< "Wrong traces file format!" << endl
					<< "Dimension of traces: " << readNumber << " (PD_n: " << PD_n << ")." << endl
					<< "Problem ID: " << PD_currentProblem << endl;
				*success = false;
				return;
			}

			if (PD_tracesNumber < 1) {
				if(BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "[" << BSF_sv_mpiRank << "] :"
					<< "Wrong trace for problem ID: " << pid << endl
					<< "Traces number: " << PD_tracesNumber << " (must be > 2)." << endl;
				PD_currentProblem++;
				PD_currentTrace = 0;
				continue;
			}

			for (int i = 0; i < PD_n; i++) {
				if (fscanf(PD_stream_traceFile, "%f", &buf) == 0) {
					cout << "[" << BSF_sv_mpiRank << "]: Unexpected end of TRACE file (first point)" << endl;
					*success = false;
					return;
				}
				PD_nextPoint[i] = buf; // Read first point in the trace
			}

			if(PD_tracesNumber < 2 && BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "[" << BSF_sv_mpiRank << "] :"
				<< "Wrong trace for problem ID: " << pid << endl
				<< "Traces number: " << PD_tracesNumber << " (must be > 2)." << endl;
#endif //PP_DATABASE_OUTPUT
			PD_currentProblem++;
			PD_currentTrace = 1;
		} while (PD_tracesNumber < 2);

	case PP_STATE_NEW_POINT:
#ifdef PP_DATABASE_OUTPUT
		if (PD_currentTrace < PD_tracesNumber) {
			vector<double> _buff;
			vector<double> _indices;
			//auto currentPoint = storage.get<SurfacePoint>(PD_traceIds[PD_currentTrace + 1]);
			auto currentPoint = PD_traces[PD_currentTrace];
			auto previousPoint = PD_traces[PD_currentTrace - 1];

			_buff = charToDouble(currentPoint.coefficients);
			for (unsigned i = 0; i < PD_n; i++) {
				PD_z[i] = PD_nextPoint[i]; // Move Z to next point in the trace
				PD_nextPoint[i] = _buff[i]; // And read forward next point in the trace
				PD_answerVector[i] = PD_nextPoint[i] - PD_z[i];
			}
			_indices = charToDouble(previousPoint.face_numbers);
			PD_hyperplaneIndices.clear();
			for (double index : _indices)
				PD_hyperplaneIndices.push_back(int(index));
		}
		else {
			for (unsigned i = 0; i < PD_n; i++) {
				PD_z[i] = PD_nextPoint[i]; // Move Z to next point in the trace
				PD_nextPoint[i] = PD_z[i]; // And read forward next point in the trace
				PD_answerVector[i] = 0.;
			}
		}
#else
		for (int i = 0; i < PD_n; i++) {
			PD_z[i] = PD_nextPoint[i]; // Move Z to next point in the trace
			if (fscanf(PD_stream_traceFile, "%f", &buf) == 0) {
				cout << '[' << BSF_sv_mpiRank << "]: Unexpected end of TRACE file (next point)" << endl;
				*success = false;
				return;
			}
			PD_nextPoint[i] = buf; // And read forward next point in the trace
			PD_answerVector[i] = PD_nextPoint[i] - PD_z[i];
		}
#endif //PP_DATABASE_OUTPUT
		basis_Init();
		if (PD_currentTrace < PD_tracesNumber) {
			fieldProjection(PD_nextPoint, PD_fieldVector);
			coordinateAngles(PD_fieldVector, PD_cosVector);
			norm_Vector(PD_answerVector);
		}
		else
			for (unsigned i = 0; i < PD_n; i++) {
				PD_fieldVector[i] = 0.;
				PD_cosVector[i] = 0.;
				PD_answerVector[i] = 0.;
			}
		PD_currentTrace++;
	}
//	basis_Init();
	PD_K = 1;
	if (!PP_CROSSFILED)
		for (int i = 0; i < PD_n - 1; i++)
			PD_K *= (2 * PP_ETA + 1);
	else
		PD_K = 2 * PP_ETA * (PD_n - 1) + 1;

	if(PP_IMAGE_OUT && BSF_sv_mpiRank == BSF_sv_mpiMaster)
		PD_I = new PT_float_T[PD_K];

	if(PP_RECEPTIVE_FIELD_OUT && BSF_sv_mpiRank == BSF_sv_mpiMaster)
		PD_field = new PT_float_T* [PD_K];
}

void PC_bsf_SetListSize(int* listSize) {
	//*listSize = (int)PD_m;
	*listSize = (int)PD_K;
}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	elem->k = i;
}

// 0. Pseudo-pojection
void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success // 1 - reduceElem was produced successfully; 0 - otherwise
) {
	PT_vector_T target;
	PT_float_T distance = -FLT_MAX;
	PT_float_T newDistance = -FLT_MAX;
	int k = mapElem->k;
	G(mapElem->k, PD_g, PP_CROSSFILED);
	//targetProjection(i, PD_g, target);
	//for (int i = 0; i < PD_m; i++) {
	for (int i : PD_hyperplaneIndices) {
		if (dotproduct_Vector(PD_A[i], PD_c) > 0/* && isInnerPoint(target)*/) {
			//if (i == 2) {
				//reduceElem->objectiveDistance = objectiveDistance(target);
			newDistance = bias(i);
		}
		if (newDistance > distance)
			distance = newDistance;
	}
	try {
		reduceElem->distances = new list<PT_float_T>();
		reduceElem->distances->push_back(distance);
	}
	catch (...)	{
		cout << "[" << BSF_sv_mpiRank << "] : " << "failure on distance pushing." << endl;
	}
}

// 1. CheckIn
void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem, int* success) {

}

void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem, int* success) {
	// not used
}

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem, int* success) {
	// not used
}

// 0. Pseudo-pojection
void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	try {
		x->distances->splice(x->distances->end(), *(y->distances));
		delete y->distances;
	}
	catch (const std::exception& e) {
		cout << "[" << BSF_sv_mpiRank << "] : " << e.what();
	}
}

// 1. CheckIn
void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {

}

void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {
	// not used
}

void PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {
	// not used
}

//0. Start
void PC_bsf_ProcessResults(
	PT_bsf_reduceElem_T* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	if (PP_IMAGE_OUT) {
		std::list<PT_float_T>::iterator it = reduceResult->distances->begin();
		if (reduceResult->distances->size() != PD_K)
			cout << "Error: reduceResult->distances.size() != PD_K" << endl;
		for (PT_integer_T k = 0; k < PD_K; k++) {
			PD_I[k] = *it++;
		}
	}

	if (PP_RECEPTIVE_FIELD_OUT) {
		for (PT_integer_T k = 0; k < PD_K; k++) {
			PD_field[k] = new PT_float_T[PD_n];
			G(k, PD_field[k], PP_CROSSFILED);
		}
	}

	delete reduceResult->distances;
//	parameter->k += 1;
//	if (parameter->k >= PD_K) {
		if (PD_currentTrace - 1 < PD_tracesNumber) {
			parameter->state = PP_STATE_NEW_POINT;
			*nextJob = BD_JOB_RESET;
		}
		else if (PD_currentProblem < PD_problemsNumber) {
		//else if (PD_currentProblem < 3) {
			parameter->state = PP_STATE_NEW_PROBLEM;
			*nextJob = BD_JOB_RESET;
		}
		else {
			storage.insert_range(PD_images.begin(), PD_images.end());
			PD_images.clear();
			*exit = true;
		}
//	}
}

// 1. Movement on Polytope  ========================================================
void PC_bsf_ProcessResults_1(
	PT_bsf_reduceElem_T_1* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {

}

void PC_bsf_ProcessResults_2(
	PT_bsf_reduceElem_T_2* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
}

void PC_bsf_ProcessResults_3(
	PT_bsf_reduceElem_T_3* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
}

void PC_bsf_JobDispatcher(
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* job,
	bool* exit,
	double t
) {
	*job = 0;
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
//	cout << "Problem " << PD_currentProblem << " of " << PD_problemsNumber << ", trace " << PD_currentTrace << " of " << PD_tracesNumber << endl;
#ifdef PP_BSF_ITER_OUTPUT
	cout << "Number of receptive points: " << PD_K << endl;
	cout << "Trace point: ";
	for (int i = 0; i < PD_n; i++) {
		cout << PD_z[i] << " ";
	}
	cout << endl;
	cout << "Objective function coordinates: ";
	for (int i = 0; i < PD_n; i++) {
		cout << PD_c[i] << " ";
	}
	cout << endl;
	basis_Print();
	cout << endl;
#endif // PP_BSF_ITER_OUTPUT
	/*	cout << "============================================== Problem parameters ===============================================" << endl;
		cout << "Problem ID: " << PD_currentProblem << endl;
		cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;
	#ifdef PP_BSF_OMP
	#ifdef PP_BSF_NUM_THREADS
		cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
	#else
		cout << "Number of Threads: " << omp_get_num_procs() << endl;
	#endif // PP_BSF_NUM_THREADS
	#else
		cout << "OpenMP is turned off!" << endl;
	#endif // PP_BSF_OMP
		cout << "Number of problems: " << PD_problemsNumber << endl;
		cout << "Dimensions: " << PD_n << ", max = " << log(PP_MAX_K) / log(2 * PP_ETA + 1) << endl;
		cout << "Number of inequalities: " << PD_m << endl;
		cout << "Receptive field rank: " << PP_ETA << endl;
		cout << "Receptive field density: " << PP_DELTA << endl;
		cout << "Maximum number of points: " << PD_K << endl;
		cout << "Receptive field coordinates: " << endl;
		for (int i = 0; i < PD_n; i++) {
			cout << PD_z[i] << " ";
		}
		cout << endl;
		cout << "Objective function coordinates: " << endl;
		for (int i = 0; i < PD_n; i++) {
			cout << PD_c[i] << " ";
		}
		cout << endl;
		basis_Print();
		cout << "Matrix of inequalities: " << endl;
		for (int i = 0; i < PD_m; i++) {
			print_Vector(PD_A[i]); cout << endl;
		}
		/**/
}

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	parameterOutP->k = parameterIn.k;
}

// 0. Start
void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob) {

}

// 1. Movement on Polytope
void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{

}

// 2.
void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// not used
}


void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// not used
}

// 0. Start
void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// Output precedents
	if (PD_currentProblem % 100 == 0) {
		storage.insert_range(PD_images.begin(), PD_images.end());
		PD_images.clear();
		cout << PD_currentProblem << " problems processed. Time: " << (clock() - PD_time) / CLOCKS_PER_SEC << endl;
	}
	if (PP_IMAGE_OUT) {
#ifdef PP_DATABASE_OUTPUT
		Image newImage;
		vector<double> _buff(PD_K);

		newImage.id = 0;
		newImage.precedent_id = PD_traces[PD_currentTrace - 2].id;
		newImage.density = PP_DELTA;
		newImage.rank = PP_ETA;

		_buff.resize(PD_n);
		for (int i = 0; i < PD_n; i++)
			_buff[i] = PD_answerVector[i];
		newImage.answer_vector = doubleToChar(_buff);

		_buff.resize(PD_n);
		//cout << "problem = " << PD_currentProblem << " trace = " << PD_currentTrace << endl;
		//cout << "cos vector: ";
		for (int i = 0; i < PD_n; i++) {
			_buff[i] = PD_cosVector[i];
			//cout << PD_cosVector[i] << ' ';
		}
		//cout << endl;
		newImage.cosine_vector = doubleToChar(_buff);

		newImage.num_of_points = PD_K;

		_buff.resize(PD_K);
#ifdef PP_NORMALS_OUTPUT
		unsigned _LENGTH = PD_K + PD_n * (PD_hyperplaneIndices.size() + 1);
		_buff.resize(_LENGTH);
#endif
		for (int i = 0; i < PD_K; i++)
			_buff[i] = PD_I[i];
#ifdef PP_NORMALS_OUTPUT
		int j = 0;
		for (int i = 0; i < PD_n; i++)
			_buff[PD_K + j++] = PD_c[i];
		for(int hyperplaneIndex : PD_hyperplaneIndices)
			for (int i = 0; i < PD_n; i++)
				_buff[PD_K + j++] = PD_A[hyperplaneIndex][i];
#endif
		newImage.data = doubleToChar(_buff);
		delete[] PD_I;

		if (PP_RECEPTIVE_FIELD_OUT) {
			_buff.resize(PD_K * PD_n);
			int index = 0;
			for (int i = 0; i < PD_K; i++) {
				for (int j = 0; j < PD_n; j++)
					_buff[index++] = PD_field[i][j];
				delete[] PD_field[i];
			}
			newImage.field_points = doubleToChar(_buff);
			delete[] PD_field;
		}

		PD_images.push_back(newImage);
		//if (newImage.id)
		//	cout << "New image ID = " << newImage.id << " is successfully saved." << endl;
	}
#else
		if (fprintf(PD_stream_outFile, "%d;%d;%d", PD_id, PD_currentProblem, PD_currentTrace) == 0)
			cout << "Error writing to " << PD_outFilename << " on problem " << PD_currentProblem << ", trace " << PD_currentTrace << endl;
		for (int i = 0; i < PD_K; i++)
			if (fprintf(PD_stream_outFile, ";%f", PD_I[i]) == 0)
				cout << "Error writing to " << PD_outFilename << " on problem " << PD_currentProblem << ", trace " << PD_currentTrace << ", PD_I index" << i << endl;
		for (int i = 0; i < PD_n; i++)
			if (fprintf(PD_stream_outFile, ";%f", PD_answerVector[i]) == 0)
				cout << "Error writing to " << PD_outFilename << " on problem " << PD_currentProblem << ", trace " << PD_currentTrace << ", PD_answerVector index" << i << endl;
		for (int i = 0; i < PD_n; i++)
			if (fprintf(PD_stream_outFile, ";%f", PD_cosVector[i]) == 0)
				cout << "Error writing to " << PD_outFilename << " on problem " << PD_currentProblem << ", trace " << PD_currentTrace << ", PD_cosVector index" << i << endl;
		fprintf(PD_stream_outFile, "\n");
#ifdef PP_BSF_ITER_OUTPUT
	cout << "End of writing to " << PD_outFilename << endl;
#endif // PP_BSF_ITER_OUTPUT
	}
	PD_id++;

	if (PP_RECEPTIVE_FIELD_OUT) {
		// Outpur retinas
		if (fprintf(PD_stream_retFile, "%d;%d;%d", PD_id - 1, PD_K, PD_n) == 0)
			cout << "Error writing to " << PD_retFilename << " on problem " << PD_currentProblem << ", trace " << PD_currentTrace << endl;
		for (int i = 0; i < PD_K; i++) {
			for (int j = 0; j < PD_n; j++)
				if (fprintf(PD_stream_retFile, ";%.14f", PD_field[i][j]) == 0)
					cout << "Error writing to " << PD_retFilename << " on problem " << PD_currentProblem << ", trace " << PD_currentTrace << ", PD_I index" << i << endl;
		}
		fprintf(PD_stream_retFile, "\n");
#ifdef PP_BSF_ITER_OUTPUT
		cout << "End of writing to " << PD_retFilename << endl;
#endif // PP_BSF_ITER_OUTPUT
	}
#endif //PP_DATABASE_OUTPUT
}

// 1. Movement on Polytope
void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {

}

void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// not used
}

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// not used
}

//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; };
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; };
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; };
void PC_bsfAssignMpiMaster(int value) { BSF_sv_mpiMaster = value; };
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; };
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; };
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; };
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; };

//----------------------------- User functions -----------------------------
inline void basis_Init() {
	//PD_c
	int j;
	PT_float_T length;
	PT_float_T tailSum;
	PT_vector_T PD_c2;
	copy_Vector(PD_c2, PD_c);
	square_Vector(PD_c2);

	copy_Vector(PD_E[0], PD_c);
	for (int i = 1; i < PD_n; i++) {
		for (j = 0; j < i; j++)	PD_E[i][j] = 0;
		tailSum = vector_Sum(PD_c2, i);
		if (tailSum == 0) {
			PD_E[i][i - 1] = 0;
			PD_E[i][i] = 1;
			j++;
			for (; j < PD_n; j++) { PD_E[i][j] = 0; }
		}
		else if (PD_c[i - 1] == 0.) {
			PD_E[i][i - 1] = 1.;
			for (; j < PD_n; j++) { PD_E[i][j] = 0; }
		}
		else {
			PD_E[i][i - 1] = (PT_float_T)((-1. * tailSum) / PD_c[i - 1]);
			for (; j < PD_n; j++) { PD_E[i][j] = PD_c[j]; }
		}
		length = sqrt(dotproduct_Vector(PD_E[i], PD_E[i]));
		divide_Vector(PD_E[i], length);
	}
}
inline void print_Vector(PT_vector_T x) {
	for (int i = 0; i < PD_n; i++)
		cout << x[i] << " ";
}
inline void add_Vector(PT_vector_T To, PT_vector_T From) {
	for (int i = 0; i < PD_n; i++)
		To[i] += From[i];
}
inline void subtract_Vector(PT_vector_T To, PT_vector_T From) {
	for (int i = 0; i < PD_n; i++)
		To[i] -= From[i];
}
inline void copy_Vector(PT_vector_T To, PT_vector_T From) {
	for (int i = 0; i < PD_n; i++)
		To[i] = From[i];
}
inline void multiply_Vector(PT_vector_T To, PT_float_T C) {
	for (int i = 0; i < PD_n; i++)
		To[i] *= C;
}
inline PT_float_T dotproduct_Vector(PT_vector_T x, PT_vector_T y) {
	PT_float_T result = 0.0f;
	for (int i = 0; i < PD_n; i++) {
		result += x[i] * y[i];
	}
	return result;
}
inline void divide_Vector(PT_vector_T To, PT_float_T C) {
	for (int i = 0; i < PD_n; i++)
		To[i] /= C;
}
inline void square_Vector(PT_vector_T To) {
	for (int i = 0; i < PD_n; i++)
		To[i] *= To[i];
}
inline PT_float_T vector_Sum(PT_vector_T v, int start) {
	PT_float_T result = 0.0f;
	for (int i = start; i < PD_n; i++) {
		result += v[i];
	}
	return result;
}
inline void norm_Vector(PT_vector_T To) {
	PT_float_T length = 0;
	for (int i = 0; i < PD_n; i++)
		length += To[i] * To[i];
	length = sqrt(length);
	for (int i = 0; i < PD_n; i++)
		To[i] /= length;
}
inline void basis_Print() {
	for (int i = 1; i < PD_n; i++) {
		print_Vector(PD_E[i]);
		cout << endl;
	}
}

inline void G(PT_integer_T k, PT_vector_T out, bool cross) {
	PT_vector_T tempPoint;
	PT_vector_T coordinate;
	PT_integer_T dimensionPointsNumber;
	int i[PP_MAX_N], currentDimension;
	PT_integer_T pointNo = k;
	
	if (!cross) {
		for (int j = PD_n - 1; j > 0; j--) {
			dimensionPointsNumber = (PT_integer_T)powf(2 * PP_ETA + 1, (PT_float_T)j - 1); //Possible overfilling!
			i[j - 1] = (int)(pointNo / dimensionPointsNumber);
			pointNo = pointNo % dimensionPointsNumber;
		}
	}
	else
	{
		dimensionPointsNumber = 2 * PP_ETA;
		currentDimension = (int)(pointNo / dimensionPointsNumber);
		pointNo = (pointNo % dimensionPointsNumber);
		for (int j = 0; j < PD_n - 1; j++) {
			if (j != currentDimension)
				i[j] = 0;
			else
				if (pointNo < PP_ETA)
					i[j] = pointNo - PP_ETA;
				else
					i[j] = pointNo - PP_ETA + 1;
		}
	}
	copy_Vector(tempPoint, PD_z);
	for (int j = 1; j < PD_n; j++) {
		copy_Vector(coordinate, PD_E[j]);
		if(!cross)
			multiply_Vector(coordinate, (PT_float_T)(i[j - 1] * PP_DELTA - PP_ETA * PP_DELTA));
		else
			multiply_Vector(coordinate, (PT_float_T)(i[j - 1] * PP_DELTA));
		add_Vector(tempPoint, coordinate);
	}
	for (int i = 0; i < PD_n; i++)
		out[i] = tempPoint[i];
};



inline bool isInnerPoint(PT_vector_T point) {
	bool result = true;
	for (int i = 0; i < PD_m; i++) {
		if (dotproduct_Vector(PD_A[i], point) > PD_b[i])
			result = false;
	}
	return result;
}

inline void targetProjection(int i, PT_vector_T _In, PT_vector_T _Out) {
	PT_vector_T projection;
	PT_vector_T temp;

	//------------ Computing target projection gamma_i ----------//
	copy_Vector(temp, PD_c);
	multiply_Vector(temp, (PT_float_T)((dotproduct_Vector(PD_A[i], _In) - PD_b[i]) / dotproduct_Vector(PD_A[i], PD_c)));
	copy_Vector(projection, _In);
	subtract_Vector(projection, temp);
	copy_Vector(_Out, projection);
}

inline void fieldProjection
(PT_vector_T _In, PT_vector_T _Out) {
	PT_vector_T projection;
	PT_vector_T temp;
	PT_vector_T distance;
	
	//------------ Computing target projection pi_c ----------//
	copy_Vector(temp, _In);
	subtract_Vector(temp, PD_z);
	copy_Vector(distance, PD_c);
	multiply_Vector(distance, (PT_float_T)((dotproduct_Vector(PD_c, temp)) / dotproduct_Vector(PD_c, PD_c)));
	copy_Vector(projection, _In);
	subtract_Vector(projection, distance);
	copy_Vector(_Out, projection);
}

inline void coordinateAngles(PT_vector_T _In, PT_vector_T _Out)
{
	PT_vector_T field_vector;
	PT_float_T length = 0;
	copy_Vector(field_vector, _In);
	subtract_Vector(field_vector, PD_z);

	//------------ Computing coordinate angles cos alpha_i ----------//
	for (int i = 0; i < PD_n; i++)
		length += field_vector[i] * field_vector[i];
	length = sqrt(length);

	for(int i = 0; i < PD_n; i++)
		_Out[i] = dotproduct_Vector(PD_E[i], field_vector) / length;
}

inline PT_float_T objectiveDistance(PT_vector_T x) {
	PT_vector_T temp;

	//------------ Computing target distance rho_c -------------//
	copy_Vector(temp, PD_z);
	subtract_Vector(temp, x);

	return (PT_float_T)(dotproduct_Vector(PD_c, temp) / sqrt(dotproduct_Vector(PD_c, PD_c)));
}

inline PT_float_T bias(int i) {
	PT_float_T result = (PT_float_T)((dotproduct_Vector(PD_A[i], PD_g) - PD_b[i]) / dotproduct_Vector(PD_A[i], PD_c)) * sqrt(dotproduct_Vector(PD_c, PD_c));
	return result;
}

#ifdef PP_DATABASE_OUTPUT
std::vector<double> charToDouble(std::vector<char> _In) {
	std::vector<double> _Out;
	size_t size = _In.size();
	_Out.resize(size / sizeof(double));
	std::memcpy(_Out.data(), _In.data(), size);
	return _Out;
}

std::vector<char> doubleToChar(std::vector<double> _In) {
	std::vector<char> _Out;
	size_t size = _In.size() * sizeof(double);
	_Out.resize(size);
	std::memcpy(_Out.data(), _In.data(), size);
	return _Out;
}

void printLppForm(Problem problem, std::vector<Inequality> inequalities) {
	int M = 2 * problem.N + inequalities.size();
	int width = 10;
	std::vector<double> _vec;
	std::cout << problem.id << '\t' << problem.N << '\t' << M << std::endl;
	for (int i = 0; i < problem.N; i++) {
		for (int j = 0; j < problem.N; j++)
			if (i == j) std::cout << std::setw(width) << double(1);
			else std::cout << std::setw(width) << double(0);
		std::cout << std::setw(width) << problem.high << std::endl;
	}
	for (int i = 0; i < inequalities.size(); i++) {
		_vec = charToDouble(inequalities[i].coefficients);
		for (int j = 0; j < problem.N; j++)
			std::cout << std::setw(width) << _vec[j];
		std::cout << std::setw(width) << inequalities[i].b << std::endl;
	}
	for (int i = 0; i < problem.N; i++) {
		for (int j = 0; j < problem.N; j++)
			if (i == j) std::cout << std::setw(width) << double(-1);
			else std::cout << std::setw(width) << double(0);
		std::cout << std::setw(width) << problem.low << std::endl;
	}
	_vec = charToDouble(problem.c);
	std::copy(_vec.begin(), _vec.end(), std::ostream_iterator<double>(std::cout, "\t"));
	std::cout << std::endl;
}
#endif // PP_DATABASE_OUTPUT
