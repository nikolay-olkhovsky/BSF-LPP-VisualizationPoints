/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Visualization of Linear Programming Problem (ViLiPP)
Module: Problem-Data.h (Problem Data)
Prefix: PD
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "Problem-Types.h"			// Problem Parameters 
using namespace std;
//========================== Problem variables ====================================
static int				PD_n;		// Space dimension
static int				PD_m;		// Number of inequalities
static PT_integer_T		PD_K;		// Number of receptive fied points
static PT_integer_T		PD_id;		// Number of current precedent
static PT_float_T		PD_time;	// Time counter
static PT_integer_T		PD_problemsNumber;
static PT_integer_T		PD_currentProblem;
static PT_integer_T		PD_tracesNumber;
static PT_integer_T		PD_currentTrace;
static PT_integer_T		PD_testSum;
static PT_integer_T		PD_read_state;	// State of problem trace processing
static bool				PD_initState;	// Flag to initial state when parameter is not set
//========================== Problem data structures ==============================
static PT_matrix_T PD_A;	// Matrix of the system Ax <= b
static PT_column_T PD_b;	// Column of the constant terms of the system Ax <= b
static PT_vector_T PD_c;	// Coefficients of the objective function <c,x>
static PT_matrix_T PD_E;	// Matrix of vectors e(i) forming basis othogonal to objective function
static PT_vector_T PD_g;	// Point of retina
static PT_vector_T PD_z;	// Center of retina
static vector<int> PD_hyperplaneIndices; // Indices of active hyperplanes
static PT_vector_T PD_nextPoint;	// Next point in the trace
static PT_vector_T PD_answerVector;	// Answer for the precedent
static PT_vector_T PD_fieldVector;	// Answer for the precedent
static PT_vector_T PD_cosVector;	// Answer for the precedent
static PT_image_T PD_I; // Retina
static PT_field_T PD_field;
//========================== Files ================================================
static string PD_problemName;
#ifdef PP_DATABASE_OUTPUT
static auto storage = sqlite_orm::make_storage("C:/HS/dataset100000_density0001.sqlite3",
    sqlite_orm::make_table("problems",
        sqlite_orm::make_column("id", &Problem::id, sqlite_orm::primary_key()),
        sqlite_orm::make_column("N", &Problem::N),
        sqlite_orm::make_column("seed", &Problem::seed),
        sqlite_orm::make_column("high", &Problem::high),
        sqlite_orm::make_column("low", &Problem::low),
        sqlite_orm::make_column("c", &Problem::c)
    ),
    sqlite_orm::make_table("inequalities",
        sqlite_orm::make_column("id", &Inequality::id, sqlite_orm::primary_key()),
        sqlite_orm::make_column("coefficients", &Inequality::coefficients),
        sqlite_orm::make_column("b", &Inequality::b),
        sqlite_orm::make_column("problem_id", &Inequality::problem_id),
        sqlite_orm::foreign_key(&Inequality::problem_id).references(&Problem::id)
    ),
    sqlite_orm::make_table("surface_points",
        sqlite_orm::make_column("id", &SurfacePoint::id, sqlite_orm::primary_key()),
        sqlite_orm::make_column("coefficients", &SurfacePoint::coefficients),
        sqlite_orm::make_column("problem_id", &SurfacePoint::problem_id),
        sqlite_orm::foreign_key(&SurfacePoint::problem_id).references(&Problem::id)
    ),
    sqlite_orm::make_table("precedents",
        sqlite_orm::make_column("id", &Precedent::id, sqlite_orm::primary_key().autoincrement()),
        sqlite_orm::make_column("coefficients", &Precedent::coefficients),
        sqlite_orm::make_column("face", &Precedent::face),
        sqlite_orm::make_column("d", &Precedent::d),
        sqlite_orm::make_column("shift", &Precedent::shift),
        sqlite_orm::make_column("face_count", &Precedent::face_count),
        sqlite_orm::make_column("face_numbers", &Precedent::face_numbers),
        sqlite_orm::make_column("problem_id", &Precedent::problem_id),
        sqlite_orm::foreign_key(&Precedent::problem_id).references(&Problem::id)
    ),
    sqlite_orm::make_table("images",
        sqlite_orm::make_column("id", &Image::id, sqlite_orm::primary_key().autoincrement()),
        sqlite_orm::make_column("density", &Image::density),
        sqlite_orm::make_column("rank", &Image::rank),
        sqlite_orm::make_column("answer_vector", &Image::answer_vector),
        sqlite_orm::make_column("cosine_vector", &Image::cosine_vector),
        sqlite_orm::make_column("num_of_points", &Image::num_of_points),
        sqlite_orm::make_column("data", &Image::data),
        sqlite_orm::make_column("precedent_id", &Image::precedent_id),
        sqlite_orm::make_column("field_points", &Image::field_points),
        sqlite_orm::foreign_key(&Image::precedent_id).references(&Precedent::id)
    )
);
static std::vector<unsigned> PD_problemIds;
static std::vector<unsigned> PD_traceIds;
std::vector<Problem> PD_problems;
std::vector<Inequality> PD_inequalities;
std::vector<Precedent> PD_points;
std::vector<Precedent> PD_traces;
std::vector<Image> PD_images;
#else
static FILE* PD_stream_lppFile;
static string PD_lppFilename; /* LPP file in the following format:
------------ begin of file -------------
numOfProblems
m n
A_11 A_12 ... A_1n b_1
A_21 A_22 ... A_2n b_2
...
A_m1 A_m2 ... A_mn b_m
c_1 c_2 ... c_n
...
------------ end of file----------------*/
static FILE* PD_stream_traceFile;
static string PD_traceFilename; /* TRACE file in the following format:
------------ begin of file -------------
numOfTraces
m n
A_11 A_12 ... A_1n b_1
A_21 A_22 ... A_2n b_2
...
A_m1 A_m2 ... A_mn b_m
...
------------ end of file----------------*/
static FILE* PD_stream_outFile;
static string PD_outFilename;

static FILE* PD_stream_retFile;
static string PD_retFilename; /* OUT file in the following format:
------------ begin of file -------------
numOfRetinas
m n
R_1 R_2 ... R_n
...
------------ end of file----------------*/
#endif // PP_DATABASE_OUTPUT
