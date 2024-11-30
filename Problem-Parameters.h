/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Visualization of Linear Programming Problem (ViLiPP)
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author: Nikolay A. Olkhovsky

This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#define PP_FILE_INI "config.ini"
//=========================== Problem Parameters =========================
static int PP_N;
#define PP_DATABASE_OUTPUT
//#define PP_NORMALS_OUTPUT
#define PP_MAX_M	61								// Maximal number of inequalities (PP_N + 1 + PP_NUM_OF_RND_INEQUALITIES)
#define PP_MAX_N	81									// Maximal  space dimension equals (2*PP_N + 1 + PP_NUM_OF_RND_INEQUALITIES)
#define PP_MAX_K	1000000000
//#define PP_ETA		30									// Rank of receptive field 		
//#define PP_DELTA	0.2									// Density of receptive field
static int			PP_ETA;
static double		PP_DELTA;
static std::string	PP_PATH;
static std::string	PP_LPP_FILE;
//#define PP_TRACE_FILE "lp_dataset_tr.mtx.pack"
//#define PP_RETINA_FILE "lp_dataset_ret.csv" // File with retina generated for every point
//#define PP_IMAGE_FILE "lp_dataset_im.csv"		// File with output results
static bool PP_IMAGE_OUT;					// Flag to save list of target distances to the file
static bool PP_RECEPTIVE_FIELD_OUT;			// Flag to save receptive field points to the file
static bool PP_CROSSFILED;			// Form of the receptive field (false - hypercube, true - cross)

//------------------------- Matrix format ----------------
static std::string PP_PROBLEM_NAME;
static std::string PP_MTX_PREFIX;
static std::string PP_MTX_POSTFIX_TR;
static std::string PP_MTX_POSTFIX_RET;
static std::string PP_MTX_POSTFIX_IM;

//-------------------------- Macroses ---------------------------
#define PF_MIN(x,y) (x<y?x:y)
#define PF_MAX(x,y) (x>y?x:y)

//-------------------------- States ---------------------------
#define PP_STATE_NEW_PROBLEM		0
#define PP_STATE_NEW_POINT			1
