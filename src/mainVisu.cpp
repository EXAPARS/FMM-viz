/*
  Copyright 2015 - UVSQ
  Authors list: Nathalie Möller, Eric Petit

  This file is part of the FMM-viz.

  FMM-viz is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later version.

  FMM-viz is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License along with
  the FMM-viz. If not, see <http://www.gnu.org/licenses/>.
*/

#include "mpi.h"
#include "GASPI.h"
#include <GL/glut.h>
#include <unistd.h>


#include "Node.hpp"
#include "Particles.hpp"
#include "Decomposition.hpp"
#include "Gaspi_communicator.hpp"

#include <iostream>
#include <iomanip>
#include <string.h>
#include <algorithm>
#include <tuple>

#include <limits.h>
#include <unistd.h>

#include "LBHistApprox.hpp"
#include "LBHistExact.hpp"
#include "LBMortonSyncMPI.hpp"	

#include "Camera.hpp"
#include "visuTools.hpp"
#include "LB.hpp"

// Program settings
int freeze = -1;
int level = 0;		
int state = 0;
int targetHeight;
int info = 1;
int neighbors = 1;
int bb = 1;

// MPI
int nbParticles;

// Trees
Node<Particles> * treeORIGIN = nullptr;
Node<Particles> * treeAPPROX = nullptr;
Node<Particles> * treeEXACT = nullptr;
Node<Particles> * treeOCTREE = nullptr;
Node<Particles> * treeMORTON = nullptr;
Node<Particles> * treeMORTONasync = nullptr;
Node<Particles> * treeMORTONSync = nullptr;

// Load Balancing
int wrank;
int wsize;
LB lb;
decompo nb1ers;
double dist;

// Camera
Camera camera;
bool btState32 = false;
bool btState128 = false;
bool btState256 = false;

// OpenGL variables
float rotY = 0.0f;
float rotX = 0.0f;
float rotZ = 0.0f;
float dx = 0.0f;
float dy = 0.0f;
float dz = 0.0f;
float ratio;
float fovy = 70.0f;

// Bounding Boxes
vector< vector<double> > vBB;			// bounding boxes
vector< vector<double> > vBBMORTON;
vector< vector<double> > vBBMORTONasync;
vector< vector<double> > vBBMORTONSync;
vector <vector< vector<double> > >vBBSpectreMorton; // Bounding boxes by octree level
vector <vector< vector<double> > >vBBMortonByLevel;

// Morton lines
vector<vec3D> MortonLineBITS; 						// Morton
vector<vec3D> MortonLineDFS; 						// Morton
vector<vec3D> MortonLineOnMortonTreeSync; 			// Morton
vector<vec3D> MortonLineOnMortonTreeAsync;
vector<vec3D> centersLine;
vector < vector <vec3D> > centersLinePerLevel;

// Octree
vector<vec3D> vCoordsOctree;	// Coordinates
vector<int> vNbParticlesOCTREE;			// nbParticles
vector<int> vIdxParticlesOCTREE;		// indexes

// Coordinates
vector<vec3D> vCoordsEXACT;			// Exact
vector<vec3D> vCoordsAPPROX;			// Approx
vector<vec3D> vCoordsORIGIN;			// Original
vector<vec3D> vCoordsMORTONSync;			// Morton Mpi
vector<vec3D> vCoordsMORTONasync;	// Morton Gaspi
vector<vec3D> vCoordsSPECTRE;		// Spectre
vector<vec3D> vCoordsSPECTREneighbors;		// Spectre

// spectre LB Morton 
vector<vec3D> vCoordsSPECTREmortoned;		// Spectre
vector<vec3D> vCoordsSPECTREneighborsMortoned;		// Spectre

// spectre LB HIST 
vector<vec3D> vCoordsSPECTREhist;
vector<vec3D> vCoordsSPECTREhistNeighb;


/* ------------------  BOUNDING BOX  ------------------------------ */
float box[24] = {
	0, 				    0, 				    0,
	float(COORDMAX), 	0, 				    0,
	float(COORDMAX), 	float(COORDMAX), 	0,
	0, 				    float(COORDMAX), 	0,
	0, 				    0, 				    float(COORDMAX),
	float(COORDMAX), 	0, 				    float(COORDMAX),
	float(COORDMAX), 	float(COORDMAX), 	float(COORDMAX),
	0, 				    float(COORDMAX), 	float(COORDMAX)		
};

unsigned char boxIndexes[] = {0, 1, 1, 2, 2, 3, 3, 0, 0, 4, 1, 5, 2, 6, 3, 7, 4, 5, 5, 6, 6, 7, 7, 4 };


/* ------------------ UPDATE DATA FUNCTIONS ----------------------- */

void rechargeOctree()
{
	vBB.clear();
	vNbParticlesOCTREE.clear();
	vIdxParticlesOCTREE.clear();
	rechargeOctreeParticles(treeOCTREE, level, vNbParticlesOCTREE, vIdxParticlesOCTREE);	
	rechargeOctreeBB(treeOCTREE, level, vBB);	
}

void rechargeMortonLines()
{
	// clear
	MortonLineBITS.clear();
	MortonLineDFS.clear();
	MortonLineOnMortonTreeSync.clear();
	MortonLineOnMortonTreeAsync.clear();
		
	// update
	computeMortonLineBITS(vBB, MortonLineBITS);
	computeMortonLineDFS(treeOCTREE, MortonLineDFS, level);
	computeMortonLineDFS(treeMORTON, MortonLineOnMortonTreeSync, level);	
	computeMortonLineDFS(treeMORTONasync, MortonLineOnMortonTreeAsync, level);		
}

void rechargeGlobalTabs()
{
	rechargeOctree();
	rechargeMortonLines();
}

/* ------------------  DISPLAY FUNCTIONS  ------------------ */
void displayVectorOfParticles(vector<vec3D> v, int r, int g, int b)
{
	glPointSize(2.0f);	
	glColor3f(r,g,b);	
	glVertexPointer( 3, GL_DOUBLE, 0, v.data() ) ;
	glDrawArrays(GL_POINTS, 0, v.size());	
}

void displayBoundingBoxes(vector< vector<double> >  v)
{	
	
	glLineWidth(1);
	glColor3f(0.5,0.5,0.5);	
	for (decltype(vBB.size()) i=0; i< v.size(); ++i)
	{
		glVertexPointer( 3, GL_DOUBLE, 0, v[i].data() );
		glDrawElements( GL_LINES, 24, GL_UNSIGNED_BYTE, boxIndexes );		
	}
}

void displayBoxesParticles()
{	
	glColor3f(1,1,1);
	for (decltype(vBB.size()) i=0; i< vBB.size(); ++i)
	{
		switch (i%8)
		{
			case 0: glColor3f(1,0.5,0); break;
			case 1: glColor3f(1,1,1); break;
			case 2: glColor3f(1,1,0); break;
			case 3: glColor3f(1,0,1); break;
			case 4: glColor3f(1,0,0); break;
			case 5: glColor3f(0,1,1); break;
			case 6: glColor3f(0,1,0); break;
			case 7: glColor3f(0,0,1); break;
			default: cerr << "Unexpected case" << endl; exit(0);			
		}		
		glVertexPointer( 3, GL_DOUBLE, 0, vCoordsOctree.data() ) ;
		glDrawArrays(GL_POINTS, vIdxParticlesOCTREE[i], vNbParticlesOCTREE[i]);
	}	
}

void displayLine(vector<vec3D> line, int width, int r, int g, int b)
{
	glLineWidth(width);
	glColor3f(r,g,b);			
	glVertexPointer( 3, GL_DOUBLE, 0, line.data() );
	glDrawArrays(GL_LINE_STRIP, 0, line.size());	
}

void displayTree()
{
	switch(state) 
	{
		/** Construction Octree + Exemple Morton **/
		case 0 : // Octree
			displayBoundingBoxes(vBB);
			break;

		case 1 : // Theory - Morton - Bits Interleaving
			displayBoundingBoxes(vBB);
			displayLine(MortonLineBITS, 5,1,0,0);
			break;

		case 2 : // Display UAV
			displayBoxesParticles();	
			break;	
		
		/** Construction Octree + LB avec UAV importé **/			
		case 3 : // EXACT HISTOGRAM			---> ok
			displayVectorOfParticles(vCoordsEXACT,0,1,0);
			displayBoundingBoxes(vBB);
			break;
		
		case 4 : // APPROX HISTOGRAM		---> ok
			displayVectorOfParticles(vCoordsAPPROX,1,0,0);
			displayBoundingBoxes(vBB);
			break;
		
		case 5 : // MORTON - GASPI - Async                      --- Avec création récursive de l'octree ET Morton DFS +++>> OK		
			displayVectorOfParticles(vCoordsMORTONasync,1,0,0);
			displayBoundingBoxes(vBBMORTONasync);
			displayLine(MortonLineOnMortonTreeAsync,5,1,0,0);
			break;
		
		/** Application au code SPECTRE, avec octree de SPECTRE **/						
		case 6 : // Original LB SPECTRE
			displayVectorOfParticles(vCoordsSPECTRE,0,0,0);
			if (neighbors > 0)
			{
				glPointSize(4.0f);
				displayVectorOfParticles(vCoordsSPECTREneighbors,1,0,0);
				glPointSize(2.0f);
			}
			displayBoundingBoxes(vBBSpectreMorton[level]);
			break;		
		
		case 7 : // Spectre + LIB Morton 
			displayVectorOfParticles(vCoordsSPECTREmortoned,0,0,0);
			if (neighbors > 0)
			{
				glPointSize(4.0f);
				displayVectorOfParticles(vCoordsSPECTREneighborsMortoned,1,0,0);
				glPointSize(2.0f);
			}
			displayBoundingBoxes(vBBSpectreMorton[level]);
			break;

		case 8 : // Spectre + Approx Histogram on spectre octree
			displayVectorOfParticles(vCoordsSPECTREhist,0,0,0);
			if (neighbors > 0)
			{
				glPointSize(4.0f);
				displayVectorOfParticles(vCoordsSPECTREhistNeighb,1,0,0);
				glPointSize(2.0f);
			}
			displayBoundingBoxes(vBBSpectreMorton[level]);			
			break;

		default :
			cout << "default case, state : " << state << endl;
	}
}

/* ------------------  DISPLAY FUNCTION---------------------------- */
void displayFunc()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	
	glLoadIdentity();
	gluLookAt( camera.getPosition().x, camera.getPosition().y, camera.getPosition().z,
	     camera.getPosition().x - camera.getDirection().x, camera.getPosition().y - camera.getDirection().y, camera.getPosition().z - camera.getDirection().z,
	     camera.getUp().x, camera.getUp().y, camera.getUp().z
	     );	

	glRotatef(rotY, 0, 1, 0);
	glRotatef(rotX, 1, 0, 0);
	glRotatef(rotZ, 0, 0, 1);
	
	glTranslatef(- 0.000000002 * COORDMAX/2, - 0.000000002 * COORDMAX/2, - 0.000000002 * COORDMAX/2);
	glScalef(0.000000002, 0.000000002, 0.000000002);

	glPointSize(2.0f);
	displayTree();
	
	glutSwapBuffers();
}

void reshapeFunc( int w, int h )
{
  glViewport( 0, 0, w, h );
  ratio = (float)w/(float)h;  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(fovy, ratio, 0.1f, 1000.0f);
  glMatrixMode(GL_MODELVIEW);
}

/* ------------------  INTERACT FUNCTIONS ------------------------- */
void Idle()
{	
	if (freeze != -1)
	{
		rotY += 0.05;
		glutPostRedisplay();
	}
}

void updateState(int oldState, int newState, int oLevel)
{
	// test if Morton
	int willBeMorton;
	if (newState == 3 || newState == 4 || newState == 6 ||newState == 7)
		willBeMorton = 1;
	else
		willBeMorton = 0;
	
	// update state			
	state = newState;
	if (oLevel != 0)
		level = 0;
	
	// recharge necessary data
	if (willBeMorton)
		rechargeMortonLines();
	else 
		rechargeOctree();		
}

void keyboardFunc(unsigned char key, int x, int y) 
{

	switch(key) 
	{
		// rotation		
		case 'a' : rotX += 1; break;
		case 'z' : rotX -= 1; break;	  
		case 'q' : rotY += 1; break;
		case 's' : rotY -= 1; break;
		case 'w' : rotZ += 1; break;
		case 'x' : rotZ -= 1; break;
		case 'i' : rotX = rotY = rotZ = dx = dy = dz = 0; break;
		case ' ' : freeze *= -1; break;	
		
		// neighbors
		case 'n' : neighbors *= -1; break;				
		case 'b' : bb *= -1; break;
		
		// tree level
		case '+'  : if (level < 8)
					{
						level++;					   
						rechargeGlobalTabs();
						cout << "level : " << level << endl;				   					   
				   }
				   break;
		case '-'  : if (level > 0)
				   {
					   level--; 
					   rechargeGlobalTabs();	 
					   cout << "level : " << level << endl;
				   }
				   break; 
		
		// display settings
		case '0' : updateState(state, 0, level); printInfoTerminal(state, level); break;	
		case '1' : updateState(state, 1, level); printInfoTerminal(state, level); break;	
		case '2' : updateState(state, 2, level); printInfoTerminal(state, level); break;
		case '3' : updateState(state, 3, level); printInfoTerminal(state, level); break;		
		case '4' : updateState(state, 4, level); printInfoTerminal(state, level); break;
		case '5' : updateState(state, 5, level); printInfoTerminal(state, level); break;
		case '6' : updateState(state, 6, level); printInfoTerminal(state, level); break;
		case '7' : updateState(state, 7, level); printInfoTerminal(state, level); break;
		case '8' : updateState(state, 8, level); printInfoTerminal(state, level); break;
		case '9' : updateState(state, 9, level); printInfoTerminal(state, level); break; 

		
		// quit
		case 0x1B: exit(0); break;
	}
			
	glutPostRedisplay();
}

/* classes : camera,
*/
void specialFunc(int key, int x, int y)
{
	switch( key ) {
		case GLUT_KEY_LEFT:
			camera.move( Camera::LEFT );
			break;
		case GLUT_KEY_RIGHT:
			camera.move( Camera::RIGHT );
			break;
		case GLUT_KEY_UP:
			camera.move( Camera::FORWARD );
			break;
		case GLUT_KEY_DOWN:
			camera.move( Camera::BACK );
			break;
		case GLUT_KEY_PAGE_UP:
			camera.move( Camera::UP );
			break;
		case GLUT_KEY_PAGE_DOWN:
			camera.move( Camera::DOWN );
			break;
	}
	glutPostRedisplay();	
}

void joystickFunc( unsigned int buttons, int x, int y, int z)
{	
	float scale = 0.02f;
	camera.setStep(scale);

	int X = 0;
	int Y = 0;
	
	// Camera translations
	if (x == -1000)
		camera.move (Camera::LEFT);	  
	if (x == 1000)
		camera.move (Camera::RIGHT);
	if (y == -1000)
		camera.move (Camera::UP);
	if (y == 1000)
		camera.move (Camera::DOWN);
	if( buttons & 16 )
		camera.move (Camera::BACK);	  
	if( buttons & 64 )
		camera.move (Camera::FORWARD);
	
	// Camera rotations
	if( buttons & 2 ) {
	  X = -1;
	}
	if( buttons & 8 ) {
	  X = 1;
	}
	if( buttons & 1 ) {
	  Y = -1;
	}
	if( buttons & 4 ) {
	  Y = 1;
	}  
	camera.yaw( X *0.5 * 3.14159f/180.0f );
	camera.pitch( Y *0.5 * 3.14159f/180.0f );
	 
	// Reset
	if( buttons & 512 )
	{
		camera.setPosition( glm::vec3( 0.0f, 0.0f, 6.0f ) );
		camera.setDirection( glm::vec3( 0.0f, 0.0f, 1.0f ) );
		rotY = 0;								
	}  
	
	// Model Rotation
	if( buttons & 256 )
		btState256 = true;
				
	if( btState256 && !(buttons & 256) )
	{
		btState256 = false;
			freeze *= -1;
	}						
		 
	// Load Balancing strategy
	if( buttons & 32 )
		btState32 = true;
				
	if( btState32 && !(buttons & 32) )
	{
		btState32 = false;
		if (state < 9)
		{
			updateState(state, state+1, level);
			printInfoTerminal(state, level);
		}
	}
	
	if( buttons & 128 )
		btState128 = true;
				
	if( btState128 && !(buttons & 128) )
	{
		btState128 = false;
		if (state > 0)
		{	
			updateState(state, state-1, level);
			printInfoTerminal(state, level);
		}
	}
		
	glutPostRedisplay();  
}

/* ------------------  MAIN --------------------------------------- */

int main(int argc, char* argv[]) 
{
	/**
	 *  MPI initialization
	 **/	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&wrank);
	MPI_Comm_size(MPI_COMM_WORLD,&wsize);	
			
	/**
	 *  OpenGL initialization
	 **/		
	// glut
	glutInit(&argc, argv);
	int w = 380;
	int h = 380;
	glutInitWindowSize(w, h);

	if (wrank < 3)
		glutInitWindowPosition(500 + (wrank*w), 0);
	else if (wrank < 6)
		glutInitWindowPosition(/*1600 +*/ ((wrank-3)*w), 400);
	else if (wrank < 9)
		glutInitWindowPosition(/*1600 +*/ ((wrank-6)*w), 900);
	else
		glutInitWindowPosition(/*1600 +*/ ((wrank-9)*w), 1350);			
			
	string titre = "rank " + to_string(wrank);
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
	glutCreateWindow(titre.c_str());
	
	// OpenGL primitives
	glutIdleFunc(Idle);	
	glutKeyboardFunc(keyboardFunc);
	glutSpecialFunc(specialFunc);
	glutJoystickFunc(joystickFunc, 5);	
	glutDisplayFunc(displayFunc);
	glutReshapeFunc(reshapeFunc);	

	// OpenGL			
	glEnable( GL_DEPTH_TEST );  
	glEnableClientState( GL_VERTEX_ARRAY );
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);	
	camera.setPosition( glm::vec3( 0.0f, 0.0f, 6.0f ) );	
	
	/**
	 * Program settings
	 **/

	SUCCESS_OR_DIE(gaspi_proc_init(GASPI_BLOCK));
	gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

	string dir = "/path/to/lib/";

	/** input **/
	/* Dronera ORIGINAL in spectre*/	
	string file = 	dir + "Dronera/"
					+ to_string(wsize)
					+"_MPI/elemts_before_rank_" 
					+ to_string(wrank) 
					+ "_oct_1.txt";
	string fileSPECTRE = dir  + "Dronera/"
					+ to_string(wsize)
					+ "_MPI/elemts_after__rank_" 
					+ to_string(wrank)
					+ "_oct_1.txt";					 
	string fileSPECTREneighbors	 =	dir  + "Dronera/"
					+ to_string(wsize)
					+ "_MPI/Neighbors_rank_" 
					+ to_string(wrank) 
					+ "_oct_1.txt";

	/* Dronera via library MORTON in spectre */
	string fileSPECTREMortoned = dir + "Dronera_LIB_Morton/"
					+ to_string(wsize)
					+ "_MPI/elemts_after__rank_" 
					+ to_string(wrank)
					+ "_oct_1.txt";					 
	string fileSPECTREneighborsMortoned	 =	dir + "Dronera_LIB_Morton/"
					+ to_string(wsize)
					+ "_MPI/Neighbors_rank_" 
					+ to_string(wrank) 
					+ "_oct_1.txt";

	/* Dronera via library HIST_APPROX in spectre */
	string fileSPECTREHapprox = dir + "Dronera_LIB_Hist/"
					+ to_string(wsize)
					+ "_MPI/elts_after__rank_" 
					+ to_string(wrank)
					+ "_oct_1.txt";	
									 
	string fileSPECTREHapproxneighb	 =	dir +"Dronera_LIB_Hist/"
					+ to_string(wsize)
					+ "_MPI/elts_after__Neighbors_rank_" 
					+ to_string(wrank) 
					+ "_oct_1.txt";					

	/** Drone UAV parameters **/
	int spectreOctreeLevels = 9;
	double maxEdge = 5234.847260;


	// Load original particles and update nbParticles
	vec3D * coords;	
	loadCoordinatesASCIIWithoutQuantity(nbParticles, coords, file);
						
	/** Load Balancing parameters **/
	double dist = 50;
	double tol = 0.0001;
	nb1ers.create(wsize);
	int first = 0;
	int last = nbParticles -1;
	Gaspi_communicator * gComm = nullptr;
	
	/// 0 - ORIGINAL 
	origin(coords, nbParticles, vCoordsORIGIN);

	/// 1 - EXACT Histograms	
	LBB<HistExact>(coords, nbParticles, vCoordsEXACT, treeEXACT, nb1ers, dist, tol, first, last, gComm);
	
	/// 2 - APPROX Histograms
	LBB<HistApprox>(coords, nbParticles, vCoordsAPPROX, treeAPPROX, nb1ers, dist, tol, first, last, gComm);

	/// Calcul de l'OCTREE (division NM)
	targetHeight = 5;
	LBoctree(coords, targetHeight, nbParticles, treeOCTREE, vCoordsOctree);
	rechargeOctreeParticles(treeOCTREE, level, vNbParticlesOCTREE, vIdxParticlesOCTREE);	
	rechargeOctreeBB(treeOCTREE, level, vBB);

	/// 3 - MORTON LINE - ASYNC - GASPI
	LBmortonAsync(file, nbParticles, vCoordsMORTONasync, treeMORTONasync, nb1ers, dist, tol, first, last, gComm);
	updateMortonBB (treeMORTONasync, vBBMORTONasync);	

	/// 4 - SPECTRE LOAD BALANCING RESULTS
	// Spectre's UAV Octree
	centersLinePerLevel.resize(spectreOctreeLevels);
	for (int i=0; i<spectreOctreeLevels; i++)
		readCentersFromFile(centersLinePerLevel[i],"/home/nath/Documents/Dassault-UVSQ/2_These/FMM/FMM-input-DA/Dronera_LIB_Morton/centers/centers__rank_0_level_" + to_string(i) + ".txt");
	computeBoundingBoxesFromCentersByLevel(centersLinePerLevel, vBBSpectreMorton, spectreOctreeLevels, maxEdge*641017);	
	// Orig Spectre LB
	LBspectre(fileSPECTRE, vCoordsSPECTRE);
	LBspectre(fileSPECTREneighbors, vCoordsSPECTREneighbors);

	/// 5 - Résultats de l'utilisation de LB Morton, sur octree SPECTRE.
	LBspectre(fileSPECTREMortoned, vCoordsSPECTREmortoned);
	LBspectre(fileSPECTREneighborsMortoned, vCoordsSPECTREneighborsMortoned);

	/// 6 - Résultats de l'utilisation de LB Histo APPROX, sur octree SPECTRE. (Exact ne collerait pas à l'octree)
	LBspectreHist(fileSPECTREHapprox, vCoordsSPECTREhist);
	LBspectreHist(fileSPECTREHapproxneighb, vCoordsSPECTREhistNeighb);

	printInfoTerminal(state, level);

	// stop gaspi
	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

	// OpenGL display LOOP
	glutMainLoop();
	
	// dealloc
	delete treeORIGIN;
	delete treeAPPROX;
	delete treeEXACT;
	delete treeOCTREE;	
	
	// End MPI
	MPI_Finalize();
	
	return 0;
}
