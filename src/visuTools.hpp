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

#ifndef VISU_TOOLS_HPP
#define VISU_TOOLS_HPP

#include "Node.hpp"
#include "Particles.hpp"

// Vectors updating
void rechargeOctreeBB(Node<Particles> * n, int & height, vector <vector <double> > & vBB);
void rechargeOctreeParticles(Node<Particles> * n, int & height, vector<int> & vNbParticles, 
	vector<int> & vIdxParticles);

// Morton computation	
void computeMortonLineBITS(vector<vector<double>> vBB, vector<vec3D> & morton);
void computeMortonLineDFS(Node<Particles> * treeOCTREE, vector<vec3D> & MortonLineDFS, const int & level);
void readCentersFromFile(vector<vec3D> & centersline, const string & file);
void computeBoundingBoxesFromCentersByLevel(vector <vector<vec3D>> centersPerLevel, vector <vector< vector<double> > > & vBB, int nbLevels, double maxEdge);
void updateBB_BFS(Node<Particles> * treeMORTON, vector <vector< vector<double>>> & vBB, int level);


void updateMortonLineBITS(vector<double> bb, vector<vec3D> & morton);
ui64 splitBy3(unsigned int a);
ui64 computeMorton3DKey(vec3D coord);
bool mortonInf(vec3D c1, vec3D c2);
void sortMortonLineBITSVector(vector<vec3D> & morton);	


void updateMortonBB(Node<Particles> * n, vector <vector <double> > & vBBM);

void printInfoTerminal(int state, int level);

void loadCoordinatesASCIIWithoutQuantity(int & nbParticles, vec3D *& coords, const string & file);
void printKeyboardInfos();
void printMenuInfos();

#endif
