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

#include "visuTools.hpp"
#include <queue>
#include <set>

/**
 * Parcours l'arbre en BFS et appel à la fonction compBB pour chaque node
 **/
void rechargeOctreeBB(Node<Particles> * n, int & height, vector <vector <double> > & vBB)
{		
	queue<Node<Particles> *> fifo;	// BFS fifo queue
	int depth, nbChilds;

	// add the current node to the fifo list
	fifo.push(n);
	
	// while fifo list not empty
	while(!fifo.empty())
	{
		// take the first node
		Node<Particles> * current = fifo.front();
		fifo.pop();
		
		// add to the global vectors
		depth = current->getDepth();
		(current->getChildren()) ? nbChilds=current->getNbChildren() : nbChilds=0;

		if( 
			((depth < height) && (!nbChilds)) || 	// only leaves above height
			(depth == height) 						// all nodes at height
		  )
			vBB.push_back(current->compBB());			
				
		// add childs to fifo list
		for (int i=0; i<nbChilds; ++i)			
			fifo.push(current->getChildren()[i]);
	}
}

void rechargeOctreeParticles(Node<Particles> * n, int & height, vector<int> & vNbParticles, 
	vector<int> & vIdxParticles)
{		
	queue<Node<Particles> *> fifo;	// BFS fifo queue
	int depth, nbChilds;

	// add the current node to the fifo list
	fifo.push(n);
	
	// while fifo list not empty
	while(!fifo.empty())
	{
		// take the first node	
		Node<Particles> * current = fifo.front();
		fifo.pop();
		
		// add to the global vectors
		depth = current->getDepth();
		(current->getChildren()) ? nbChilds=current->getNbChildren() : nbChilds=0;

		if( 		
			((depth < height) && (!nbChilds)) || 	// only leaves above height
			(depth == height) 						// all nodes at height			
		  )
		{ 
			vNbParticles.push_back(current->getContent().getNbParticles());
			vIdxParticles.push_back(current->getContent().getFirstIndex());			
		}				
		// add childs to fifo list
		for (int i=0; i<nbChilds; ++i)			
			fifo.push(current->getChildren()[i]);
	}
}


/**
 *  Morton Ordering - Method with BITS INTERLEAVING
 **/ 
 
void computeMortonLineBITS(vector<vector<double>> vBB, vector<vec3D> & morton)
{
	// add all BB centers
	for (unsigned int i=0; i<vBB.size(); i++)
		updateMortonLineBITS(vBB[i], morton);
		
	// sort
	sortMortonLineBITSVector(morton);
} 
 
void updateMortonLineBITS(vector<double> bb, vector<vec3D> & morton)
{
	// computer bb center coordinates
	double Ox = (bb[0] + bb[3] )/2;
	double Oy = (bb[1] + bb[10])/2;
	double Oz = (bb[2] + bb[14])/2;
	
	// insert into the morton coordinates vector
	morton.push_back(vec3D(Ox, Oy, Oz));
}


// method to seperate bits from a given integer 3 positions apart
ui64 splitBy3(unsigned int a)
{
    ui64 x = a & 0x1fffff; // we only look at the first 21 bits
    x = (x | x << 32) & 0x1f00000000ffff;  // shift left 32 bits, OR with self, and 00011111000000000000000000000000000000001111111111111111
    x = (x | x << 16) & 0x1f0000ff0000ff;  // shift left 16 bits, OR with self, and 00011111000000000000000011111111000000000000000011111111
    x = (x | x << 8) & 0x100f00f00f00f00f; // shift left 8 bits, OR with self, and 0001000000001111000000001111000000001111000000001111000000000000
    x = (x | x << 4) & 0x10c30c30c30c30c3; // shift left 4 bits, OR with self, and 0001000011000011000011000011000011000011000011000011000100000000
    x = (x | x << 2) & 0x1249249249249249;    
    return x;
}

ui64 computeMorton3DKey(vec3D coord)
{
	// cast coordinates to usigned int and right shift by 11
	ui32 x = static_cast<ui32>(coord.x) >> 11;
	ui32 y = static_cast<ui32>(coord.y) >> 11;
	ui32 z = static_cast<ui32>(coord.z) >> 11;
    
    // compute the morton key
    ui64 key = 0;
    key |= splitBy3(x) | splitBy3(y) << 1 | splitBy3(z) << 2;
	return key;
}

bool mortonInf(vec3D c1, vec3D c2)
{
	return(computeMorton3DKey(c1) < computeMorton3DKey(c2));
}

void sortMortonLineBITSVector(vector<vec3D> & morton)
{
	sort(morton.begin(), morton.end(), mortonInf);	
}

/**
 *  Morton Line - DEPTH FIRST SEARCH on Octree
 **/ 
 
void computeMortonLineDFS(Node<Particles> * treeOCTREE, vector<vec3D> & MortonLineDFS, const int & level)
{
	if (treeOCTREE != NULL)
	{
		if ( (treeOCTREE->getDepth() < level) && (treeOCTREE->getNbChildren()) )
			for (int i=0; i<(treeOCTREE->getNbChildren()); ++i)
				computeMortonLineDFS((treeOCTREE->getChildren()[i]), MortonLineDFS, level);
		else
		{
			if (treeOCTREE->getNbItems() > 0)
			{
				vector<double> bb = treeOCTREE->compBB();
				updateMortonLineBITS(bb, MortonLineDFS);
			}
		}
	}
}

/**
 * Morton Load Balancing, update structures
 **/

void updateMortonBB(Node<Particles> * n,  vector <vector <double> > & vBBM)
{
	if ( n->getChildren() > 0 )
	{
		// recursive call
		for (int i=0; i<n->getNbChildren(); ++i)
			updateMortonBB(n->getChildren()[i], vBBM);
	}
	else
	{
		if (n->getNbItems() > 0)
			vBBM.push_back(n->compBB());
	}
}
/**
 * Morton line - By reading centers from file
 **/
void readCentersFromFile(vector<vec3D> & centersline, const string & file)
{
	/* FIXME ! les coefficients de scaling en dur correspondent au drone */
	double coeff = 641017;
	double translate = 1436.4;
	
	fstream in;
	in.exceptions(ifstream::failbit|ifstream::badbit);
	try
	{
		in.open(file, ifstream::in);
		
		// get number of nodes
		int nbPoints = std::count(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>(), '\n');			
		centersline.resize(nbPoints); 
		
		// load into _coordinates
		in.seekg(0, ios::beg);
		int index = 0;
		while (index < nbPoints)
		{
			in >> centersline[index].x;
			in >> centersline[index].y;
			in >> centersline[index].z;			
			index++;
		}
		in.close();		
	
		// scaling 
		for (int i=0; i<nbPoints; i++)
		{
			centersline[i].x = (centersline[i].x + translate) * coeff;
			centersline[i].y = (centersline[i].y + translate) * coeff;
			centersline[i].z = (centersline[i].z + translate) * coeff;
		}
	}
	catch(ifstream::failure &e)
	{
		cerr << "[erreur] " << e.what() << endl;
		exit(0);
	}
}

/**
 * Calcule les bounding boxes à partir du vector de lines par niveau
 **/
void computeBoundingBoxesFromCentersByLevel(vector <vector<vec3D>> centersPerLevel, vector <vector< vector<double> > > & vBB, int nbLevels, double maxEdge)
{

	cout << setprecision(4);
//	cout << "maxEdge " << maxEdge << endl;
//	cout << "COORDMAX " << COORDMAX << endl;
	
	// Resize nb Levels of bounding boxes
	vBB.resize(nbLevels);
	double edge;
	
	// For each level
	for (int i=0; i<nbLevels; i++)
	{
		edge = maxEdge / (1 << i);
//		cout << fixed <<"edge " <<  edge << endl;
//		cout << "pow " << (1 << i) << endl;
		// Resize level of bounding boxes : = nbCenters at this level
		vBB[i].resize(centersPerLevel[i].size());
				
		// For each center
		for (uint32_t j=0; j<centersPerLevel[i].size(); j++)
		{
			// corners
			double xMIN = centersPerLevel[i][j].x - (0.5 * edge);
			double xMAX = centersPerLevel[i][j].x + (0.5 * edge);
			double yMIN = centersPerLevel[i][j].y - (0.5 * edge);
			double yMAX = centersPerLevel[i][j].y + (0.5 * edge);
			double zMIN = centersPerLevel[i][j].z - (0.5 * edge);
			double zMAX = centersPerLevel[i][j].z + (0.5 * edge);
			
			// bb
			vector<double> bb;
			bb.resize(24); // 8 x 3D coordinates
			bb[0]  = xMIN;	bb[1]  = yMIN; 	bb[2]  = zMIN;
			bb[3]  = xMAX;	bb[4]  = yMIN; 	bb[5]  = zMIN;
			bb[6]  = xMAX;	bb[7]  = yMAX; 	bb[8]  = zMIN;
			bb[9]  = xMIN;	bb[10] = yMAX; 	bb[11] = zMIN;
			bb[12] = xMIN;	bb[13] = yMIN; 	bb[14] = zMAX;
			bb[15] = xMAX;	bb[16] = yMIN; 	bb[17] = zMAX;
			bb[18] = xMAX;	bb[19] = yMAX; 	bb[20] = zMAX;
			bb[21] = xMIN;	bb[22] = yMAX; 	bb[23] = zMAX;					
		
			// insert bb into vBB
			vBB[i][j]=bb;		
		}		
	}
}

void updateBB_BFS(Node<Particles> * n, vector <vector< vector<double> > > & vBB, int level)
{
	queue<Node<Particles> *> fifo;	// BFS fifo queue

	
	vBB.resize(100);
	
	// add the current node to the fifo list
	fifo.push(n);
	
	// while fifo list not empty
	while(!fifo.empty())
	{
		// take the first node
		Node<Particles> * current = fifo.front();
		fifo.pop();
		
		/** do the job **/
		int level = current->getDepth();				
		//cout << "level :" << level << endl;
		
		// non empty boxes only
		if (current->getNbItems()>0)
		{	
			vBB[level].push_back(current->compBB());
				
			// add childs to fifo list
			int nbChilds = current->getNbChildren();
			if (nbChilds > 0)
				for (int i=0; i<nbChilds; ++i)			
					fifo.push(current->getChildren()[i]);
		}
	}
}



/**
 * Display
 **/

void printKeyboardInfos()
{
	cout << endl << "--- Keyboard Functions ------------------  " << endl;
	cout << "a | z \t\t\t: rotate around X-axis" << endl;
	cout << "q | s \t\t\t: rotate around Y-axis" << endl;
	cout << "w | x \t\t\t: rotate around Z-axis" << endl;
	cout << "up | down \t\t: translate along Z-axis" << endl;
	cout << "left | right \t\t: translate along X -axis" << endl;
	cout << "page up | page down \t: translate along Y -axis" << endl;
	cout << "+ | - \t\t\t: explore levels" << endl;
	cout << "i \t\t\t: restore initial view" << endl;
}

void printMenuInfos()
{

	cout << endl << "--- Display Choices------------------  " << endl;	
	cout << "0 - Octree " << endl;
	cout << "1 - Morton Traversal" << endl;			
	cout << "2 - UAV Test Case presentation" << endl;			
	cout << "3 - Exact Histograms " << endl;
	cout << "4 - Approximate Histograms " << endl;
	cout << "5 - Morton" << endl;
	cout << "6 - Original LB " << endl;
	cout << "7 - Morton " << endl;
	cout << "8 - Histogram " << endl;	
}

void printInfoTerminal(int state, int level)
{

	
	switch(state) 
	{
		/** Bases : Octree + Morton **/
		case 0 : 
			cout << "\n\n\n\n\n\n---------------------" << endl;
			cout << "|   Current State   |" << endl;
			cout << "---------------------" << endl;		
			cout << "0 - Octree " << endl;
			cout << "Basic Octree Construction" << endl ;
			cout <<"===> Press +/- to explore levels" << endl;				
			printKeyboardInfos();
			printMenuInfos();
			break;
		case 1 :
			cout << "\n\n\n\n\n\n---------------------" << endl;
			cout << "|   Current State   |" << endl;
			cout << "---------------------" << endl;
			cout << "1 - Morton traversal" << endl;			
			cout << "Basic Morton octree Traversal" << endl;			
			cout <<"===> Press +/- to explore levels" << endl;	
			printKeyboardInfos();
			printMenuInfos();
			break;		
		case 2 :
			cout << "\n\n\n\n\n\n---------------------" << endl;
			cout << "|   Current State   |" << endl;
			cout << "---------------------" << endl;
			cout << "2 - UAV Test Case presentation" << endl;			
			printKeyboardInfos();
			printMenuInfos();
			break;			
		/** Demos : LB FMM-lib **/		
		case 3 :
			cout << "\n\n\n\n\n\n---------------------" << endl;
			cout << "|   Current State   |" << endl;
			cout << "---------------------" << endl;
			cout << endl <<"3 - Exact Histograms " << endl;
			cout << "DEMO : Exact Histogram on UAV Test case with octree construction" << endl; 
			cout <<"  - Octree : constructed" << endl;			
			cout <<"===> Press +/- to explore levels" << endl;	
			printKeyboardInfos();
			printMenuInfos();
			break;
		case 4 :
			cout << "\n\n\n\n\n\n---------------------" << endl;
			cout << "|   Current State   |" << endl;
			cout << "---------------------" << endl;
			cout << endl << "4 - Approximate Histograms " << endl;
			cout << "DEMO : Approximate Histogram on UAV Test case with octree construction" << endl; 
			cout <<"  - Octree : constructed" << endl;			
			cout <<"===> Press +/- to explore levels" << endl;	
			printKeyboardInfos();
			printMenuInfos();
			break;
		case 5 : 
			cout << "\n\n\n\n\n\n---------------------" << endl;
			cout << "|   Current State   |" << endl;
			cout << "---------------------" << endl;
			cout << endl << "5 - Morton" << endl;
			cout << "DEMO : Morton on UAV Test case with octree construction" << endl; 			
			cout <<"  - Octree : constructed" << endl;
			cout <<"===> Press +/- to display Morton Line and explore levels" << endl;	
			printKeyboardInfos();
			printMenuInfos();
			break;
		/** Demos : With Spectre **/				
		case 6 : 
			cout << "\n\n\n\n\n\n---------------------" << endl;
			cout << "|   Current State   |" << endl;
			cout << "---------------------" << endl;
			cout << endl << "6 - Original LB " << endl;			
			cout << "ORIGINAL Load Balancing on UAV Test case" << endl; 			 
			cout <<"  - Octree : Imported" << endl;			
			cout <<"===> Press +/- to explore levels" << endl;	
			printKeyboardInfos();
			printMenuInfos();
			break;
		case 7 :
			cout << "\n\n\n\n\n\n---------------------" << endl;
			cout << "|   Current State   |" << endl;
			cout << "---------------------" << endl;
			cout << endl << "7 - Morton " << endl;
			cout << "RESULT : Morton on UAV Test case with imported octree" << endl; 						 
			cout <<"  - Octree : Imported" << endl;			
			cout <<"===> Press +/- to explore levels" << endl;	
			printKeyboardInfos();
			printMenuInfos();
			break;
		case 8 :			
			cout << "\n\n\n\n\n\n---------------------" << endl;
			cout << "|   Current State   |" << endl;
			cout << "---------------------" << endl;
			cout << endl << "8 - Histogram " << endl;
			cout << "RESULT : Histogram on UAV Test case with imported octree" << endl; 						 			
			cout <<"  - Octree : Imported" << endl;			
			cout <<"===> Press +/- to explore levels" << endl;	
			printKeyboardInfos();
			printMenuInfos();			
			break;
		default :
			cout << "default case, no display : " << state << endl;
			printMenuInfos();

	}
}

/**
 * Load Coordinates from file
 **/

void loadCoordinatesASCIIWithoutQuantity(int & nbParticles, vec3D *& coords, const string & file)
{
	fstream in;
	in.exceptions(ifstream::failbit|ifstream::badbit);
	try
	{
		in.open(file, ifstream::in);
		
		// get number of nodes
		nbParticles = std::count(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>(), '\n');
		coords = new vec3D[nbParticles];
		
		// load into _coordinates
		in.seekg(0, ios::beg);
		int index = 0;
		while (index < nbParticles)
		{
			in >> coords[index].x;
			in >> coords[index].y;
			in >> coords[index].z;			
			index++;
		}
		in.close();		
	}
	catch(ifstream::failure &e)
	{
		cerr << "[erreur] " << e.what() << endl;
		exit(0);
	}
}
