/*************************************************************************
	> File Name: main.cpp
	> Author: 
	> Mail: 
	> Created Time: 2016年12月14日 星期三 16时44分11秒
 ************************************************************************/

#include<iostream>
#include<vector>
#include<string>
#include<cmath>
#include "vertex_model.h"
using namespace std;


int main(int argc, char *argv[])
{
	
	string str0 = argv[1];
	string str3 = argv[2];
    int total_step = atoi(argv[3]);
    double seedi = atof(argv[4]);
	int cell_numbers = atoi(argv[5]);
	
	int nn = 0;
	double lx = 1.0, ly = 1.0;
    vertex_model vertex(cell_numbers);

    vertex.setRatio(lx, ly);
    vertex.initial(2);
    vertex.setNoise( seedi);

//	string str0("E:\\work\\codes\\data\\votex_model2d\\");
	nn = vertex.Honda(nn, str0, str3, 1);
	
    for( int i = 0; i< 5; i++ )
    {
        lx = lx + 0.1;
        vertex.setRatio(lx, ly);
        nn = vertex.Honda(nn, str0, str3, 1);
    }
    vertex.setThreshold0(vertex.getThreshold());
    
    nn = 0;
    vertex.setVolume0();
	for( int i = 0; i< total_step; i++)
	{

		cout << "========================================== Division  "<<  i<<endl;
		nn = vertex.Honda(nn, str0, str3, 1);
		int elongation_sig = vertex.Division();
        vertex.getLength_Width( &lx, &ly);
        cout << "  lx is " << lx <<"  ly is "<<ly<<endl;
        vertex.setRatio(lx, ly);
	}

}



