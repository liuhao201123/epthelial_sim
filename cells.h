/*************************************************************************
	> File Name: cells.h
	> Author: 
	> Mail: 
	> Created Time: 2016年12月14日 星期三 19时00分13秒
 ************************************************************************/

#ifndef _CELLS_H
#define _CELLS_H
#endif

#include<iostream>
#include<iomanip>
using namespace std;

class cells
{
    public:
        cells(double a, double b)
        {
            x = a;
            y = b;
        }
        void set(double a, double b)
        {
            x = a;
            y = b;
        }
        void get()
        {
            cout << left<< setw(12)<<x;
            cout << right<<setw(12)<<y<<endl;
        }
		double getLongx();
		double getLongy();
		double getLongLength();
		void setLongx(double a);
		void setLongy(double a);
		void setLongLength(double a);
		
        double getx();
        double gety();
		void setArea(double a);
		void setPerimeter(double a);
		void setArea0(double a);
		void setPerimeter0(double a);
		double getArea();
		double getArea0();
		double getPerimeter();
		double getPerimeter0();
    private:
        double x, y, area, perimeter;
		double Longx, Longy, LongLength;
		double area0, perimeter0;
};

double cells::getLongx()
{
	return Longx;
}
double cells::getLongy()
{
	return Longy;
}
double cells::getLongLength()
{
	return LongLength;
}

void cells::setLongx(double a)
{
	Longx = a;
}
void cells::setLongy(double a)
{
	Longy = a;
}
void cells::setLongLength(double a)
{
	LongLength = a;
}

void cells::setPerimeter(double a)
{
	perimeter = a;
}
double cells::getPerimeter()
{
	return perimeter;
}
void cells::setPerimeter0(double a)
{
	perimeter0 = a;
}
double cells::getPerimeter0()
{
	return perimeter0;
}
void cells::setArea0(double a)
{
	area0 = a;
}
double cells::getArea0()
{
	return area0;
}
void cells::setArea(double a)
{
	area = a;
}
double cells::getArea()
{
	return area;
}
double cells::getx()
{
    return x;
}
double cells::gety()
{
    return y;
}


