#define  _CRT_SECURE_NO_WARNINGS 

#include <array>
#include <vector>
#include <math.h>
#include <cstdio>

#include "lag.hpp"
#include "gif.h"
#include "draw.hpp"



class MyClass : public Lagrangian<2, double>
{
public:
	MyClass(std::array<double, 2> _position, std::array<double, 2> _velocity) : Lagrangian(_position, _velocity) {};
	~MyClass() = default;
	const double L = 1;
	const double m = 1;
	const double g = 9.8;
	double calc_lagrangian() override {
		return 1 / 6. * m * L * L * (\
			velocity[1] * velocity[1] + 4 * velocity[0] * velocity[0]\
			+ 3 * velocity[0] * velocity[1] * cos(position[0] - position[1]) \
			) \
			+ 1 / 2.*m*g*L*(3 * cos(position[0]) + cos(position[1]));
	};
};


class MyClass2 : public Lagrangian<3, double>
{
public:
	MyClass2(std::array<double, 3> _position, std::array<double, 3> _velocity) : Lagrangian(_position, _velocity) {};
	~MyClass2() = default;
	const double L1 = 1;
	const double L2 = 1;
	const double m1 = 1;
	const double m2 = 1;
	const double g = 9.8;

	double& a = position[0];
	double& b = position[1];
	double& c = position[2];
	double& va = velocity[0];
	double& vb = velocity[1];
	double& vc = velocity[2];

#define Power pow
#define Cos cos
#define Sin sin
	double T1()
	{
		return m1 * ((L1*Power(va, 2)) / 2. + (Power(L1, 3)*Power(vb, 2)) / 6. +
			(Power(L1, 2)*va*vb*Cos(b)) / 2.);
	}

	double T2()
	{
		return m2 * ((L2*Power(va, 2)) / 2. + L1 * L2*va*vb*Cos(b) +
			(Power(L1, 2)*L2*Power(vb, 2)*Power(Cos(b), 2)) / 2. +
			(Power(L2, 2)*va*vc*Cos(c)) / 2. +
			(L1*Power(L2, 2)*vb*vc*Cos(b)*Cos(c)) / 2. +
			(Power(L2, 3)*Power(vc, 2)*Power(Cos(c), 2)) / 6. +
			(Power(L1, 2)*L2*Power(vb, 2)*Power(Sin(b), 2)) / 2. +
			(L1*Power(L2, 2)*vb*vc*Sin(b)*Sin(c)) / 2. +
			(Power(L2, 3)*Power(vc, 2)*Power(Sin(c), 2)) / 6.);
	}

	double V()
	{
		return -1 / 2.* g *(m1*L1*cos(b) + m2 * (L1*cos(b) + L2 * cos(c)));
	}

	double calc_lagrangian() override {
		return T1() + T2() - V();
	};
};

int main()
{

	std::vector<uint8_t> one_frame(width * height * 4, 255);

	auto fileName = "E:\\Double-Pendulum53.gif";
	int delay = 10;
	GifWriter g;
	GifBegin(&g, fileName, width, height, delay);

	MyClass m({ 1.5 , 1.7 }, { 0,0 });

	int x1, y1, x2, y2, x3 = width / 2, y3 = height/2;
	//int x10, x20;
	//x10 = x3 + m.L1 * 100 * sin(m.position[1]);
	//x20 = x10 + m.L2 * 100 * sin(m.position[2]);

	for (size_t i = 0; i < 100000; i++)
	{
		if (i % 100 == 0)
		{
			x1 = x3 + 100 * sin(m.position[0]);
			y1 = y3 + 100 * cos(m.position[0]);
			x2 = x1 +  100 * sin(m.position[1]);
			y2 = y1 +  100 * cos(m.position[1]);

			//x3 = width /2 + 100 * m.position[0]; 
//error: m.position[0] is always zero
//so compute it by other method
//-------------------------------------
			//x3 = width / 2  - 1 / 2.*(m.m1*(x1 -  x10) + m.m2 * ((x2-x20) - (x1-x10)));
			//y3 = height / 2;

			draw_line(one_frame, width, x3, y3, x1, y1, 0, 0, 0, 255);
			draw_line(one_frame, width, x1, y1, x2, y2, 0, 0, 0, 255);
			GifWriteFrame(&g, one_frame.data(), width, height, delay);
			draw_line(one_frame, width, x3, y3, x1, y1, 255, 255, 255, 255);
			draw_line(one_frame, width, x1, y1, x2, y2, 255, 255, 255, 255);
			//printf("%f,%f,%f,%f,%f,%f\n", m.position[0], m.position[1], m.position[2], m.velocity[0], m.velocity[1], m.velocity[2]);
		}
		m.update();
	}
	GifEnd(&g);
	return 0;
}