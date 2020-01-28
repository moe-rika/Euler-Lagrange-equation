#define  _CRT_SECURE_NO_WARNINGS 

#include <array>
#include <vector>
#include <math.h>
#include <cstdio>

#include "lag.hpp"
#include "gif.h"
#include "draw.hpp"

#define Power pow
#define Cos cos
#define Sin sin
#define position0 position[0]
#define position1 position[1]
#define velocity0 velocity[0]
#define velocity1 velocity[1]


class MyClass : public Lagrangian<2, double>
{
public:
	MyClass(std::array<double, 2> _position, std::array<double, 2> _velocity) : Lagrangian(_position, _velocity) {};
	~MyClass() = default;
	const double L = 1;
	const double m = 1;
	const double g = 9.8;
	double T() const
	{
		return  1 / 6. * m * L * L * (\
			velocity[1] * velocity[1] + 4 * velocity[0] * velocity[0]\
			+ 3 * velocity[0] * velocity[1] * cos(position[0] - position[1]) \
			);
	}
	double V() const
	{
		return - 1 / 2.*m*g*L*(3 * cos(position[0]) + cos(position[1]));
	}

	double calc_lagrangian() const override {
		return T() - V();
	};
private:
	double partial_derivative_p(const int& i) override
	{
		switch (i)
		{
		case 0:
			return -1.5*g*L*m*Sin(position0) - 0.5*Power(L, 2)*m*velocity0*velocity1*Sin(position0 - position1);
		case 1:
			return 0.5*Power(L, 2)*m*velocity0*velocity1*Sin(position0 - position1) - 0.5*g*L*m*Sin(position1);
		default:
			return 0;
		}
	}
	double partial_derivative_v(const int& i) override
	{
		switch (i)
		{
		case 0:
			return 0.16666666666666666*Power(L, 2)*m*(8 * velocity0 + 3 * velocity1*Cos(position0 - position1));
		case 1:
			return 0.16666666666666666*Power(L, 2)*m*(2 * velocity1 + 3 * velocity0*Cos(position0 - position1));
		default:
			return 0;
		}
	}
	double partial_sq_derivative_v(const int& i) override
	{
		switch (i)
		{
		case 0:
			return 0. + 1.3333333333333333*Power(L, 2)*m;
		case 1:
			return 0. + 0.3333333333333333*Power(L, 2)*m;
		default:
			return 0;
		}
	}
};

//Mathematica code showed below

/*
l1x[a_, b_, c_, r_] := a + r Sin[b];
l1y[a_, b_, c_, r_] := r  Cos[b];
l2x[a_, b_, c_, r_] := a + L1 Sin[b] + r Sin[c];
l2y[a_, b_, c_, r_] := L1 Cos[b] + r Cos[c];*/
/*
Integrate[((D[l1x[a[t], b[t], c[t], r], t])^2 + (D[
		 l1y[a[t], b[t], c[t], r], t])^2)/2, {r, 0, L1}] /. {a[t] ->
	a, b[t] -> b, c[t] -> c, a'[t] -> va, b'[t] -> vb,
   c'[t] -> vc} // CForm
*/

/*
Integrate[((D[l2x[a[t], b[t], c[t], r], t])^2 + (D[
		 l2y[a[t], b[t], c[t], r], t])^2)/2, {r, 0, L2}] /. {a[t] ->
	a, b[t] -> b, c[t] -> c, a'[t] -> va, b'[t] -> vb,
   c'[t] -> vc} // CForm
*/

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


	double T1() const
	{
		return m1 * ((L1*Power(va, 2)) / 2. + (Power(L1, 3)*Power(vb, 2)) / 6. +
			(Power(L1, 2)*va*vb*Cos(b)) / 2.);
	}

	double T2() const
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

	double V() const
	{
		return -1 / 2.* g *(m1*L1*cos(b) + m2 * (L1*cos(b) + L2 * cos(c)));
	}

	double calc_lagrangian() const override {
		return T1() + T2() - V();
	};
};

int main()
{

	std::vector<uint8_t> one_frame(width * height * 4, 255);

	auto fileName = "E:\\Double-Pendulum59.gif";
	int delay = 10;
	GifWriter g;
	GifBegin(&g, fileName, width, height, delay);

	MyClass m({ 1.9 , 1.5 }, { 0,0 });

	int x1, y1, x2, y2, x3 = width / 2, y3 = height/2;

	for (size_t i = 0; i < 100000000*3; i++)
	{
		if (i % 1000000 == 0)
		{
			//x3 = width / 2 + 100 * m.position[0];
			//x1 = x3 + 100 * sin(m.position[1]);
			//y1 = y3 + 100 * cos(m.position[1]);
			//x2 = x1 +  100 * sin(m.position[2]);
			//y2 = y1 +  100 * cos(m.position[2]);


			//draw_line(one_frame, width, x3, y3, x1, y1, 0, 0, 0, 255);
			//draw_line(one_frame, width, x1, y1, x2, y2, 0, 0, 0, 255);
			//GifWriteFrame(&g, one_frame.data(), width, height, delay);
			//draw_line(one_frame, width, x3, y3, x1, y1, 255, 255, 255, 255);
			//draw_line(one_frame, width, x1, y1, x2, y2, 255, 255, 255, 255);
			//printf("%f,%f,%f,%f,%f,%f,%f,%f,%f\n", m.position[0], m.position[1], m.position[2], 
			//	m.velocity[0], m.velocity[1], m.velocity[2], 
			//	m.T1() + m.T2(), m.V(), m.T1() + m.T2() + m.V()); //m.T1() + m.T2() + m.V() must be a constant value
			//printf("%d\n",x1+2*x2+x3);// must be a constant value

			x1 = x3 + 100 * sin(m.position[0]);
			y1 = y3 + 100 * cos(m.position[0]);
			x2 = x1 + 100 * sin(m.position[1]);
			y2 = y1 + 100 * cos(m.position[1]);


			draw_line(one_frame, width, x3, y3, x1, y1, 0, 0, 0, 255);
			draw_line(one_frame, width, x1, y1, x2, y2, 0, 0, 0, 255);
			GifWriteFrame(&g, one_frame.data(), width, height, delay);
			draw_line(one_frame, width, x3, y3, x1, y1, 255, 255, 255, 255);
			draw_line(one_frame, width, x1, y1, x2, y2, 255, 255, 255, 255);
			printf("%f,%f,%f,%f,%f,%f,%f\n", m.position[0], m.position[1],
				m.velocity[0], m.velocity[1],
				m.T(), m.V(), m.T() + m.V()); //m.T1() + m.T2() + m.V() must be a constant value
			//printf("%d\n", x1 + 2 * x2 + x3);// must be a constant value
		}
		m.update();
	}
	GifEnd(&g);
	return 0;
}