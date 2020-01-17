#define  _CRT_SECURE_NO_WARNINGS 

#include <array>
#include <vector>
#include <math.h>
#include <cstdio>

#include "gif.h"

template <unsigned int Dim, typename DataType>
class Lagrangian
{
public:
	Lagrangian(std::array<DataType, Dim> _position,
		std::array<DataType, Dim> _velocity) 
		:position(_position), velocity(_velocity) {};
	~Lagrangian() = default;


	std::array<DataType, Dim> position;
	std::array<DataType, Dim> velocity;

	virtual DataType calc_lagrangian() { return DataType(0); };

	void update_p()
	{
		for (size_t i = 0; i < Dim; i++)
		{
			position[i] += velocity[i] * time_q;
		}
	}

	void update_v()
	{
		for (size_t i = 0; i < Dim; i++)
		{
			find_v(partial_derivative_v(i) + partial_derivative_p(i)*time_q, i);
		}
	}

private:
	DataType partial_derivative_p(int i)
	{
		auto temp = position[i];
		auto a = calc_lagrangian();
		position[i] += delta;
		auto b = calc_lagrangian();
		position[i] = temp;
		return (b - a) / delta;
	}
	DataType partial_derivative_v(int i)
	{
		auto temp = velocity[i];
		auto a = calc_lagrangian();
		velocity[i] += delta;
		auto b = calc_lagrangian();
		velocity[i] = temp;
		return (b - a) / delta;
	}

	DataType partial_sq_derivative_v(int i)
	{
		auto temp = velocity[i];
		auto a = calc_lagrangian();
		velocity[i] += delta;
		auto b = calc_lagrangian();
		velocity[i] += delta;
		auto c = calc_lagrangian();
		velocity[i] = temp;
		return (a - 2 * b + c) / delta / delta;
	}

	void find_v(DataType target_result, int i)
	{
		while (true)
		{
			auto pd = partial_sq_derivative_v(i);
			auto f = partial_derivative_v(i) - target_result;
			velocity[i] -= f / pd;
			if (abs(f) < precise)
				break;
		}
	}
	const DataType delta = 10E-5;
	const DataType precise = 10E-4;
	const DataType time_q = 10E-4;
};

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



int main()
{
	int width = 400;
	int height = 400;
	std::vector<uint8_t> one_frame(width * height * 4, 255);



	auto fileName = "E:\\bwgif.gif";
	int delay = 10;
	GifWriter g;
	GifBegin(&g, fileName, width, height, delay);

	MyClass m({ 0.5, -0.8 }, { 0,0 });

	one_frame[200 * width * 4 + 200 * 4] = 0;
	one_frame[200 * width * 4 + 200 * 4 + 1] = 0;
	one_frame[200 * width * 4 + 200 * 4 + 2] = 0;
	one_frame[200 * width * 4 + 200 * 4 + 3] = 0;

	for (size_t i = 0; i < 100000; i++)
	{
		if (i % 100 == 0)
		{	
			int x1, y1, x2, y2;
			x1 = 200 + 100 * cos(m.position[0]);
			y1 = 200 + 100 * sin(m.position[0]);
			x2 = x1 + 100 * cos(m.position[1]);
			y2 = y1 + 100 * sin(m.position[1]);
			std::swap(x1, y1);
			std::swap(x2, y2);
			one_frame[y1 * width * 4 + x1 * 4] = 0;
			one_frame[y1 * width * 4 + x1 * 4 + 1] = 0;
			one_frame[y1 * width * 4 + x1 * 4 + 2] = 0;
			one_frame[y1 * width * 4 + x1 * 4 + 3] = 0;
			one_frame[y2 * width * 4 + x2 * 4] = 0;
			one_frame[y2 * width * 4 + x2 * 4 + 1] = 0;
			one_frame[y2 * width * 4 + x2 * 4 + 2] = 0;
			one_frame[y2 * width * 4 + x2 * 4 + 3] = 0;
			GifWriteFrame(&g, one_frame.data(), width, height, delay);
			one_frame[y1 * width * 4 + x1 * 4] = 255;
			one_frame[y1 * width * 4 + x1 * 4 + 1] = 255;
			one_frame[y1 * width * 4 + x1 * 4 + 2] = 255;
			one_frame[y1 * width * 4 + x1 * 4 + 3] = 255;
			one_frame[y2 * width * 4 + x2 * 4] = 255;
			one_frame[y2 * width * 4 + x2 * 4 + 1] = 255;
			one_frame[y2 * width * 4 + x2 * 4 + 2] = 255;
			one_frame[y2 * width * 4 + x2 * 4 + 3] = 255;
		}
		m.update_p();
		m.update_v();
	}
	GifEnd(&g);
	return 0;
}