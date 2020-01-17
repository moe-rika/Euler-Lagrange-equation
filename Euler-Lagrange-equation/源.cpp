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

void draw_pixel(std::vector<uint8_t>& v,int _width,int x,int y,uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
	v[x * _width * 4 + y * 4] = r;
	v[x * _width * 4 + y * 4 + 1] = g;
	v[x * _width * 4 + y * 4 + 2] = b;
	v[x * _width * 4 + y * 4 + 3] = a;
}

void draw_point(std::vector<uint8_t>& v, int _width, int ra, int x, int y, 
	uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
	for (int i = -ra; i < ra; i++)
	{
		for (int j = -ra; j < ra; j++)
		{
			if(i*i + j*j <= ra*ra)
				draw_pixel(v, _width, x + i, y + j, r, g, b, a);
		}
	}
}

// 交换整数 a 、b 的值
inline void swap_int(int *a, int *b) {
	*a ^= *b;
	*b ^= *a;
	*a ^= *b;
}

// Bresenham's line algorithm
void draw_line(std::vector<uint8_t>& v, int _width, int x1, int y1, int x2, int y2, uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
	// 参数 c 为颜色值
	int dx = abs(x2 - x1),
		dy = abs(y2 - y1),
		yy = 0;

	if (dx < dy) {
		yy = 1;
		swap_int(&x1, &y1);
		swap_int(&x2, &y2);
		swap_int(&dx, &dy);
	}

	int ix = (x2 - x1) > 0 ? 1 : -1,
		iy = (y2 - y1) > 0 ? 1 : -1,
		cx = x1,
		cy = y1,
		n2dy = dy * 2,
		n2dydx = (dy - dx) * 2,
		d = dy * 2 - dx;

	if (yy) { // 如果直线与 x 轴的夹角大于 45 度
		while (cx != x2) {
			if (d < 0) {
				d += n2dy;
			}
			else {
				cy += iy;
				d += n2dydx;
			}
			draw_point(v, _width, 4, cy, cx, r, g, b, a);
			cx += ix;
		}
	}
	else { // 如果直线与 x 轴的夹角小于 45 度
		while (cx != x2) {
			if (d < 0) {
				d += n2dy;
			}
			else {
				cy += iy;
				d += n2dydx;
			}
			draw_point(v, _width, 4, cx, cy, r, g, b, a);
			cx += ix;
		}
	}
}


int main()
{
	int width = 400;
	int height = 400;
	std::vector<uint8_t> one_frame(width * height * 4, 255);



	auto fileName = "E:\\Double-Pendulum.gif";
	int delay = 10;
	GifWriter g;
	GifBegin(&g, fileName, width, height, delay);

	MyClass m({ 0.5, -0.8 }, { 0,0 });

	
	for (size_t i = 0; i < 100000; i++)
	{
		if (i % 100 == 0)
		{	
			int x1, y1, x2, y2;
			x1 = 200 + 100 * cos(m.position[0]);
			y1 = 200 + 100 * sin(m.position[0]);
			x2 = x1 + 100 * cos(m.position[1]);
			y2 = y1 + 100 * sin(m.position[1]);
			draw_line(one_frame, width, 200, 200, x1, y1, 0, 0, 0, 255);
			draw_line(one_frame, width, x1, y1, x2, y2, 0, 0, 0, 255);
			GifWriteFrame(&g, one_frame.data(), width, height, delay);
			draw_line(one_frame, width, 200, 200, x1, y1, 255, 255, 255, 255);
			draw_line(one_frame, width, x1, y1, x2, y2, 255, 255, 255, 255);
		}
		m.update_p();
		m.update_v();
	}
	GifEnd(&g);
	return 0;
}