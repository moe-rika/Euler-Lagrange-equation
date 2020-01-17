#define  _CRT_SECURE_NO_WARNINGS 




#include <array>
#include <math.h>
#include <cstdio>

template <unsigned int Dim, typename DataType>
class Lagrangian
{
public:
	Lagrangian(std::array<DataType, Dim> _position, std::array<DataType, Dim> _velocity) :position(_position), velocity(_velocity) {};
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
	MyClass(std::array<double, 2> _position, std::array<double, 2> _velocity) :Lagrangian(_position, _velocity) {};
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
	MyClass m({ 0.5, -0.8 }, { 0,0 });
	FILE* pFile = fopen("E:\\test.txt", "w");
	for (size_t i = 0; i < 100000; i++)
	{
		if (i % 100 == 0)
			fprintf(pFile,"%f %f\n", m.position[0], m.position[1]);
		m.update_p();
		m.update_v();
	}
	return 0;
}