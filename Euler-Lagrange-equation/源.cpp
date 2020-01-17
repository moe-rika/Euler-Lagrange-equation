#include <array>
#include <math.h>
#include <cstdio>

template <unsigned int Dim,typename DataType>
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

	void find_v(DataType target_result,int i)
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
	const DataType delta = 10E-7;
	const DataType precise = 10E-4;
	const DataType time_q = 10E-3;
};

class MyClass : public Lagrangian<1,double>
{
public:
	MyClass(std::array<double, 1> _position, std::array<double, 1> _velocity) :Lagrangian(_position, _velocity) {};
	~MyClass() = default;
	const double L = 1;
	const double m = 1;
	const double g = 9.8;
	double calc_lagrangian() override {
		double T = m * L * velocity[0] * L*velocity[0] / 2;
		double V = -L * m * g * sin(position[0]);
		return T - V; 
	};
};

int main()
{
	MyClass m({ 0 }, { 0 });
	for (size_t i = 0; i < 10000; i++)
	{	
		printf("%f , %f\n",m.position[0],m.velocity[0]);
		m.update_p();
		m.update_v();
	}
	return 0;
}