#pragma once
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
	std::array<DataType, Dim> old_velocity;

	virtual DataType calc_lagrangian() { return DataType(0); };

	void update()
	{
		old_velocity = velocity;
		update_v();
		update_p();
	}

private:
	DataType partial_derivative_p(int i)
	{
		auto temp = position[i];
		auto a = calc_lagrangian();
		position[i] += delta;
		auto b = calc_lagrangian();
		position[i] = temp;
		//if (i == 2)
		//	printf("%f\n", (b - a) / delta);
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

	void update_p()
	{
		for (size_t i = 0; i < Dim; i++)
		{
			position[i] += old_velocity[i] * time_q;
		}
	}

	void update_v()
	{
		for (size_t i = 0; i < Dim; i++)
		{
			find_v(partial_derivative_v(i) + partial_derivative_p(i)*time_q, i);
		}
	}

};