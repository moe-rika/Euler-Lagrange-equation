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


	virtual DataType calc_lagrangian() const { return DataType(0); };


	void update()
	{
		for (size_t i = 0; i < Dim; i++)
		{
			m_target[i] = partial_derivative_p(i)*time_q + partial_derivative_v(i);// do not change position[i]
		}
		for (size_t i = 0; i < Dim; i++)
		{
			position[i] += velocity[i] * time_q;// do not change m_target[i]
		}
		while (true)
		{
			//	double s = 0;
			//	for (size_t i = 0; i < Dim; i++)
			//	{
			//		find_v(m_target[i], i);
			//	}
			//	for (int i = Dim - 1; i >= 0; i--)
			//	{
			//		find_v(m_target[i], i);
			//		s += (partial_derivative_v(i) - m_target[i])*(partial_derivative_v(i) - m_target[i]);
			//	}
			//	if (s < precise*precise)
			//		break;
			//}
			double s = 0;
			for (size_t i = 0; i < Dim; i++)
			{
				auto pd = partial_sq_derivative_v(i);
				auto f = partial_derivative_v(i) - m_target[i];
				velocity[i] -= f / pd;
				s += f * f;
			}
			if (s < precise*precise)
				break;
		}
	}

protected:

	virtual DataType partial_derivative_p(const int& i)
	{
		auto temp = position[i];
		auto a = calc_lagrangian();
		position[i] += delta;
		auto b = calc_lagrangian();
		position[i] = temp;
		return (b - a) / delta;
	}
	virtual DataType partial_derivative_v(const int& i)
	{
		auto temp = velocity[i];
		auto a = calc_lagrangian();
		velocity[i] += delta;
		auto b = calc_lagrangian();
		velocity[i] = temp;
		return (b - a) / delta;
	}

	virtual DataType partial_sq_derivative_v(const int& i)
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

private:

	//void find_v(const DataType& target_result, int i)
	//{
	//	while (true)
	//	{
	//		auto pd = partial_sq_derivative_v(i);
	//		auto f = partial_derivative_v(i) - target_result;
	//		velocity[i] -= f / pd;
	//		if (abs(f) < precise/10)
	//			break;
	//	}
	//}
	const DataType delta = 10E-5;
	const DataType precise = 10E-10;
	const DataType time_q = 10E-8;

	std::array<DataType, Dim> m_target;
};