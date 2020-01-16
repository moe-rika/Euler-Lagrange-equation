#include <array>

template <unsigned int Dim,typename DataType>
class Lagrangian
{
public:
	Lagrangian() = default;
	~Lagrangian() = default;
	virtual DataType calc_lagrangian() { return DataType(0); };

private:
	DataType partial_derivative_p()
	{

	}
	DataType partial_derivative_v()
	{

	}
	std::array<DataType, Dim> position;
	std::array<DataType, Dim> velocity;
};
