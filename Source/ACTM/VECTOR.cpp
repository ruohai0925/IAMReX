//定义了一矢量类，三维的，包括加减和赋值运算
#include "VECTOR.H"

VECTOR VECTOR::operator+(const VECTOR& V) const
{
    VECTOR ans;
    ans.x = this->x + V.x;
    ans.y = this->y + V.y;
    ans.z = this->z + V.z;
    return ans;
}

VECTOR& VECTOR::operator=(const VECTOR &V)
{
    x = V.x;
    y = V.y;
    z = V.z;
    return *this;
}

VECTOR VECTOR::operator-(const VECTOR &V) const
{
    VECTOR ans;
    ans.x = this->x - V.x;
    ans.y = this->y - V.y;
    ans.z = this->z - V.z;
    return ans;
}

VECTOR VECTOR::operator-() const
{
    return VECTOR(-x,-y,-z);
}

// 矢量与标量的乘法
VECTOR VECTOR::operator*(double scalar) const {
    return VECTOR(x * scalar, y * scalar, z * scalar);
}
//点乘
double VECTOR::operator*(const VECTOR &V) const
{
    return x * V.x + y * V.y + z * V.z;
}

VECTOR VECTOR::operator^(const VECTOR &V) const
{
    return VECTOR(y * V.z - z * V.y, z * V.x - x * V.z, x * V.y - y * V.x);
}
bool VECTOR::operator==(const VECTOR &V) const
{
    return (x == V.x) && (y == V.y) && (z == V.z);
}

bool VECTOR::operator!=(const VECTOR &V) const
{
    return !(*this == V);
}

double abs(const VECTOR &vec)
{
    return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

// 矢量与标量的乘法
VECTOR operator*(double scalar, const VECTOR& vec) {
    return VECTOR(vec.x * scalar, vec.y * scalar, vec.z * scalar);
}

VECTOR::VECTOR(double _x, double _y, double _z)
{
    x=_x;
    y=_y;
    z=_z;
}

VECTOR VECTOR::Unit() const
{
    double mag = abs(*this);
    return (*this) * (1.0 / (mag ? mag : 1.0));
}

std::ostream& operator<<(std::ostream& os, const VECTOR& vec) {
    os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
    return os;
}
