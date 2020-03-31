#include "quaternion.hpp"
#include <cmath>

namespace mbgl {

Quaternion Quaternion::identity = Quaternion(0.0, 0.0, 0.0, 1.0);

Quaternion Quaternion::conjugate() const {
    return Quaternion(-x, -y, -z, s);
}

Quaternion Quaternion::fromAxisAngle(const vec3& axis, double angleRad) {
    const double coss = std::cos(0.5 * angleRad);
    const double sins = std::sin(0.5 * angleRad);

    Quaternion q;
    q.x = sins * axis[0];
    q.y = sins * axis[1];
    q.z = sins * axis[2];
    q.s = coss;
    return q;
}

Quaternion Quaternion::fromEulerAngles(double x, double y, double z) {
    double cy = std::cos(z * 0.5);
    double sy = std::sin(z * 0.5);
    double cp = std::cos(y * 0.5);
    double sp = std::sin(y * 0.5);
    double cr = std::cos(x * 0.5);
    double sr = std::sin(x * 0.5);

    Quaternion q;
    q.s = cy * cp * cr + sy * sp * sr;
    q.x = cy * cp * sr - sy * sp * cr;
    q.y = sy * cp * sr + cy * sp * cr;
    q.z = sy * cp * cr - cy * sp * sr;

    return q;
}

Quaternion Quaternion::multiply(const Quaternion& o) const {
    Quaternion q;
    q.s = s * o.s - x * o.x - y * o.y - z * o.z;
    q.x = s * o.x + x * o.s + y * o.z - z * o.y;
    q.y = s * o.y + y * o.s + z * o.x - x * o.z;
    q.z = s * o.z + z * o.s + x * o.y - y * o.x;
    return q;
}

vec3 Quaternion::transform(const vec3& v) const {
    const Quaternion src = {v[0], v[1], v[2], 0.0};
    const Quaternion res = multiply(src).multiply(conjugate());
    return {res.x, res.y, res.z};
}

mat4 Quaternion::toRotationMatrix() const 
{
    mat4 mat;
    matrix::identity(mat);

    const double tx  = 2.0 * x;
    const double ty  = 2.0 * y;
    const double tz  = 2.0 * z;
    const double twx = tx * s;
    const double twy = ty * s;
    const double twz = tz * s;
    const double txx = tx * x;
    const double txy = ty * x;
    const double txz = tz * x;
    const double tyy = ty * y;
    const double tyz = tz * y;
    const double tzz = tz * z;

    mat[0] = 1.0 - (tyy + tzz);
    mat[1] = txy - twz;
    mat[2] = txz + twy;
    mat[4] = txy + twz;
    mat[5] = 1.0 - (txx + tzz);
    mat[6] = tyz - twx;
    mat[8] = txz - twy;
    mat[9] = tyz + twx;
    mat[10] = 1.0 - (txx + tyy);

    return mat;
}

}