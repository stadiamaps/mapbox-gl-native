#pragma once

#include "mat3.hpp"
#include "mat4.hpp"

namespace mbgl {

struct Quaternion{
    union {
        vec4 m;
        struct {
            double x, y, z, s;
        };
    };

    Quaternion()
        : Quaternion(0.0, 0.0, 0.0, 0.0)
    { }

    Quaternion(double x_, double y_, double z_, double s_)
        : x(x_), y(y_), z(z_), s(s_)
    { }

    Quaternion conjugate() const;
    Quaternion multiply(const Quaternion& o) const;
    vec3 transform(const vec3& v) const;
    mat4 toRotationMatrix() const;

    static Quaternion fromAxisAngle(const vec3& axis, double angleRad);
    static Quaternion fromEulerAngles(double x, double y, double z);

    static Quaternion identity;
};

}