#include "camera.hpp"
#include <mbgl/util/constants.hpp>
#include <mbgl/util/geo.hpp>
#include <mbgl/util/projection.hpp>
#include <assert.h>
#include <cmath>

namespace mbgl {
namespace util {

static double* getColumn(mat4& matrix, int col) {
    assert(col >= 0 && col < 4);
    return &matrix[col * 4];
}

static mat4 updateCameraTransform(const Quaternion& orientation, const double* translation) {
    // Construct rotation matrix from orientation
    mat4 m = orientation.toRotationMatrix();

    // Apply translation to the matrix
    double* col = getColumn(m, 3);

    col[0] = translation[0];
    col[1] = translation[1];
    col[2] = translation[2];

    return m;
}

static void updateTransform(const mat4& camera, bool flipY, double pixelsPerMeters, mat4& outTransform, mat4& outInvTransform) {
    mat4 flipMatrix;
    mat4 zScaleMatrix;
    
    matrix::identity(flipMatrix);
    matrix::scale(flipMatrix, flipMatrix, 1.0, flipY ? 1.0 : -1.0, 1.0);

    // TODO: Scale from meters to pixels
    matrix::identity(zScaleMatrix);
    matrix::scale(zScaleMatrix, zScaleMatrix, 1.0, 1.0, pixelsPerMeters);

    // camera matrix transform camera from local coordinate space to world.
    // Elevation scale matrix transforms renderables to world space as well.
    // This means that for full world->camera transformation we need:
    // outInvTransform = (flip * cam^-1 * zScale)^-1 => (zScale^-1 * cam * flip^-1)
    
    mat4 invCameraTransform;
    matrix::invert(invCameraTransform, camera);

    matrix::identity(outInvTransform);
    matrix::multiply(outInvTransform, outInvTransform, flipMatrix);
    matrix::multiply(outInvTransform, outInvTransform, invCameraTransform);
    matrix::multiply(outInvTransform, outInvTransform, zScaleMatrix);
    matrix::invert(outTransform, outInvTransform);
}

Camera::Camera() : size(0, 0), fovy(1.0), orientation(Quaternion::identity), flippedY(false), pixelsPerMeter(1.0) {
    matrix::identity(cameraTransform);
    matrix::identity(transform);
    matrix::identity(invTransform);
    matrix::identity(projection);
    matrix::identity(invProjection);
}

void Camera::perspective(double fovY, double aspectRatio, double nearZ, double farZ) {
    fovy = fovY;
    matrix::perspective(projection, fovY, aspectRatio, nearZ, farZ);
    matrix::invert(invProjection, projection);
}

void Camera::setFlippedY(bool flipped) {
    flippedY = flipped;
    updateTransform(cameraTransform, flippedY, pixelsPerMeter, transform, invTransform);
}

void Camera::setOrientation(const Quaternion& orientation_) {
    orientation = orientation_;
    cameraTransform = updateCameraTransform(orientation, getColumn(cameraTransform, 3));
    updateTransform(cameraTransform, flippedY, pixelsPerMeter, transform, invTransform);
}

void Camera::setPosition(const vec3& position, double pixelsPerMeter_) {
    pixelsPerMeter = pixelsPerMeter_;
    cameraTransform = updateCameraTransform(orientation, position.data());
    updateTransform(cameraTransform, flippedY, pixelsPerMeter, transform, invTransform);
}

void Camera::setPosition(const LatLng& position, double zoom) {
    const double distanceToMap = 0.5 * size.height / std::tan(fovy / 2.0);
    // Convert LatLng to pixel coordinates
    const double worldSize = std::pow(2.0, zoom) * util::tileSize;
    const double bc = worldSize / util::DEGREES_MAX;
    const double cc = worldSize / util::M2PI;

    const double m = 1 - 1e-15;
    const double f = util::clamp(std::sin(util::DEG2RAD * position.latitude()), -m, m);

    vec3 pixelPos = {
        -position.longitude() * bc,
        0.5 * cc * std::log((1 + f) / (1 - f)),
        distanceToMap
    };

    const double dx = pixelPos[0] - 0.5 * worldSize;
    const double dy = pixelPos[1] - 0.5 * worldSize;
    pixelPos[0] = -dx;
    pixelPos[1] = -dy;
    //printf("%f %f\n", dx, dy);
    pixelsPerMeter = 1.0 / Projection::getMetersPerPixelAtLatitude(position.latitude(), zoom);

    cameraTransform = updateCameraTransform(orientation, pixelPos.data());
    updateTransform(cameraTransform, flippedY, pixelsPerMeter, transform, invTransform);
}

void Camera::setSize(const Size& size_) {
    size = size_;
}

}
}