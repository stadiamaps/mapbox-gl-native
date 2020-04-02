#include "camera.hpp"
#include <mbgl/util/constants.hpp>
#include <mbgl/util/geo.hpp>
#include <mbgl/util/projection.hpp>
#include <assert.h>
#include <cmath>

namespace mbgl {
namespace util {

//export function latFromMercatorY(y: number) {
//    const y2 = 180 - y * 360;
//    return 360 / Math.PI * Math.atan(Math.exp(y2 * Math.PI / 180)) - 90;
//}

static double latFromMercatorY(double mercatorY) {
    return 2.0 * std::atan(std::exp(M_PI - mercatorY * 2.0 * M_PI)) - M_PI_2;
}

static double* getColumn(mat4& matrix, int col) {
    assert(col >= 0 && col < 4);
    return &matrix[col * 4];
}

static const double* getColumn(const mat4& matrix, int col) {
    assert(col >= 0 && col < 4);
    return &matrix[col * 4];
}

static double computePixelsPerMeter(double mercatorY, double zoom) {
    const double latitude = latFromMercatorY(mercatorY);
    return 1.0 / Projection::getMetersPerPixelAtLatitude(latitude, zoom);
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

static void updateTransform(const mat4& camera, bool flipY, mat4& outTransform, mat4& outInvTransform) {
    mat4 flipMatrix;
    mat4 zScaleMatrix;
    
    matrix::identity(flipMatrix);
    matrix::scale(flipMatrix, flipMatrix, 1.0, flipY ? 1.0 : -1.0, 1.0);

    // TODO: Remove scale from here!
    matrix::identity(zScaleMatrix);
    matrix::scale(zScaleMatrix, zScaleMatrix, 1.0, 1.0, 1.0);

    // camera matrix transform camera from local coordinate space to world.
    // Elevation scale matrix transforms renderables to world space as well.
    // This means that for full camera->world transformation we need:
    // outTransform = (flip * cam^-1 * zScale)^-1 => (zScale^-1 * cam * flip^-1)
    
    mat4 invCameraTransform;
    matrix::invert(invCameraTransform, camera);

    matrix::identity(outInvTransform);
    matrix::multiply(outInvTransform, outInvTransform, flipMatrix);
    matrix::multiply(outInvTransform, outInvTransform, invCameraTransform);
    matrix::multiply(outInvTransform, outInvTransform, zScaleMatrix);
    matrix::invert(outTransform, outInvTransform);
}

Camera::Camera() : size(0, 0), fovy(1.0), orientation(Quaternion::identity), flippedY(false) {
    matrix::identity(cameraTransform);
    matrix::identity(transform);
    matrix::identity(invTransform);
    matrix::identity(projection);
    matrix::identity(invProjection);
}

const mat4 Camera::getCameraToWorld(double zoom) const {
    // transform-matrix will transform points from camera space to mercator space
    // x, y and z coordinates have to be transformed to pixels
    mat4 cameraToWorld;

    const double scale = std::pow(2.0, zoom);
    const double worldSize = Projection::worldSize(scale);
    const double pixelsPerMeter = computePixelsPerMeter(getColumn(cameraTransform, 3)[1], zoom);

    // Mercator -> pixel coordinates
    mat4 mercatorToWorld;
    mat4 pixelsToMeters;
    matrix::identity(mercatorToWorld);
    matrix::identity(pixelsToMeters);

    matrix::scale(mercatorToWorld, mercatorToWorld, worldSize, worldSize, worldSize);
    matrix::multiply(cameraToWorld, mercatorToWorld, transform);

    // Elevation from pixel -> meters
    matrix::scale(pixelsToMeters, pixelsToMeters, 1.0, 1.0, 1.0 / pixelsPerMeter);
    matrix::multiply(cameraToWorld, pixelsToMeters, cameraToWorld);

    return cameraToWorld;
}

const mat4 Camera::getWorldToCamera(double zoom) const {
    mat4 matrix;// = getCameraToWorld();
    matrix::invert(matrix, getCameraToWorld(zoom));
    return matrix;
}

void Camera::perspective(double fovY, double aspectRatio, double nearZ, double farZ) {
    fovy = fovY;
    matrix::perspective(projection, fovY, aspectRatio, nearZ, farZ);
    matrix::invert(invProjection, projection);
}

void Camera::setFlippedY(bool flipped) {
    flippedY = flipped;
    updateTransform(cameraTransform, flippedY, transform, invTransform);
}

void Camera::setOrientation(const Quaternion& orientation_) {
    orientation = orientation_;
    cameraTransform = updateCameraTransform(orientation, getColumn(cameraTransform, 3));
    updateTransform(cameraTransform, flippedY, transform, invTransform);
}

void Camera::setPosition(const vec3& mercatorPosition) {
    cameraTransform = updateCameraTransform(orientation, mercatorPosition.data());
    updateTransform(cameraTransform, flippedY, transform, invTransform);
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
    //pixelsPerMeter = 1.0 / Projection::getMetersPerPixelAtLatitude(position.latitude(), zoom);

    cameraTransform = updateCameraTransform(orientation, pixelPos.data());
    updateTransform(cameraTransform, flippedY, transform, invTransform);
}

void Camera::setSize(const Size& size_) {
    size = size_;
}

}
}