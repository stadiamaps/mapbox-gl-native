#include "camera.hpp"
#include <mbgl/util/constants.hpp>
#include <mbgl/util/geo.hpp>
#include <mbgl/util/projection.hpp>
#include <mbgl/math/log2.hpp>
#include <assert.h>
#include <cmath>

namespace mbgl {
namespace util {

static double mercatorXfromLng(double lng) {
    return (180.0 + lng) / 360.0;
}

static double mercatorYfromLat(double lat) {
    return (180.0 - (180.0 / M_PI * std::log(std::tan(M_PI_4 + lat * M_PI / 360.0)))) / 360.0;
}

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

static double handednessFactor(const mat4& projectionMatrix) {
    // Return z-scale value of the projection matrix.
    // This is -1 for right handed projection matrices and
    // +1 for left handed matrices
    return getColumn(projectionMatrix, 2)[3];
}

static vec3 toMercator(const LatLng& location, double elevationMeters) {
    const double pixelsPerMeter = 1.0 / Projection::getMetersPerPixelAtLatitude(location.latitude(), 0.0);
    const double worldSize = Projection::worldSize(std::pow(2.0, 0.0));

    return {
         mercatorXfromLng(location.longitude()),
         mercatorYfromLat(location.latitude()),
         elevationMeters * pixelsPerMeter / worldSize
    };
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

Camera::Camera() : size(0, 0), fovy(1.0), orientation(Quaternion::identity), flippedY(false), zoom(0.0) {
    matrix::identity(cameraTransform);
    matrix::identity(projection);
    matrix::identity(invProjection);
}

vec3 Camera::getPosition() const {
    const double* p = getColumn(cameraTransform, 3);
    return { p[0], p[1], p[2] };
}

double Camera::getZoom() const {
    return zoom;
}

void Camera::setZoom(double zoom_) {
    zoom = zoom_;
}

const mat4 Camera::getCameraToWorld() const {
    mat4 cameraToWorld;
    matrix::invert(cameraToWorld, getWorldToCamera());
    return cameraToWorld;
}

const mat4 Camera::getWorldToCamera() const {

    // transformation chain from world space to camera space:
    // 1. Height value (z) of renderables is in meters. Scale z coordinate by pixelsPerMeter
    // 2. Transform from pixel coordinates to camera space with cameraMatrix^-1
    // 3. flip Y if required

    // worldToCamera: flip * cam^-1 * zScale
    // cameraToWorld: (flip * cam^-1 * zScale)^-1 => (zScale^-1 * cam * flip^-1)
    mat4 flipMatrix;
    mat4 zScaleMatrix;

    const double scale = std::pow(2.0, zoom);
    const double worldSize = Projection::worldSize(scale);
    const double latitude = latFromMercatorY(getColumn(cameraTransform, 3)[1]);
    const double pixelsPerMeter = 1.0 / Projection::getMetersPerPixelAtLatitude(latitude, zoom);
    
    // Position of the camera is stored in mercator coordinates. Scale it to pixel coordinates
    mat4 camera = cameraTransform;
    getColumn(camera, 3)[0] *= worldSize;
    getColumn(camera, 3)[1] *= worldSize;
    getColumn(camera, 3)[2] *= worldSize;

    matrix::identity(flipMatrix);
    matrix::scale(flipMatrix, flipMatrix, 1.0, flippedY ? 1.0 : -1.0, 1.0);

    matrix::identity(zScaleMatrix);
    matrix::scale(zScaleMatrix, zScaleMatrix, 1.0, 1.0, pixelsPerMeter);

    mat4 invCamera;
    mat4 result;

    matrix::invert(invCamera, camera);
    matrix::identity(result);

    matrix::multiply(result, result, flipMatrix);
    matrix::multiply(result, result, invCamera);
    matrix::multiply(result, result, zScaleMatrix);

    return result;
}

vec3 Camera::forward() const {
    const double* column = getColumn(cameraTransform, 2);
    // Front direction is determined by the handedness of the projection matrix.
    // The default direction should be towards the map, [0, 0, -1]
    const double scale = handednessFactor(projection);
    return {
        column[0] * scale,
        column[1] * scale,
        column[2] * scale
    };
}

vec3 Camera::right() const {
    const double* column = getColumn(cameraTransform, 0);
    return {
        column[0],
        column[1],
        column[2]
    };
}

vec3 Camera::up() const {
    const double* column = getColumn(cameraTransform, 1);
    // Up direction has to be flipped as [0, 0, -1] is pointing upwards in default orientation
    return {
        -column[0],
        -column[1],
        -column[2]
    };
}

void Camera::perspective(double fovY, double aspectRatio, double nearZ, double farZ) {
    fovy = fovY;
    matrix::perspective(projection, fovY, aspectRatio, nearZ, farZ);
    matrix::invert(invProjection, projection);
}

void Camera::lookAtPoint(const LatLng& location) {
    const vec3 mercator = toMercator(location, 0.0);
    const double* position = getColumn(cameraTransform, 3);

    const double dx = mercator[0] - position[0];
    const double dy = mercator[1] - position[1];
    const double dz = mercator[2] - position[2];

    const double rotZ = std::atan2(-dy, dx) - M_PI_2;
    const double rotX = std::atan2(std::sqrt(dx * dx + dy * dy), -dz);

    setOrientation(rotX, rotZ);
}

void Camera::setFlippedY(bool flipped) {
    flippedY = flipped;
}

void Camera::setOrientation(float pitch, float bearing) {
    // Both angles have to be negated to achieve CW rotation around the axis of rotation
    Quaternion rotBearing = Quaternion::fromEulerAngles(0.0, 0.0, -bearing);
    Quaternion rotPitch = Quaternion::fromEulerAngles(-pitch, 0.0, 0.0);

    orientation = rotBearing.multiply(rotPitch);
    cameraTransform = updateCameraTransform(orientation, getColumn(cameraTransform, 3));
}

void Camera::setPosition(const vec3& mercatorLocation) {
    const double pixelDistance = 0.5 * size.height / std::tan(fovy / 2.0);
    const double scale = pixelDistance / (mercatorLocation[2] * util::tileSize);

    zoom = util::log2(scale);
    cameraTransform = updateCameraTransform(orientation, mercatorLocation.data());
}

void Camera::setPosition(const LatLng& location, double elevationMeters) {
    setPosition(toMercator(location, elevationMeters));
}

void Camera::setPositionZoom(const LatLng& location, double zoom_) {
    // Pixel distance from camera to map center is always constant.
    // Use this fact to compute elevation in meters at the provided location
    const double pixelDistance = 0.5 * size.height / std::tan(fovy / 2.0);
    const double metersPerPixel = Projection::getMetersPerPixelAtLatitude(location.latitude(), zoom_);
    const double meterElevation = pixelDistance * metersPerPixel;
    
    return setPosition(location, meterElevation);
}

void Camera::setSize(const Size& size_) {
    size = size_;
}

}
}