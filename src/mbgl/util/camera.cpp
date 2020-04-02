#include "camera.hpp"
#include <mbgl/util/constants.hpp>
#include <mbgl/util/geo.hpp>
#include <mbgl/util/projection.hpp>
#include <assert.h>
#include <cmath>

namespace mbgl {
namespace util {

// export function mercatorXfromLng(lng: number) {
//     return (180 + lng) / 360;
// }

// export function mercatorYfromLat(lat: number) {
//     return (180 - (180 / Math.PI * Math.log(Math.tan(Math.PI / 4 + lat * Math.PI / 360)))) / 360;
// }

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

static vec3 toMercator(const LatLng& location, double elevationMeters) {
    const double pixelsPerMeter = 1.0 / Projection::getMetersPerPixelAtLatitude(location.latitude(), 0.0);
    const double worldSize = Projection::worldSize(std::pow(2.0, 0));

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

Camera::Camera() : size(0, 0), fovy(1.0), orientation(Quaternion::identity), flippedY(false) {
    matrix::identity(cameraTransform);
    matrix::identity(projection);
    matrix::identity(invProjection);
}

const mat4 Camera::getCameraToWorld(double zoom) const {
    mat4 cameraToWorld;
    matrix::invert(cameraToWorld, getWorldToCamera(zoom));
    return cameraToWorld;
}

const mat4 Camera::getWorldToCamera(double zoom) const {

    // transformation chain from world space to camera space:
    // 1. multiply elevation with pixelsPerMeter
    // 2. Transform from pixel coordinates to camera space with cameraMatrix^1
    // 3. flip Y if required

    // worldToCamera: flip * cam^-1 * zScale
    // cameraToWorld: (flip * cam^-1 * zScale)^-1 => (zScale^-1 * cam * flip^-1)
    mat4 flipMatrix;
    mat4 zScaleMatrix;

    const double scale = std::pow(2.0, zoom);
    const double worldSize = Projection::worldSize(scale);
    const double latitude = latFromMercatorY(getColumn(cameraTransform, 3)[1]);
    const double pixelsPerMeter = 1.0 / Projection::getMetersPerPixelAtLatitude(latitude, zoom);
    
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

void Camera::perspective(double fovY, double aspectRatio, double nearZ, double farZ) {
    fovy = fovY;
    matrix::perspective(projection, fovY, aspectRatio, nearZ, farZ);
    matrix::invert(invProjection, projection);
}

void Camera::lookAtPoint(const LatLng& location) {
    // TODO: replace euler angles!
    const vec3 mercator = toMercator(location, 0.0);

    const double dx = mercator[0] - getColumn(cameraTransform, 3)[0];
    const double dy = mercator[1] - getColumn(cameraTransform, 3)[1];
    const double dz = mercator[2] - getColumn(cameraTransform, 3)[2];

    const double rotZ = std::atan2(-dy, dx) - M_PI_2;
    const double rotX = std::atan2(std::sqrt(dx * dx + dy * dy), -dz);
    (void)rotZ;
    (void)rotX;

    Quaternion rotBearing = Quaternion::fromEulerAngles(0.0, 0.0, rotZ);
    Quaternion rotPitch = Quaternion::fromEulerAngles(rotX, 0.0, 0.0);
    Quaternion rotation = rotPitch.multiply(rotBearing);

    setOrientation(rotation);
}

void Camera::setFlippedY(bool flipped) {
    flippedY = flipped;
}

void Camera::setOrientation(const Quaternion& orientation_) {
    orientation = orientation_;
    cameraTransform = updateCameraTransform(orientation, getColumn(cameraTransform, 3));
}

void Camera::setPosition(const vec3& mercatorLocation) {
    cameraTransform = updateCameraTransform(orientation, mercatorLocation.data());
}

void Camera::setPosition(const LatLng& location, double elevationMeters) {
    vec3 position = toMercator(location, elevationMeters);
    cameraTransform = updateCameraTransform(orientation, position.data());
}

void Camera::setPositionZoom(const LatLng& location, double zoom) {
    // Pixel distance from camera to map center is always constant. Use this fact to
    // compute elevation in meters at provided location
    const double pixelElevation = 0.5 * size.height / std::tan(fovy / 2.0);
    const double metersPerPixel = Projection::getMetersPerPixelAtLatitude(location.latitude(), zoom);
    const double meterElevation = pixelElevation * metersPerPixel;
    
    return setPosition(location, meterElevation);
}

void Camera::setSize(const Size& size_) {
    size = size_;
}

}
}