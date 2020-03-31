#pragma once

#include "quaternion.hpp"

namespace mbgl {

class LatLng;

namespace util {

class Camera {
public:
    Camera();

    // Sets perspective projection for the camera
    void perspective(double fovY, double aspectRatio, double nearZ, double farZ);

    const Quaternion& getOrientation() const { return orientation; }
    const mat4& getCameraToWorld() const { return transform; }
    const mat4& getWorldToCamera() const { return invTransform; }
    const mat4& getCameraToClip() const { return projection; }

    void setFlippedY(bool flipped);

    void setOrientation(const Quaternion& orientation_);

    // Set position in pixel coordinates
    void setPosition(const vec3& position, double pixelsPerMeter_);

    // Set position in LatLng and elevation in meters
    void setPosition(const LatLng& position, double pixelElevation, double zoom);

private:
    Quaternion orientation;
    mat4 cameraTransform;       // world transformation of the camera. Position and orientation
    mat4 transform;             // World to camera space. Combination of camera and map transforms
    mat4 invTransform;
    mat4 projection;
    mat4 invProjection;
    bool flippedY;
    double pixelsPerMeter;
};

}
}