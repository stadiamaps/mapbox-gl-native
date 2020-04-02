#pragma once

#include "quaternion.hpp"
#include <mbgl/util/size.hpp>

namespace mbgl {

class LatLng;

namespace util {

class Camera {
public:
    Camera();

    // Update screen size
    void setSize(const Size& size_);

    // Sets perspective projection for the camera
    void perspective(double fovY, double aspectRatio, double nearZ, double farZ);

    const Quaternion& getOrientation() const { return orientation; }
    const mat4 getCameraToWorld(double zoom) const;
    const mat4 getWorldToCamera(double zoom) const;
    const mat4& getCameraToClip() const { return projection; }

    void setFlippedY(bool flipped);

    void setOrientation(const Quaternion& orientation_);

    // Set position in mercator coordinates
    void setPosition(const vec3& mercatorPosition);

    // Set position in LatLng and elevation in meters
    void setPosition(const LatLng& position, double zoom);

private:
    Size size;
    double fovy;
    Quaternion orientation;
    mat4 cameraTransform;       // Position (mercator) and orientation of the camera
    mat4 transform;             // Full transformation from camera to world
    mat4 invTransform;          // Full transformation from world to camera
    mat4 projection;
    mat4 invProjection;
    bool flippedY;
};

}
}