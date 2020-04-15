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

    vec3 getPosition() const;
    const Quaternion& getOrientation() const { return orientation; }
    const mat4 getCameraToWorld() const;
    const mat4 getWorldToCamera() const;
    const mat4& getCameraToClip() const { return projection; }

    double getZoom() const;
    void setZoom(double zoom);

    // Forward direction of the camera. Default value with no rotation is [0, 0, -1]
    vec3 forward() const;
    vec3 right() const;
    vec3 up() const;

    void lookAtPoint(const LatLng& location);

    void setFlippedY(bool flipped);

    void setOrientation(float pitch, float bearing);

    // Set position in mercator coordinates
    void setPosition(const vec3& mercatorLocation);

    // Set location in LatLng coordinates and elevation in meters
    void setPosition(const LatLng& location, double elevationMeters);

    // Set location in LatLng. Elevation is deduced from the zoom value
    void setPositionZoom(const LatLng& location, double zoom_);

private:
    Size size;
    double fovy;
    Quaternion orientation;
    mat4 cameraTransform;       // Position (mercator) and orientation of the camera
    mat4 projection;
    mat4 invProjection;
    bool flippedY;
    double zoom;
};

}
}