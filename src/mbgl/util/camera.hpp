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

    vec3 getPosition() const;
    mat4 getCameraToWorld(double zoom, bool flippedY) const;
    mat4 getWorldToCamera(double zoom, bool flippedY) const;
    mat4 getCameraToClipPerspective(double nearZ, double farZ) const;

    double getFovY() const;
    void setFovY(double fov);

    vec3 forward() const;
    vec3 right() const;
    vec3 up() const;

    void lookAtPoint(const LatLng& location);

    const Quaternion& getOrientation() const { return orientation; }
    void setOrientation(float pitch, float bearing);

    // Set position in mercator coordinates
    void setPosition(const vec3& mercatorLocation);

    // Set location in LatLng coordinates and altitude in meters
    void setPosition(const LatLng& location, double altitudeMeters);

    // Set location in LatLng. Elevation is deduced from the zoom value
    void setPositionZoom(const LatLng& location, double zoom);

private:
    Size size;
    double fovy;
    Quaternion orientation;
    mat4 cameraTransform;       // Position (mercator) and orientation of the camera
};

}
}