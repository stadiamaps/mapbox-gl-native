#include <mbgl/renderer/raster_bucket.hpp>
#include <mbgl/renderer/painter.hpp>

using namespace mbgl;

RasterBucket::RasterBucket(const std::shared_ptr<Texturepool> &texturepool, const StyleBucketRaster& properties)
: properties(properties),
  texture(properties),
  raster(texturepool) {
}

void RasterBucket::render(Painter &painter, std::shared_ptr<StyleLayer> layer_desc, const Tile::ID &id) {
    painter.renderRaster(*this, layer_desc, id);
}

bool RasterBucket::setImage(const std::string &data) {
    return raster.load(data);
}

void RasterBucket::drawRaster(RasterShader& shader, VertexBuffer &vertices, VertexArrayObject &array) {
    raster.bind(true);
    shader.setImage(0);
    array.bind(shader, vertices, BUFFER_OFFSET(0));
    glDrawArrays(GL_TRIANGLES, 0, (GLsizei)vertices.index());
}

bool RasterBucket::hasData() const {
    return raster.isLoaded();
}
