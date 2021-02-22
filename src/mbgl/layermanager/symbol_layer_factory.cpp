#include <iostream>

#include <mbgl/layermanager/symbol_layer_factory.hpp>

#include <mbgl/layout/symbol_layout.hpp>
#include <mbgl/renderer/layers/render_symbol_layer.hpp>
#include <mbgl/style/layers/symbol_layer.hpp>
#include <mbgl/style/layers/symbol_layer_impl.hpp>

namespace mbgl {

const style::LayerTypeInfo* SymbolLayerFactory::getTypeInfo() const noexcept {
    return style::SymbolLayer::Impl::staticTypeInfo();
}

std::unique_ptr<style::Layer> SymbolLayerFactory::createLayer(const std::string& id, const style::conversion::Convertible& value) noexcept {
    optional<std::string> source = getSource(value);
    if (!source) {
        return nullptr;
    }
    return std::unique_ptr<style::Layer>(new style::SymbolLayer(id, *source));
}

std::unique_ptr<Layout> SymbolLayerFactory::createLayout(const LayoutParameters& parameters,
                                                         std::unique_ptr<GeometryTileLayer> tileLayer,
                                                         const std::vector<Immutable<style::LayerProperties>>& group) noexcept {
    try {
        return std::make_unique<SymbolLayout>(parameters.bucketParameters, group, std::move(tileLayer), parameters);
    } catch(std::runtime_error& e) {
        std::cerr << "[Mapbox GL Native] Encountered error decoding feature: " << e.what() << "\n";
        return nullptr;
    }
}

std::unique_ptr<RenderLayer> SymbolLayerFactory::createRenderLayer(Immutable<style::Layer::Impl> impl) noexcept {
    assert(impl->getTypeInfo() == getTypeInfo());
    return std::make_unique<RenderSymbolLayer>(staticImmutableCast<style::SymbolLayer::Impl>(impl));
}

} // namespace mbgl
