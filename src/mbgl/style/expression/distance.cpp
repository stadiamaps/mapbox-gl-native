#include <mbgl/style/expression/distance.hpp>

#include <mapbox/cheap_ruler.hpp>
#include <mapbox/geojson.hpp>
#include <mapbox/geometry.hpp>

#include <mbgl/style/conversion/json.hpp>
#include <mbgl/tile/geometry_tile_data.hpp>

#include <mbgl/util/logging.hpp>
#include <mbgl/util/string.hpp>

#include <rapidjson/document.h>

#include <tuple>

namespace mbgl {
namespace {

// a, b are end points for line segment1, c and d are end points for line segment2
bool lineIntersectLine(const Point<double>& a, const Point<double>& b, const Point<double>& c, const Point<double>& d) {
    const auto perp = [](const Point<double>& v1, const Point<double>& v2) { return (v1.x * v2.y - v1.y * v2.x); };

    // check if two segments are parallel or not
    // precondition is end point a, b is inside polygon, if line a->b is
    // parallel to polygon edge c->d, then a->b won't intersect with c->d
    auto vectorP = Point<double>(b.x - a.x, b.y - a.y);
    auto vectorQ = Point<double>(d.x - c.x, d.y - c.y);
    if (perp(vectorQ, vectorP) == 0) return false;

    // check if p1 and p2 are in different sides of line segment q1->q2
    const auto twoSided =
        [](const Point<double>& p1, const Point<double>& p2, const Point<double>& q1, const Point<double>& q2) {
            double x1;
            double y1;
            double x2;
            double y2;
            double x3;
            double y3;

            // q1->p1 (x1, y1), q1->p2 (x2, y2), q1->q2 (x3, y3)
            x1 = p1.x - q1.x;
            y1 = p1.y - q1.y;
            x2 = p2.x - q1.x;
            y2 = p2.y - q1.y;
            x3 = q2.x - q1.x;
            y3 = q2.y - q1.y;
            auto ret1 = (x1 * y3 - x3 * y1);
            auto ret2 = (x2 * y3 - x3 * y2);
            return ((ret1 > 0 && ret2 < 0) || (ret1 < 0 && ret2 > 0));
        };

    // If lines are intersecting with each other, the relative location should be:
    // a and b lie in different sides of segment c->d
    // c and d lie in different sides of segment a->b
    return (twoSided(a, b, c, d) && twoSided(c, d, a, b));
}

double shortestDistanceToLine(const mapbox::geometry::point<double>& point,
                              const mapbox::geometry::line_string<double>& line,
                              mapbox::cheap_ruler::CheapRuler& ruler) {
    const auto nearestPoint = std::get<0>(ruler.pointOnLine(line, point));
    return ruler.distance(point, nearestPoint);
}

double shortestDistanceToLines(const mapbox::geometry::point<double>& point,
                               const mapbox::geometry::multi_line_string<double>& lines,
                               mapbox::cheap_ruler::CheapRuler& ruler) {
    double dist = std::numeric_limits<double>::infinity();
    for (const auto& line : lines) {
        const auto nearestPoint = std::get<0>(ruler.pointOnLine(line, point));
        auto tempDist = ruler.distance(point, nearestPoint);
        if (tempDist == 0.) return tempDist;
        dist = std::min(dist, tempDist);
    }
    return dist;
}

double shortestDistanceToPoints(const mapbox::geometry::point<double>& point,
                                const mapbox::geometry::multi_point<double>& points,
                                mapbox::cheap_ruler::CheapRuler& ruler) {
    double dist = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < points.size(); ++i) {
        auto tempDist = ruler.distance(point, points[i]);
        if (tempDist == 0.) return tempDist;
        dist = std::min(dist, tempDist);
    }
    return dist;
}

double shortestDistanceLineToLine(const mapbox::geometry::line_string<double>& line1,
                                  const mapbox::geometry::line_string<double>& line2,
                                  mapbox::cheap_ruler::CheapRuler& ruler) {
    double dist = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < line1.size() - 1; ++i) {
        const auto& p1 = line1[i];
        const auto& p2 = line1[i + 1];
        for (std::size_t j = 0; j < line2.size() - 1; ++j) {
            const auto& q1 = line2[j];
            const auto& q2 = line2[j + 1];
            if (lineIntersectLine(p1, p2, q1, q2)) return 0.;
            dist = std::min(dist, shortestDistanceToLine(p1, mapbox::geometry::line_string<double>{q1, q2}, ruler));
            dist = std::min(dist, shortestDistanceToLine(p2, mapbox::geometry::line_string<double>{q1, q2}, ruler));
            dist = std::min(dist, shortestDistanceToLine(q1, mapbox::geometry::line_string<double>{p1, p2}, ruler));
            dist = std::min(dist, shortestDistanceToLine(q1, mapbox::geometry::line_string<double>{p1, p2}, ruler));
        }
    }
    return dist;
}

double pointDistanceToGeometry(const mapbox::geometry::point<double>& point,
                               const Feature::geometry_type& geoSet,
                               mapbox::cheap_ruler::CheapRuler::Unit unit) {
    mapbox::cheap_ruler::CheapRuler ruler(point.y, unit);
    return geoSet.match([&point, &ruler](const mapbox::geometry::point<double>& p) { return ruler.distance(point, p); },
                        [&point, &ruler](const mapbox::geometry::multi_point<double>& points) {
                            return shortestDistanceToPoints(point, points, ruler);
                        },
                        [&point, &ruler](const mapbox::geometry::line_string<double>& line) {
                            return shortestDistanceToLine(point, line, ruler);
                        },
                        [&point, &ruler](const mapbox::geometry::multi_line_string<double>& lines) {
                            return shortestDistanceToLines(point, lines, ruler);
                        },
                        [](const auto&) -> double { return std::numeric_limits<double>::infinity(); });
}

double lineDistanceToGeometry(const mapbox::geometry::line_string<double>& line,
                              const Feature::geometry_type& geoSet,
                              mapbox::cheap_ruler::CheapRuler::Unit unit) {
    assert(!line.empty());
    mapbox::cheap_ruler::CheapRuler ruler(line.front().y, unit);
    return geoSet.match(
        [&line, &ruler](const mapbox::geometry::point<double>& p) { return shortestDistanceToLine(p, line, ruler); },
        [&line, &ruler](const mapbox::geometry::multi_point<double>& points) {
            double dist = std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < points.size(); ++i) {
                dist = std::min(dist, shortestDistanceToLine(points[i], line, ruler));
            }
            return dist;
        },
        [&line, &ruler](const mapbox::geometry::line_string<double>& line1) {
            return shortestDistanceLineToLine(line, line1, ruler);
        },
        [&line, &ruler](const mapbox::geometry::multi_line_string<double>& lines) {
            double dist = std::numeric_limits<double>::infinity();
            for (const auto& l : lines) {
                auto tempDist = shortestDistanceLineToLine(line, l, ruler);
                if (tempDist == 0.) return 0.;
                dist = std::min(dist, tempDist);
            }
            return dist;
        },
        [](const auto&) -> double { return -1.0; });
}

double calculateDistance(const GeometryTileFeature& feature,
                         const CanonicalTileID& canonical,
                         const Feature::geometry_type& geoSet,
                         mapbox::cheap_ruler::CheapRuler::Unit unit) {
    return convertGeometry(feature, canonical)
        .match(
            [&geoSet, &unit](const mapbox::geometry::point<double>& point) -> double {
                return pointDistanceToGeometry(point, geoSet, unit);
            },
            [&geoSet, &unit](const mapbox::geometry::multi_point<double>& points) -> double {
                double ret = std::numeric_limits<double>::infinity();
                for (const auto& p : points) {
                    auto dist = pointDistanceToGeometry(p, geoSet, unit);
                    if (dist == 0.) return dist;
                    ret = std::min(ret, dist);
                }
                return ret;
            },
            [&geoSet, &unit](const mapbox::geometry::line_string<double>& line) -> double {
                return lineDistanceToGeometry(line, geoSet, unit);
            },
            [&geoSet, &unit](const mapbox::geometry::multi_line_string<double>& lines) -> double {
                double ret = std::numeric_limits<double>::infinity();
                for (const auto& line : lines) {
                    auto dist = lineDistanceToGeometry(line, geoSet, unit);
                    if (dist == 0.) return dist;
                    ret = std::min(ret, dist);
                }
                return ret;
            },
            [](const auto&) -> double { return -1.0; });
}

struct Arguments {
    Arguments(GeoJSON& geojson_, mapbox::cheap_ruler::CheapRuler::Unit unit_)
        : geojson(std::move(geojson_)), unit(unit_) {}

    GeoJSON geojson;
    mapbox::cheap_ruler::CheapRuler::Unit unit;
};

optional<Arguments> parseValue(const style::conversion::Convertible& value, style::expression::ParsingContext& ctx) {
    if (isArray(value)) {
        // object value, quoted with ["Distance", GeoJSONObj, "unit"]
        auto length = arrayLength(value);
        if (length != 2 && length != 3) {
            ctx.error("'distance' expression requires exactly one argument, but found " +
                      util::toString(arrayLength(value) - 1) + " instead.");
            return nullopt;
        }

        mapbox::cheap_ruler::CheapRuler::Unit unit = mapbox::cheap_ruler::CheapRuler::Unit::Meters;
        if (arrayLength(value) == 3) {
            auto input = toString(arrayMember(value, 2)).value_or("Meters");
            if (input == "Meters" || input == "Metres") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Meters;
            }
            if (input == "Kilometers") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Kilometers;
            }
            if (input == "Miles") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Miles;
            }
            if (input == "Inches") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Inches;
            }
        }
        const auto& argument1 = arrayMember(value, 1);
        if (isObject(argument1)) {
            style::conversion::Error error;
            auto geojson = toGeoJSON(argument1, error);
            if (geojson && error.message.empty()) {
                return Arguments(*geojson, unit);
            }
            ctx.error(error.message);
        }
    }
    ctx.error("'distance' expression needs to be an array with one/two arguments.");
    return nullopt;
}

optional<Feature::geometry_type> getGeometry(const Feature& feature, mbgl::style::expression::ParsingContext& ctx) {
    const auto type = apply_visitor(ToFeatureType(), feature.geometry);
    if (type == FeatureType::Point || type == FeatureType::LineString) {
        return feature.geometry;
    }
    ctx.error("'distance' expression requires valid geojson source that contains Point/LineString geometry type.");
    return nullopt;
}
} // namespace

namespace style {
namespace expression {

Distance::Distance(GeoJSON geojson, Feature::geometry_type geometries_, mapbox::cheap_ruler::CheapRuler::Unit unit_)
    : Expression(Kind::Distance, type::Number),
      geoJSONSource(std::move(geojson)),
      geometries(std::move(geometries_)),
      unit(unit_) {}

Distance::~Distance() = default;

using namespace mbgl::style::conversion;

EvaluationResult Distance::evaluate(const EvaluationContext& params) const {
    if (!params.feature || !params.canonical) {
        return EvaluationError{"distance expression requirs valid feature and canonical information."};
    }
    auto geometryType = params.feature->getType();
    if (geometryType == FeatureType::Point) {
        auto distance = calculateDistance(*params.feature, *params.canonical, geometries, unit);
        return distance;
    }
    if (geometryType == FeatureType::LineString) {
        auto distance = calculateDistance(*params.feature, *params.canonical, geometries, unit);
        return distance;
    }
    return EvaluationError{"distance expression currently only supports feature with Point geometry."};
}

ParseResult Distance::parse(const Convertible& value, ParsingContext& ctx) {
    auto parsedValue = parseValue(value, ctx);
    if (!parsedValue) {
        return ParseResult();
    }

    return parsedValue->geojson.match(
        [&parsedValue, &ctx](const mapbox::geometry::geometry<double>& geometrySet) {
            if (auto ret = getGeometry(mbgl::Feature(geometrySet), ctx)) {
                return ParseResult(
                    std::make_unique<Distance>(parsedValue->geojson, std::move(*ret), parsedValue->unit));
            }
            return ParseResult();
        },
        [&parsedValue, &ctx](const mapbox::feature::feature<double>& feature) {
            if (auto ret = getGeometry(mbgl::Feature(feature), ctx)) {
                return ParseResult(
                    std::make_unique<Distance>(parsedValue->geojson, std::move(*ret), parsedValue->unit));
            }
            return ParseResult();
        },
        [&parsedValue, &ctx](const mapbox::feature::feature_collection<double>& features) {
            for (const auto& feature : features) {
                if (auto ret = getGeometry(mbgl::Feature(feature), ctx)) {
                    return ParseResult(
                        std::make_unique<Distance>(parsedValue->geojson, std::move(*ret), parsedValue->unit));
                }
            }
            return ParseResult();
        },
        [&ctx](const auto&) {
            ctx.error("'distance' expression requires valid geojson that contains LineString/Point geometries.");
            return ParseResult();
        });

    return ParseResult();
}

Value convertValue(const mapbox::geojson::rapidjson_value& v) {
    if (v.IsDouble()) {
        return v.GetDouble();
    }
    if (v.IsString()) {
        return std::string(v.GetString());
    }
    if (v.IsArray()) {
        std::vector<Value> result;
        result.reserve(v.Size());
        for (const auto& m : v.GetArray()) {
            result.push_back(convertValue(m));
        }
        return result;
    }
    if (v.IsObject()) {
        std::unordered_map<std::string, Value> result;
        for (const auto& m : v.GetObject()) {
            result.emplace(m.name.GetString(), convertValue(m.value));
        }
        return result;
    }
    // Ignore other types as valid geojson only contains above types.
    return Null;
}

mbgl::Value Distance::serialize() const {
    std::unordered_map<std::string, Value> serialized;
    rapidjson::CrtAllocator allocator;
    const mapbox::geojson::rapidjson_value value = mapbox::geojson::convert(geoJSONSource, allocator);
    if (value.IsObject()) {
        for (const auto& m : value.GetObject()) {
            serialized.emplace(m.name.GetString(), convertValue(m.value));
        }
    } else {
        mbgl::Log::Error(mbgl::Event::General,
                         "Failed to serialize 'distance' expression, converted rapidJSON is not an object");
    }
    return std::vector<mbgl::Value>{{getOperator(), *fromExpressionValue<mbgl::Value>(serialized)}};
}

bool Distance::operator==(const Expression& e) const {
    if (e.getKind() == Kind::Distance) {
        auto rhs = static_cast<const Distance*>(&e);
        return geoJSONSource == rhs->geoJSONSource && geometries == rhs->geometries && unit == rhs->unit;
    }
    return false;
}

std::vector<optional<Value>> Distance::possibleOutputs() const {
    return {nullopt};
}

std::string Distance::getOperator() const {
    return "distance";
}

} // namespace expression
} // namespace style
} // namespace mbgl
