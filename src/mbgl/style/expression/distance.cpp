#include <mbgl/style/expression/distance.hpp>

#include <mapbox/cheap_ruler.hpp>
#include <mapbox/geojson.hpp>
#include <mapbox/geometry.hpp>

#include <mbgl/style/conversion/json.hpp>
#include <mbgl/tile/geometry_tile_data.hpp>
#include <mbgl/util/geometry_util.hpp>
#include <mbgl/util/logging.hpp>
#include <mbgl/util/string.hpp>

#include <rapidjson/document.h>

#include <algorithm>
#include <deque>
#include <limits>
#include <queue>
#include <tuple>

namespace mbgl {
namespace {

using BBox = std::array<double, 4>;

// bbox[minX, minY, maxX, maxY]
void updateBBox(BBox& bbox, const mapbox::geometry::point<double>& p) {
    bbox[0] = std::min(p.x, bbox[0]);
    bbox[1] = std::min(p.y, bbox[1]);
    bbox[2] = std::max(p.x, bbox[2]);
    bbox[3] = std::max(p.y, bbox[3]);
}

void updateGeoSetBBox(BBox& bbox, const Feature::geometry_type& geoSet) {
    geoSet.match([&bbox](const mapbox::geometry::point<double>& point) { updateBBox(bbox, point); },
                 [&bbox](const mapbox::geometry::multi_point<double>& points) {
                     for (const auto& p : points) {
                         updateGeoSetBBox(bbox, p);
                     }
                 },
                 [&bbox](const mapbox::geometry::line_string<double>& line) {
                     for (const auto& p : line) {
                         updateGeoSetBBox(bbox, p);
                     }
                 },
                 [&bbox](const mapbox::geometry::multi_line_string<double>& lines) {
                     for (const auto& line : lines) {
                         updateGeoSetBBox(bbox, line);
                     }
                 },
                 [&bbox](const mapbox::geometry::polygon<double>& polygon) {
                     for (const auto& ring : polygon) {
                         for (const auto& p : ring) {
                             updateGeoSetBBox(bbox, p);
                         }
                     }
                 },
                 [&bbox](const mapbox::geometry::multi_polygon<double>& polygons) {
                     for (const auto& polygon : polygons) {
                         updateGeoSetBBox(bbox, polygon);
                     }
                 },
                 [](const auto&) {});
}

BBox getBBox(const Feature::geometry_type& geoSet) {
    BBox bbox{std::numeric_limits<double>::infinity(),
              std::numeric_limits<double>::infinity(),
              -std::numeric_limits<double>::infinity(),
              -std::numeric_limits<double>::infinity()};
    updateGeoSetBBox(bbox, geoSet);
    return bbox;
}

// Calculate the distance between two bounding boxes.
// Calculate the delta in x and y direction, and use two fake points {0, 0} and {dx, dy} to calculate the distance.
// Distance will be 0 if bounding box are overlapping.
double bboxToBBoxDistance(const BBox& bbox1, const BBox& bbox2, mapbox::cheap_ruler::CheapRuler& ruler) {
    double dx = 0.;
    double dy = 0.;
    // bbox1 in left side
    if (bbox1[2] < bbox2[0]) {
        dx = bbox2[0] - bbox1[2];
    }
    // bbox1 in right side
    if (bbox1[0] > bbox2[2]) {
        dx = bbox1[0] - bbox2[2];
    }
    // bbox1 in above side
    if (bbox1[1] > bbox2[3]) {
        dy = bbox1[1] - bbox2[3];
    }
    // bbox1 in down side
    if (bbox1[3] < bbox2[1]) {
        dy = bbox2[1] - bbox1[3];
    }
    return ruler.distance(mapbox::geometry::point<double>{0., 0.}, mapbox::geometry::point<double>{dx, dy});
}

// geoSet1 could be any geometry type, the bounding box is calculated based on the whole geometry sets
// geoSet2 is required to be multi_line_string or multi_polygon, the bounding box is calculated based on
// each line_string or polygon element belong to it.
// The purpose of this function is to filter out line_string or polygon that are far away from geoSet1,
// based on the bounding box distance.
std::vector<std::size_t> getFilteredIndexes(const Feature::geometry_type& geoSet1,
                                            const Feature::geometry_type& geoSet2,
                                            mapbox::cheap_ruler::CheapRuler& ruler) {
    auto bbox1 = getBBox(geoSet1);
    std::vector<std::pair<double, std::size_t>> distanceList;
    geoSet2.match(
        [&bbox1, &distanceList, &ruler](const mapbox::geometry::multi_line_string<double>& lines) {
            distanceList.reserve(lines.size());
            for (std::size_t i = 0; i < lines.size(); ++i) {
                auto bbox2 = getBBox(lines[i]);
                auto pair = std::make_pair(bboxToBBoxDistance(bbox1, bbox2, ruler), i);
                distanceList.insert(std::upper_bound(distanceList.begin(),
                                                     distanceList.end(),
                                                     pair,
                                                     [](const std::pair<double, std::size_t>& left,
                                                        const std::pair<double, std::size_t>& right) {
                                                         return left.first < right.first;
                                                     }),
                                    std::move(pair));
            }
        },
        [&bbox1, &distanceList, &ruler](const mapbox::geometry::multi_polygon<double>& polygons) {
            distanceList.reserve(polygons.size());
            for (std::size_t i = 0; i < polygons.size(); ++i) {
                auto bbox2 = getBBox(polygons[i]);
                auto pair = std::make_pair(bboxToBBoxDistance(bbox1, bbox2, ruler), i);
                distanceList.insert(std::upper_bound(distanceList.begin(),
                                                     distanceList.end(),
                                                     pair,
                                                     [](const std::pair<double, std::size_t>& left,
                                                        const std::pair<double, std::size_t>& right) {
                                                         return left.first < right.first;
                                                     }),
                                    std::move(pair));
            }
        },
        [](const auto&) {});

    std::vector<std::size_t> ret;
    ret.reserve(distanceList.size());
    for (auto iter = distanceList.begin(); iter != distanceList.end(); ++iter) {
        if (iter == distanceList.begin() || iter->first == 0.) {
            ret.emplace_back(iter->second);
        }
    }
    return ret;
}

double pointToLineDistance(const mapbox::geometry::point<double>& point,
                           const mapbox::geometry::line_string<double>& line,
                           mapbox::cheap_ruler::CheapRuler& ruler) {
    const auto nearestPoint = std::get<0>(ruler.pointOnLine(line, point));
    return ruler.distance(point, nearestPoint);
}

double pointToPointsDistance(const mapbox::geometry::point<double>& point,
                             const mapbox::geometry::multi_point<double>& points,
                             mapbox::cheap_ruler::CheapRuler& ruler) {
    double dist = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < points.size(); ++i) {
        dist = std::min(dist, ruler.distance(point, points[i]));
        if (dist == 0.) return dist;
    }
    return dist;
}

double pointsToPointsDistance(const mapbox::geometry::multi_point<double>& points1,
                              const mapbox::geometry::multi_point<double>& points2,
                              mapbox::cheap_ruler::CheapRuler& ruler) {
    double dist = std::numeric_limits<double>::infinity();
    for (const auto& point : points1) {
        for (size_t i = 0; i < points2.size(); ++i) {
            dist = std::min(dist, ruler.distance(point, points2[i]));
            if (dist == 0.) return dist;
        }
    }
    return dist;
}

double lineToLineDistance(const mapbox::geometry::line_string<double>& line1,
                          const mapbox::geometry::line_string<double>& line2,
                          mapbox::cheap_ruler::CheapRuler& ruler) {
    double dist = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < line1.size() - 1; ++i) {
        const auto& p1 = line1[i];
        const auto& p2 = line1[i + 1];
        for (std::size_t j = 0; j < line2.size() - 1; ++j) {
            const auto& q1 = line2[j];
            const auto& q2 = line2[j + 1];

            if (GeometryUtil<double>::segmentIntersectSegment(p1, p2, q1, q2)) return 0.;
            dist = std::min(dist, pointToLineDistance(p1, mapbox::geometry::line_string<double>{q1, q2}, ruler));
            dist = std::min(dist, pointToLineDistance(p2, mapbox::geometry::line_string<double>{q1, q2}, ruler));
            dist = std::min(dist, pointToLineDistance(q1, mapbox::geometry::line_string<double>{p1, p2}, ruler));
            dist = std::min(dist, pointToLineDistance(q1, mapbox::geometry::line_string<double>{p1, p2}, ruler));
        }
    }
    return dist;
}

// Divide and conqure, the time complexity is O(n*lgn), faster than Brute force O(n*n)
// However it requires extra space, it is a trade off between space and time.
double pointSetsDistance(const mapbox::geometry::multi_point<double>& pointSet1,
                         bool set1IsLine,
                         const mapbox::geometry::multi_point<double>& pointSet2,
                         bool set2IsLine,
                         mapbox::cheap_ruler::CheapRuler& ruler) {
    using DistPair = std::tuple<double, mapbox::geometry::multi_point<double>, mapbox::geometry::multi_point<double>>;
    class Comparator {
    public:
        bool operator()(DistPair& left, DistPair& right) { return std::get<0>(left) < std::get<0>(right); }
    };
    // The priority queue will ensure the top element would always be the pair that has the biggest distance
    using DistQueue = std::priority_queue<DistPair, std::deque<DistPair>, Comparator>;

    auto miniDist = ruler.distance(pointSet1[0], pointSet2[0]);
    DistQueue distQueue;
    distQueue.push(std::make_tuple(0, pointSet1, pointSet2));

    while (!distQueue.empty()) {
        auto& distPair = distQueue.top();
        distQueue.pop();
        if (std::get<0>(distPair) > miniDist) break;
        auto& pSetA = std::get<1>(distPair);
        auto& pSetB = std::get<2>(distPair);
        static const std::size_t MinClusterSize = 20;
        // In case the set size are relatively small, we could use brute-force directly
        if (pSetA.size() <= MinClusterSize || pSetB.size() <= MinClusterSize) {
            if (set1IsLine && set2IsLine) {
                miniDist = std::min(miniDist, lineToLineDistance(pSetA, pSetB, ruler));
                if (miniDist == 0.) return 0.;
            } else if (!set1IsLine && !set2IsLine) {
                miniDist = std::min(miniDist, pointsToPointsDistance(pSetA, pSetB, ruler));
                if (miniDist == 0.) return 0.;
            } else if (!set1IsLine && set2IsLine) {
                for (const auto& p : pSetA) {
                    miniDist = std::min(miniDist, pointToLineDistance(p, pSetB, ruler));
                    if (miniDist == 0.) return 0.;
                }
            } else if (set1IsLine && !set2IsLine) {
                for (const auto& p : pSetB) {
                    miniDist = std::min(miniDist, pointToLineDistance(p, pSetA, ruler));
                    if (miniDist == 0.) return 0.;
                }
            }
        } else {
            auto size1 = pSetA.size() / 2;
            auto pSetA1 = mapbox::geometry::multi_point<double>(pSetA.begin(), pSetA.begin() + size1);
            // If set is a line string the cutting point needs to be taken into both of the evolved sets
            auto pSetA2 = mapbox::geometry::multi_point<double>(pSetA.begin() + size1 - set1IsLine, pSetA.end());

            auto size2 = pSetB.size() / 2;
            auto pSetB1 = mapbox::geometry::multi_point<double>(pSetB.begin(), pSetB.begin() + size2);
            auto pSetB2 = mapbox::geometry::multi_point<double>(pSetB.begin() + size2 - set2IsLine, pSetB.end());

            const auto updateQueue = [&distQueue, &miniDist, &ruler](
                                         const mapbox::geometry::multi_point<double>& set1,
                                         const mapbox::geometry::multi_point<double>& set2) {
                auto tempDist = bboxToBBoxDistance(getBBox(set1), getBBox(set2), ruler);
                // Insert new pair to the queue if the bbox distance is less or equal to miniDist,
                // otherwise the pair could be abandoned as its smallest distance is bigger than mininDist.
                // The pair with biggest distance will be at the top
                if (tempDist <= miniDist) distQueue.push(std::make_tuple(tempDist, set1, set2));
            };

            updateQueue(pSetA1, pSetB1);
            updateQueue(pSetA1, pSetB2);
            updateQueue(pSetA2, pSetB1);
            updateQueue(pSetA2, pSetB2);
        }
    }
    return miniDist;
}

double pointsToLinesDistance(const mapbox::geometry::multi_point<double>& points,
                             const mapbox::geometry::multi_line_string<double>& lines,
                             mapbox::cheap_ruler::CheapRuler& ruler) {
    auto indexes = getFilteredIndexes(points, lines, ruler);
    double dist = std::numeric_limits<double>::infinity();
    for (const auto& index : indexes) {
        dist = std::min(dist,
                        pointSetsDistance(points, false /*isLineString*/, lines[index], true /*isLineString*/, ruler));
        if (dist == 0.) return dist;
    }
    return dist;
}

double lineToLinesDistance(const mapbox::geometry::line_string<double>& line,
                           const mapbox::geometry::multi_line_string<double>& lines,
                           mapbox::cheap_ruler::CheapRuler& ruler) {
    auto indexes = getFilteredIndexes(line, lines, ruler);
    double dist = std::numeric_limits<double>::infinity();
    for (const auto& index : indexes) {
        dist =
            std::min(dist, pointSetsDistance(line, true /*isLineString*/, lines[index], true /*isLineString*/, ruler));
        if (dist == 0.) return dist;
    }
    return dist;
}

double pointsToGeometryDistance(const mapbox::geometry::multi_point<double>& points,
                                const Feature::geometry_type& geoSet,
                                mapbox::cheap_ruler::CheapRuler::Unit unit) {
    mapbox::cheap_ruler::CheapRuler ruler(points.front().y, unit);
    return geoSet.match(
        [&points, &ruler](const mapbox::geometry::point<double>& p) { return pointToPointsDistance(p, points, ruler); },
        [&points, &ruler](const mapbox::geometry::multi_point<double>& points1) {
            return pointSetsDistance(points, false /*isLineString*/, points1, false /*isLineString*/, ruler);
        },
        [&points, &ruler](const mapbox::geometry::line_string<double>& line) {
            return pointSetsDistance(points, false /*isLineString*/, line, true /*isLineString*/, ruler);
        },
        [&points, &ruler](const mapbox::geometry::multi_line_string<double>& lines) {
            return pointsToLinesDistance(points, lines, ruler);
        },
        [](const auto&) { return std::numeric_limits<double>::quiet_NaN(); });
}

double lineToGeometryDistance(const mapbox::geometry::line_string<double>& line,
                              const Feature::geometry_type& geoSet,
                              mapbox::cheap_ruler::CheapRuler::Unit unit) {
    assert(!line.empty());
    mapbox::cheap_ruler::CheapRuler ruler(line.front().y, unit);
    return geoSet.match(
        [&line, &ruler](const mapbox::geometry::point<double>& p) { return pointToLineDistance(p, line, ruler); },
        [&line, &ruler](const mapbox::geometry::multi_point<double>& points) {
            return pointSetsDistance(points, false /*isLineString*/, line, true /*isLineString*/, ruler);
        },
        [&line, &ruler](const mapbox::geometry::line_string<double>& line1) {
            return pointSetsDistance(line, true /*isLineString*/, line1, true /*isLineString*/, ruler);
        },
        [&line, &ruler](const mapbox::geometry::multi_line_string<double>& lines) {
            return lineToLinesDistance(line, lines, ruler);
        },
        [](const auto&) { return std::numeric_limits<double>::quiet_NaN(); });
}

double calculateDistance(const GeometryTileFeature& feature,
                         const CanonicalTileID& canonical,
                         const Feature::geometry_type& geoSet,
                         mapbox::cheap_ruler::CheapRuler::Unit unit) {
    return convertGeometry(feature, canonical)
        .match(
            [&geoSet, &unit](const mapbox::geometry::point<double>& point) -> double {
                return pointsToGeometryDistance(mapbox::geometry::multi_point<double>{point}, geoSet, unit);
            },
            [&geoSet, &unit](const mapbox::geometry::multi_point<double>& points) -> double {
                return pointsToGeometryDistance(points, geoSet, unit);
            },
            [&geoSet, &unit](const mapbox::geometry::line_string<double>& line) -> double {
                return lineToGeometryDistance(line, geoSet, unit);
            },
            [&geoSet, &unit](const mapbox::geometry::multi_line_string<double>& lines) -> double {
                double dist = std::numeric_limits<double>::infinity();
                for (const auto& line : lines) {
                    dist = std::min(dist, lineToGeometryDistance(line, geoSet, unit));
                    if (dist == 0.) return dist;
                }
                return dist;
            },
            [](const auto&) -> double { return std::numeric_limits<double>::quiet_NaN(); });
}

struct Arguments {
    Arguments(GeoJSON& geojson_, mapbox::cheap_ruler::CheapRuler::Unit unit_)
        : geojson(std::move(geojson_)), unit(unit_) {}

    GeoJSON geojson;
    mapbox::cheap_ruler::CheapRuler::Unit unit;
};

optional<Arguments> parseValue(const style::conversion::Convertible& value, style::expression::ParsingContext& ctx) {
    if (isArray(value)) {
        // object value, quoted with ["Distance", GeoJSONObj, "unit(optional)"]
        auto length = arrayLength(value);
        if (length != 2 && length != 3) {
            ctx.error("'distance' expression requires either one argument or two arguments, but found " +
                      util::toString(arrayLength(value) - 1) + " instead.");
            return nullopt;
        }

        // Parse Unit info for distance calculation, "Meters" by default
        mapbox::cheap_ruler::CheapRuler::Unit unit = mapbox::cheap_ruler::CheapRuler::Unit::Meters;
        if (length == 3) {
            auto input = toString(arrayMember(value, 2));
            if (input == nullopt) {
                ctx.error("Failed to parse unit argument from 'distance' expression");
                return nullopt;
            }
            if (*input == "Meters" || *input == "Metres") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Meters;
            } else if (*input == "Kilometers") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Kilometers;
            } else if (*input == "Miles") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Miles;
            } else if (*input == "NauticalMiles") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::NauticalMiles;
            } else if (*input == "Yards") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Kilometers;
            } else if (*input == "Feet") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Miles;
            } else if (*input == "Inches") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Inches;
            } else {
                ctx.error(
                    "'distance' expression only accepts following Units:  'Kilometers', 'Miles', 'NauticalMiles', "
                    "'Meters', 'Metres', 'Yards', 'Feet', 'Inches'.");
                return nullopt;
            }
        }
        // Parse geometry info
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
    ctx.error(
        "'distance' expression needs to be an array with format [\"Distance\", GeoJSONObj, \"unit(optional, Meters by "
        "default)\"].");
    return nullopt;
}

optional<Feature::geometry_type> getGeometry(const Feature& feature, mbgl::style::expression::ParsingContext& ctx) {
    const auto type = apply_visitor(ToFeatureType(), feature.geometry);
    if (type == FeatureType::Point || type == FeatureType::LineString) {
        return feature.geometry;
    }
    ctx.error("'distance' expression requires valid geojson object with valid geometry type: Point/LineString.");
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
    if (geometryType == FeatureType::Point || geometryType == FeatureType::LineString) {
        auto distance = calculateDistance(*params.feature, *params.canonical, geometries, unit);
        if (!std::isnan(distance)) {
            return distance;
        }
    }
    return EvaluationError{"distance expression currently only evaluates Point/LineString geometries."};
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
