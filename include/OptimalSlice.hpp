/*
 * This is the implement slice algorithm from Rodrigo, 2017 (An Optimal Algorithm for 3D Triangle Mesh Slicing)
 * which is claimed to be faster than slic3r and CGAL method.
 */
#pragma once
#ifndef SLICE_PRECISION
// There's a problem of double hashing with the precision less than 1e-8 (e.g. 1e-10)
// when performed contour constructing
#define SLICE_PRECISION 1e-8
#endif
#define DOUBLE_EQ(x,y) (abs(x - y) < SLICE_PRECISION)
#define DOUBLE_GT(x,y) ((x - y) > SLICE_PRECISION)
#define DOUBLE_LT(x,y) ((y - x) > SLICE_PRECISION)
#define DOUBLE_GT_EQ(x,y) (DOUBLE_EQ(x,y) || DOUBLE_GT(x,y))
#define DOUBLE_LT_EQ(x,y) (DOUBLE_EQ(x,y) || DOUBLE_LT(x,y))
#define USE_PARALLEL

#include "Mesh.h"
#include "tbb/tbb.h"

namespace slice {
    
    enum Direction {X = 0,Y,Z};
    enum PolygonSide {OUTSIDE = 0, INSIDE};

    struct Point2d {
        double x;
        double y;
        bool operator==(const Point2d& ls) const {
            return (DOUBLE_EQ(x, ls.x) && DOUBLE_EQ(y, ls.y));
        }
        bool operator< (const Point2d& ls) const {
            if (x < ls.x) return true;
            return y < ls.y;
        }
    };
    
    struct Point3d { 
        double x;
        double y;
        double z;
        bool operator==(const Point3d& ls) const {
            return (DOUBLE_EQ(x, ls.x) && DOUBLE_EQ(y, ls.y) && DOUBLE_EQ(z, ls.z));
        }
        bool operator< (const Point3d& ls) const {
            if (x < ls.x) return true;
            else if DOUBLE_EQ(x, ls.x)
                if (y < ls.y) return true;
                else if (DOUBLE_EQ(y, ls.y)) {
                    return (z < ls.z);
                }
            return false;
        }
    };

    Point2d make_point2d(Point3d p, int direction = Direction::Z) {
        if (direction == Direction::X) return Point2d{ p.y, p.z };
        else if (direction == Direction::Y) return Point2d{ p.x, p.z };
        return Point2d{ p.x, p.y };
    }

    class Triangle {
    public:
        Point3d v[3];
        double min[3], max[3];
        Triangle(TMesh::FaceType face) {
            assert(face.VN() == 3);
            vcg::Point3d point = face.P(0);
            v[0].x = point[0];
            v[0].y = point[1];
            v[0].z = point[2];
            point = face.P(1);
            v[1].x = point[0];
            v[1].y = point[1];
            v[1].z = point[2];
            point = face.P(2);
            v[2].x = point[0];
            v[2].y = point[1];
            v[2].z = point[2];
            updateMinMax();
        }
        void updateMinMax() {
            min[0] = std::min({ v[0].x, v[1].x, v[2].x });
            min[1] = std::min({ v[0].y, v[1].y, v[2].y });
            min[2] = std::min({ v[0].z, v[1].z, v[2].z });
            max[0] = std::max({ v[0].x, v[1].x, v[2].x });
            max[1] = std::max({ v[0].y, v[1].y, v[2].y });
            max[2] = std::max({ v[0].z, v[1].z, v[2].z });
        }

        double minOf(int i, int j, int direction) {
            assert(direction >= 0 && direction <= 2);
            assert(i >= 0 && i <= 2);
            assert(j >= 0 && j <= 2);
            if (direction == 0)
                return v[i].x < v[j].x ? v[i].x : v[j].x;
            else if (direction == 1)
                return v[i].y < v[j].y ? v[i].y : v[j].y;
            else
                return v[i].z < v[j].z ? v[i].z : v[j].z;
        }
        double maxOf(int i, int j, int direction) {
            assert(direction >= 0 && direction <= 2);
            assert(i >= 0 && i <= 2);
            assert(j >= 0 && j <= 2);
            if (direction == 0)
                return v[i].x < v[j].x ? v[j].x : v[i].x;
            else if (direction == 1)
                return v[i].y < v[j].y ? v[j].y : v[i].y;
            else
                return v[i].z < v[j].z ? v[j].z : v[i].z;
        }

        bool operator== (const Triangle& ls) const {
            return (v[0] == ls.v[0]) && (v[1] == ls.v[1]) && (v[2] == ls.v[2]);
        }

        bool operator< (const Triangle& ls) const {
            if (v[0] < ls.v[0]) return true;
            else if (v[0] == ls.v[0])
                if (v[1] < ls.v[1]) return true;
                else if (v[1] == ls.v[1]) {
                    return (v[2] < ls.v[2]);
                }
            return false;
        }
    };

    class Line {
    public:
        Point3d v[2];

        Line() {}

        Line(Point3d v0, Point3d v1, size_t index) {
            v[0] = v0;
            v[1] = v1;
            sort();
        }

        void sort() {
            if (v[1] < v[0]) std::swap(v[0], v[1]);
        }

        bool operator== (const Line& ls) const {
            return (v[0] == ls.v[0]) && (v[1] == ls.v[1]);
        }

        bool operator< (const Line& ls) const {
            if (v[0] < ls.v[0]) return true;
            else if (v[0] == ls.v[0]) return (v[1] < ls.v[1]);
            return false;
        }
    };

    class SupportLine2d {
    public:
        double x, y;
        double m, theta, c;
        bool is_vertical = false;
        SupportLine2d(double _x, double _y, double _theta) : x(_x), y(_y), theta(_theta) {
            update();
        }
        SupportLine2d(Point2d p, double _theta) : x(p.x), y(p.y), theta(_theta) {
            update();
        }
        void update() {
            // adjust theta to be in range [-90, 90]
            if (theta > 90) theta -= 180;
            else if (theta < -90) theta += 180;
            if (DOUBLE_EQ(abs(theta), 90)) {
                is_vertical = true;
                m = 0;
                c = x;
            }
            else {
                is_vertical = false;
                m = tan(theta * M_PI / 180);
                c = y - m * x;
            }
        }
        void update(double _x, double _y) {
            x = _x;
            y = _y;
            update();
        }

        void update(const Point2d& p) {
            update(p.x, p.y);
        }
        
        SupportLine2d& operator+= (const double theta) {
            this->theta += theta;
            update();
            return *this;
        }
    };

    double two_line_distance(const SupportLine2d& line1, const SupportLine2d &line2) {
        if (DOUBLE_EQ(line1.theta, line2.theta) && DOUBLE_EQ(line1.m, line2.m) && line1.is_vertical == line2.is_vertical) {
            return abs(line2.c - line1.c) / sqrt(line1.m * line1.m + 1);
        }
        return 0;
    }

    double line_point_angle(const SupportLine2d& line, const Point2d& p) {
        double angle;
        if (DOUBLE_EQ(p.x, line.x)) angle = 90;
        else if (DOUBLE_EQ(p.y, line.y)) angle = 0;
        else angle = atan((p.y - line.y) / (p.x - line.x)) * 180 / M_PI;
        // calculate delta
        if (DOUBLE_GT_EQ(line.theta, angle)) angle = line.theta - angle;
        else angle = 180 + line.theta - angle;
        // normalized angle
        if (angle > 90) angle -= 180;
        return angle;
    }

    class FeretDiameter {
    public:
        bool empty = true;
        double min;
        double max;
        double perpendicularMax;
        double perimeter;
        double angleMin;
        double angleMax;
        FeretDiameter() : min(0), max(0), perpendicularMax(0), perimeter(0), angleMin(0), angleMax(0) {}
        FeretDiameter(const SupportLine2d (&lines)[4]) {
            double d1 = two_line_distance(lines[0], lines[2]);
            double d2 = two_line_distance(lines[1], lines[3]);
            if (d1 > d2) {
                min = d2; angleMin = lines[1].theta;
                max = d1; angleMax = lines[0].theta;
            }
            else {
                min = d1; angleMin = lines[0].theta;
                max = d2; angleMax = lines[1].theta;
            }
            perpendicularMax = max;
            empty = false;
        }
        void update(const SupportLine2d (&lines)[4]) {
            double d1 = two_line_distance(lines[0], lines[2]);
            double d2 = two_line_distance(lines[1], lines[3]);
            double angle1 = lines[0].theta;
            double angle2 = lines[1].theta;
            if (d1 > d2) {
                std::swap(d1, d2);
                std::swap(angle1, angle2);
            }
            if (min > d1) {
                min = d1;
                angleMin = angle1;
                perpendicularMax = d2;
            }
            if (max < d1) {
                max = d1;
                angleMax = angle1;
            }
            if (max < d2) {
                max = d2;
                angleMax = angle2;
            }
        }
    };

    typedef std::vector<double> Plane;
    typedef std::vector<Triangle> Triangles;
    typedef std::vector<Triangles> Layer;
    typedef std::vector<Line> Lines;
    typedef std::vector<Lines> Slice;
    typedef std::vector<Point2d> Polygon;
    typedef std::vector<Polygon> Polygons;
    typedef std::vector<Polygons> ContourSlice;
    typedef std::pair<Point2d, Point2d> PairPoint2d;
    typedef std::unordered_map<Point2d, PairPoint2d> ContourHash;
    typedef std::vector<PolygonSide> PolygonSides;

    std::ostream& operator<< (std::ostream& out, Point2d const& data) {
        out << "[" << data.x << "," << data.y << "]";
        return out;
    }

    std::ostream& operator<< (std::ostream& out, PairPoint2d const& data) {
        out << "(" << data.first << " " << data.second << ")";
        return out;
    }

    std::ostream& operator<< (std::ostream& out, Point3d const& data) {
        out << "[" << data.x << "," << data.y << "," << data.z << "]";
        return out;
    }

    std::ostream& operator<< (std::ostream& out, Triangle const& data) {
        out << "slice::Triangle(" << data.v[0] << " " << data.v[1] << " " << data.v[2] << ")";
        return out;
    }

    std::ostream& operator<< (std::ostream& out, Line const& data) {
        out << "slice::Line(" << data.v[0] << " " << data.v[1] << ")";
        return out;
    }
    
    std::ostream& operator<< (std::ostream& out, FeretDiameter const& data) {
        out << "slice::FeretDiameter(min:" << data.min << ", max:" << data.max << ", pmax:" << data.perpendicularMax << ")";
        return out;
    }

    std::ostream& operator<< (std::ostream& out, SupportLine2d const& data) {
        out << "slice::SupportLine2d(x:" << data.x << ", y:" << data.y << ", m:" << data.m << ", c:" << data.c << ", theta:" << data.theta << "[" << data.is_vertical <<"])";
        return out;
    }
}

namespace std {
    template<> struct hash<slice::Point2d> {
        size_t operator()(const slice::Point2d& p) const noexcept {
            size_t x = hash<long long>()(llround(p.x/SLICE_PRECISION));
            size_t y = hash<long long>()(llround(p.y/SLICE_PRECISION));
            return x ^ (y << 1);
        }
    };
}

namespace slice {
    
    inline void build_triangle_list(TMesh& mesh, size_t grid_size, Plane& P, Layer& L, int direction = Direction::Z) {
        // Uniform slicing with delta > 0
        // in this case, grid_size = k from Rodrigo paper
        assert(grid_size > 1 && direction <= 2 && direction >= 0);
        vcg::tri::UpdateBounding<TMesh>::Box(mesh);
        vcg::Box3d bbox = mesh.bbox;
        double minBBox = bbox.min[direction];
        double maxBBox = bbox.max[direction];
        vcg::Point3d dim = bbox.Dim();
        double delta = dim[direction] / (grid_size - 1);
        // build Plane vector P[0...k+1]
        P.resize(grid_size + 2);
        P[0] = minBBox - 10 * delta;
        P[1] = minBBox;
        P[grid_size + 1] = maxBBox + 10 * delta;
        for (size_t i = 2; i <= grid_size; i++)
            P[i] = P[i - 1] + delta;
        // initialize layer L[0...k+1]
        L.resize(grid_size + 2);
        for (size_t i = 0; i <= grid_size + 1; i++) L[i].clear();
        // foreach triangle in mesh
        // TODO: implement parallel version with tbb
        for (TMesh::FaceIterator it = mesh.face.begin(); it != mesh.face.end(); it++) {
            if (!it->IsD())
            {
                Triangle triangle(*it);
                size_t i = 0;
                i = size_t(ceil((triangle.min[direction] - P[1]) / delta) + 1);
                assert(i > 0 && i <= grid_size + 1);
                L[i].push_back(triangle);
            }
        }
    }

    inline Point3d compute_point_at_plane(Point3d v0, Point3d v1, double position, int direction = Direction::Z) {
        double dx = v1.x - v0.x;
        double dy = v1.y - v0.y;
        double dz = v1.z - v0.z;
        if (direction == 2) {
            assert(dz != 0);
            double frac = (position - v0.z) / dz;
            double x = frac * dx + v0.x;
            double y = frac * dy + v0.y;
            return Point3d{ x, y, position };
        }
        else if (direction == 1) {
            assert(dy != 0);
            double frac = (position - v0.y) / dy;
            double x = frac * dx + v0.x;
            double z = frac * dz + v0.z;
            return Point3d{ x, position, z };
        }
        else {
            assert(dx != 0);
            double frac = (position - v0.x) / dx;
            double y = frac * dy + v0.y;
            double z = frac * dz + v0.z;
            return Point3d{ position, y, z };
        }
    }

    // Modified version of Rodrigo (2017) and Adnan's slicing algorithm (2018)
    // (Real-time slicing algorithm for Stereolithography (STL) CAD model applied in additive manufacturing industry)
    // The ill-conditioned case will be cured
    bool compute_intersection(Triangle t, double position, Line& L, int direction = Direction::Z) {
        assert(direction >= 0 && direction <= 2);
        assert(t.min[direction] <= position && t.max[direction] >= position);
        int np = 0; // number of endpoints on the plane
        std::vector<int> found_indexs;
        found_indexs.reserve(3);
        for (int i = 0; i < 3; i++) {
            if ((direction == Direction::X && DOUBLE_EQ(t.v[i].x, position)) ||
                (direction == Direction::Y && DOUBLE_EQ(t.v[i].y, position)) ||
                (direction == Direction::Z && DOUBLE_EQ(t.v[i].z, position))) {
                np++;
                found_indexs.push_back(i);
            }
        }
        if (np == 0) {
            int k = 0;
            for (int i = 0; i < 3; i++) {
                int next_i = (i == 2) ? 0 : i + 1;
                double min = t.minOf(i, next_i, direction);
                double max = t.maxOf(i, next_i, direction);
                if (min <= position && max >= position) {
                    assert(k < 2);
                    L.v[k] = compute_point_at_plane(t.v[i], t.v[next_i], position, direction);
                    k++;
                }
            }
            assert(k == 2);
            L.sort();
            return true;
        }
        else if (np == 1 && DOUBLE_GT(t.max[direction], position) && DOUBLE_LT(t.min[direction], position)) {
            assert(found_indexs.size() == 1);
            int i = (found_indexs[0] + 1) % 3;
            int next_i = (i + 1) % 3;
            L.v[0] = t.v[found_indexs[0]];
            L.v[1] = compute_point_at_plane(t.v[i], t.v[next_i], position, direction);
            L.sort();
            return true;
        }
        else if (np == 2) {
            assert(found_indexs.size() == 2);
            L.v[0] = t.v[found_indexs[0]];
            L.v[1] = t.v[found_indexs[1]];
            L.sort();
            return true;
        }
        return false;
    }

    Slice incremental_slicing(TMesh& mesh, size_t grid_size, int direction = Direction::Z) {
        slice::Plane P;
        slice::Layer L;
        slice::build_triangle_list(mesh, grid_size, P, L, direction);
        Slice S(grid_size);
        
        Triangles A;
        for (size_t i = 1; i <= grid_size; i++) {
            if (L[i].size() > 0) {
                A.reserve(A.size() + L[i].size());
                A.insert(A.end(), L[i].begin(), L[i].end());
            }
            S[i - 1].clear();
#ifndef USE_PARALLEL
            for (Triangles::iterator t = A.begin(); t != A.end();) {
                if (t->max[direction] < P[i]) {
                    t = A.erase(t);
                }
                else {
                    Line line;
                    if (t->max[direction] >= P[i] && t->min[direction] <= P[i]) {
                        if (compute_intersection(*t, P[i], line, direction)) {
                            S[i - 1].push_back(line);
                        }
                    }
                    t++;
                }
            }
#else
            tbb::spin_mutex concatMutex;
            tbb::spin_mutex writeMutex;
            static tbb::affinity_partitioner ap;
            std::vector<size_t> deleteIndex;
            deleteIndex.clear();
            
            Triangles new_A;
            new_A.reserve(A.size());
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, A.size()),
                [&](const tbb::blocked_range<size_t> &r)
                {
                    Triangles local_A;
                    std::copy(A.begin() + r.begin(), A.begin() + r.end(), std::back_inserter(local_A));
                    for (Triangles::iterator t = local_A.begin(); t != local_A.end();) {
                        if (t->max[direction] < P[i]) {
                            t = local_A.erase(t);
                        }
                        else {
                            Line line;
                            if (t->max[direction] >= P[i] && t->min[direction] <= P[i]) {
                                if (compute_intersection(*t, P[i], line, direction)) {
                                    tbb::spin_mutex::scoped_lock lock(writeMutex);
                                    S[i - 1].push_back(line);
                                }
                            }
                            t++;
                        }
                    }
                    if (local_A.size() > 0) {
                        tbb::spin_mutex::scoped_lock lock(concatMutex);
                        new_A.insert(new_A.end(), local_A.begin(), local_A.end());
                    }
                },
                ap
            );
            A = new_A;
#endif
        }
        return S;
    }

#ifndef USE_PARALLEL
    ContourSlice contour_construct(Slice const& S, int direction = Direction::Z) {
        ContourSlice CS(S.size());
        ContourHash hash;
        for (size_t i = 0, len = S.size(); i < len; i++) {
            CS[i].clear();
            hash.clear();
            hash.reserve(S[i].size() + 1);
            for (Lines::const_iterator l = S[i].begin(); l != S[i].end(); l++) {
                Point2d u = make_point2d(l->v[0], direction);
                Point2d v = make_point2d(l->v[1], direction);
                ContourHash::iterator item = hash.find(u);
                if (item == hash.end())
                    hash.emplace(u, make_pair(v, v));
                else {
                    Point2d w = item->second.first;
                    Point2d target = item->second.second;
                    if (w == target) {
                        item->second = make_pair(w, v);
                    }
                }
                item = hash.find(v);
                if (item == hash.end())
                    hash.emplace(v, make_pair(u, u));
                else {
                    Point2d w = item->second.first;
                    Point2d target = item->second.second;
                    if (w == target) {
                        item->second = make_pair(w, u);
                    }
                }
            }
            //std::cout << " [Hash OK: " << hash.size() << "] from S[i]: " << S[i].size() << std::endl;
            //std::copy(S[i].begin(), S[i].end(), std::ostream_iterator<Line>(std::cout, " "));
            //std::cout << std::endl;
            /* TODO: remove this debug message
            for (ContourHash::const_iterator item = hash.begin(); item != hash.end(); item++) {
                std::cout << "  " << item->first << ": " << item->second << (item->second.first == item->second.second ? " [***]" : "") << std::endl;
                //std::cout << "  " << item->first << ": " << std::hash<Point2d>()(item->first) << std::endl;
            }
            //*/
            while (!hash.empty()) {
                ContourHash::const_iterator item = hash.begin();
                assert(item != hash.end());
                Polygon C;
                C.push_back(item->first);
                C.push_back(item->second.first);
                Point2d last = item->second.second;
                hash.erase(item);
                for (size_t j = 1;; j++) {
                    item = hash.find(C[j]);
                    //if (item == hash.end()) break;
                    assert(item != hash.end());
                    if (!(C[j] == last)) {
                        if (item->second.first == C[j - 1])
                            C.push_back(item->second.second);
                        else
                            C.push_back(item->second.first);
                    }
                    hash.erase(item);
                    if (C[j] == last) break;
                }
                CS[i].push_back(C);
            }
        }
        return CS;
    }
#else
    ContourSlice contour_construct(Slice const& S, int direction = Direction::Z) {
        ContourSlice CS(S.size());
        static tbb::affinity_partitioner ap;
        tbb::spin_mutex printMutex;
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, S.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i < r.end(); i++) {
                    CS[i].clear();
                    ContourHash hash;
                    hash.clear();
                    hash.reserve(S[i].size() + 1);
                    for (Lines::const_iterator l = S[i].begin(); l != S[i].end(); l++) {
                        Point2d u = make_point2d(l->v[0], direction);
                        Point2d v = make_point2d(l->v[1], direction);
                        ContourHash::iterator item = hash.find(u);
                        if (item == hash.end())
                            hash.emplace(u, make_pair(v, v));
                        else {
                            Point2d w = item->second.first;
                            Point2d target = item->second.second;
                            if (w == target) {
                                item->second = make_pair(w, v);
                            }
                        }
                        item = hash.find(v);
                        if (item == hash.end())
                            hash.emplace(v, make_pair(u, u));
                        else {
                            Point2d w = item->second.first;
                            Point2d target = item->second.second;
                            if (w == target) {
                                item->second = make_pair(w, u);
                            }
                        }
                    }
                    while (!hash.empty()) {
                        ContourHash::const_iterator item = hash.begin();
                        assert(item != hash.end());
                        Polygon C;
                        C.push_back(item->first);
                        C.push_back(item->second.first);
                        Point2d last = item->second.second;
                        hash.erase(item);
                        for (size_t j = 1;; j++) {
                            item = hash.find(C[j]);
                            assert(item != hash.end());
                            if (!(C[j] == last)) {
                                if (item->second.first == C[j - 1])
                                    C.push_back(item->second.second);
                                else
                                    C.push_back(item->second.first);
                            }
                            hash.erase(item);
                            if (C[j] == last) break;
                        }
                        CS[i].push_back(C);
                    }
                }
            },
            ap
        );
        return CS;
    }
#endif
    //  Finley DR (2007) Point-in-polygon algorithm 
    //  — determining whether a point is inside a complex polygon. 
    //  Available at: http://alienryderflex.com/polygon/.
    // return true = inside, false = outside
    bool is_point_inside_polygon(const Point2d &p, const Polygon &C) {
        bool oddNodes = 0;
        size_t size = C.size();
        for (size_t i = 0; i < size; i++) {
            size_t j = (i == size - 1) ? 0 : i + 1;
            if (((C[i].y < p.y && C[j].y >= p.y) || (C[j].y < p.y && C[i].y >= p.y)) && (C[i].x <= p.x || C[j].x <= p.x)) {
                oddNodes ^= (((p.y - C[i].y) / (C[j].y - C[i].y) * (C[j].x - C[i].x) + C[i].x) < p.x);
            }
        }
        return oddNodes;
    }

    PolygonSides contour_inside_test(const Polygons& C) {
        PolygonSides position(C.size());
        for (size_t i = 0, size = C.size(); i < size; i++) {
            PolygonSide pos = PolygonSide::OUTSIDE;
            for (size_t j = 0, size = C.size(); j < size; j++) {
                if (i == j) continue;
                if (is_point_inside_polygon(C[i][0], C[j])) {
                    pos = PolygonSide::INSIDE;
                    break;
                }
            }
            position[i] = pos;
        }
        return position;
    }

    // This formular came from the cross product of two vectors that is the area of PARALLELOGRAM
    // Then the area of polygon is 1/2 * sum of all parallelogram
    // ref: http://geomalgorithms.com/a01-_area.html
    double measure_polygon_area(const Polygon& C) {
        double A2 = 0;
        for (size_t i = 0, s = C.size(); i < s; i++) {
            size_t i_prev = (i == 0) ? s - 1 : i - 1;
            size_t i_next = (i == s - 1) ? 0 : i + 1;
            A2 += C[i].x * (C[i_next].y - C[i_prev].y);
        }
        return abs(A2) * 0.5;
    }

    double measure_point_square_distance(const Point2d& p1, const Point2d& p2) {
        return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
    }

    // To find orientation of ordered triplet (p, q, r). 
    // The function returns following values 
    // 0 --> p, q and r are colinear 
    // 1 --> Clockwise (Turn Right)
    // -1 --> Counterclockwise (Turn Left)
    int orientation(const Point2d& p, const Point2d& q, const Point2d& r) {
        double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
        if (DOUBLE_EQ(val, 0)) return 0;
        else if (val > 0) return 1;
        else if (val < 0) return -1;
    }

    struct CompareOrientation {
        Point2d origin_;
        CompareOrientation(const Point2d& origin) : origin_(origin) {}
        bool operator() (const Point2d& q, const Point2d& r) {
            int o = orientation(origin_, q, r);
            if (o == 0)
                return (measure_point_square_distance(origin_, r) >= measure_point_square_distance(origin_, q));
            return (o < 0);
        }
    };

    // Graham scan algorithm for constructing convex hull
    // https://www.geeksforgeeks.org/convex-hull-set-2-graham-scan/
    Polygon convexhull(Polygon points) {
        Polygon p;
        if (points.size() < 3) return p;
        double ymin = points[0].y;
        size_t index = 0;
        // Find the bottom-left point
        for (size_t i = 1, size = points.size(); i < size; i++) {
            if (points[i].y < ymin || (DOUBLE_EQ(points[i].y, ymin) && points[i].x < points[index].x)) {
                ymin = points[i].y;
                index = i;
            }
        }
        assert(points.size() > 3 && index >= 0 && index < points.size());
        std::swap(points[0], points[index]);
        Point2d origin = points.front();
        std::sort(points.begin() + 1, points.end(), CompareOrientation(origin));
        // delete the colinear points 
        for (Polygon::iterator p = points.begin() + 1; p != points.end();) {
            Polygon::iterator next_p = (p + 1);
            if (next_p != points.end()) {
                if (orientation(points.front(), *p, *next_p) == 0) {
                    p = points.erase(p);
                }
                else {
                    p++;
                }
            }
            else {
                break;
            }
        }
        if (points.size() < 3) return p;
        // Create convexhull polygon with points[0..2]
        p.reserve(points.size());
        p.push_back(points[0]);
        p.push_back(points[1]);
        p.push_back(points[2]);
        for (size_t i = 3, size = points.size(); i < size; i++) {
            Point2d prev = *(p.end() - 2);
            Point2d current = p.back();
            Point2d next = points[i];
            // check if clockwise (Turn right), then remove current point from convexhull polygon
            while (slice::orientation(prev, current, next) > 0) {
                p.pop_back();
                if (p.size() < 2) break;
                current = p.back();
                prev = *(p.end() - 2);
            }
            p.push_back(next);
        }
        return p;
    }

    FeretDiameter measure_polygon_feret_diameter(const Polygon& p) {
        // find max min in convexhull
        Polygon convex = convexhull(p);
        if (convex.size() < 3) return FeretDiameter();
        size_t size = size = convex.size();
        size_t index_point[4] = { 0,0,0,0 }; // store top-most, right-most, bottom-most, left-most index
        double minX = convex[0].x, minY = convex[0].y, maxX = convex[0].x, maxY = convex[0].y;
        double perimeter = 0;
        for (size_t i = 0; i < size; i++) {
            double x = convex[i].x, y = convex[i].y;
            size_t next_i = (i == size - 1) ? 0 : i + 1;
            perimeter += sqrt(measure_point_square_distance(convex[next_i], convex[i]));
            if (maxY < y || (DOUBLE_EQ(maxY, y) && x < convex[index_point[0]].x)) {
                maxY = y;
                index_point[0] = i;
            }
            if (maxX < x || (DOUBLE_EQ(maxX, x) && y > convex[index_point[1]].y)) {
                maxX = x;
                index_point[1] = i;
            }
            if (minY > y || (DOUBLE_EQ(minY, y) && x > convex[index_point[2]].x)) {
                minY = y;
                index_point[2] = i;
            }
            if (minX > x || (DOUBLE_EQ(minX, x) && y < convex[index_point[3]].y)) {
                minX = x;
                index_point[3] = i;
            }  
        }
        SupportLine2d line[4] = { 
            SupportLine2d(convex[index_point[0]], 0),
            SupportLine2d(convex[index_point[1]], 90),
            SupportLine2d(convex[index_point[2]], 0),
            SupportLine2d(convex[index_point[3]], 90)
        };
        // rotating caliber with perpendicular set of 4 supporting lines
        // the conhexhull's direction is counterclockwise and this procedure direction will be clockwise
        // the loop will finish when the initial upper-most point reached the bottom-most point
        size_t last_index = index_point[1] == 0 ? size - 1 : index_point[1] - 1;
        size_t _count = 0;
        FeretDiameter feret(line);
        feret.perimeter = perimeter;
        while (index_point[0] != last_index && _count < size) {
            // measure the next minimal angle (min_theta) of 4 lines
            double min_angle = 359;
            size_t min_index = 0;
            double angle[4];
            for (int i = 0; i < 4; i++) {
                size_t next_index = index_point[i] == 0 ? size - 1 : (index_point[i] - 1);
                angle[i] = line_point_angle(line[i], convex[next_index]);
                if (min_angle > angle[i]) {
                    min_angle = angle[i];
                    min_index = i;
                }
            }
            //std::cout << " -- min_angle: " << min_angle << " min_index:" << min_index << std::endl;
            //for (int i = 0; i < 4; i++) std::cout << angle[i] << " ";
            //std::cout << std::endl;
            for (int i = 0; i < 4; i++) {
                if (DOUBLE_EQ(min_angle, angle[i])) {
                    size_t next_index = index_point[i] == 0 ? size - 1 : (index_point[i] - 1);
                    line[i].update(convex[next_index]);
                    index_point[i] = next_index;
                }
                assert(DOUBLE_GT_EQ(min_angle, -1) && DOUBLE_LT_EQ(min_angle, 180));
                if (min_angle < 90) line[i] += -min_angle;
                else line[i] += 180 - min_angle;
            }
            //for (int i = 0; i < 4; i++) std::cout << "    " << line[i] << std::endl;
            feret.update(line);
            //std::cout << feret << std::endl;
            _count++;
        }
        return feret;
    }
}

namespace slice {
    // Export the contour to SVG
    bool write_svg(std::string filename, const Polygons &C, const int width, const int height, const int min_x, const int min_y, bool show_convexhull = false) {
        if (C.empty()) return false;
        PolygonSides p = contour_inside_test(C);
        assert(p.size() == C.size());
        std::ofstream svg(filename, std::ofstream::out);
        svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
        svg << "<!DOCTYPE svg PUBLIC \" -//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">" << std::endl;
        svg << "<svg viewBox=\"" << min_x-1 << " " << min_y-1 << " " << width+1 << " " << height+1 << "\" xmlns=\"http://www.w3.org/2000/svg\">" << std::endl;
        for (Polygons::const_iterator t = C.begin(); t != C.end(); t++) {
            if (t->size() < 2) continue;
            size_t index = (size_t)(t - C.begin());
            bool is_inside = (p[index] == PolygonSide::INSIDE);
            // Write origin polygon
            Polygon::const_iterator c = t->begin();
            svg << "<path d=\"M" << c->x << "," << c->y << " L";
            c++;
            for (; c != t->end(); c++) {
                svg << c->x << "," << c->y << " ";
            }
            svg << "z\" fill=\"transparent\" stroke=\"#" << (is_inside ? "f00" : "000") << "\" stroke-width=\"0.1\" stroke-linejoin=\"round\"/>" << std::endl;
            // Write a convex-hull polygon
            if (show_convexhull) {
                Polygon convex = convexhull(*t);
                if (convex.size() >= 3) {
                    c = convex.begin();
                    svg << "<path d=\"M" << c->x << "," << c->y << " L";
                    c++;
                    for (; c != convex.end(); c++) {
                        svg << c->x << "," << c->y << " ";
                    }
                    svg << "z\" stroke-dasharray=\"1,1\" fill=\"transparent\" stroke=\"#" << (is_inside ? "d33" : "666") << "\" stroke-width=\"0.1\" stroke-linejoin=\"round\"/>" << std::endl;
                }
            }
        }
        svg << "</svg>" << std::endl;
        svg.close();
        return true;
    }

    // Measure feret diameter
    void measure_feret_and_shape(const slice::ContourSlice& CS, std::vector<double>& minFeret, std::vector<double>& maxFeret, std::vector<double> (&shapes)[5]) {
        // Foreach slice in 3D
        for (slice::ContourSlice::const_iterator cs = CS.begin(); cs != CS.end(); cs++) {
            if (cs->empty()) continue;
            // Inside or outside testing
            slice::PolygonSides p = contour_inside_test(*cs);
            assert(p.size() == cs->size());
            // For each contour in slice, calculate min max of contour
            std::vector<double> minX, minY, maxX, maxY;
            std::vector<double> sortedMinX, sortedMinY, sortedMaxX, sortedMaxY;
            size_t n_outside = 0, n_inside = 0;
            for (slice::Polygons::const_iterator c = cs->begin(); c != cs->end(); c++) {
                // if polygon must have more-than 2 lines
                if (c->size() <= 2) continue;
                size_t index_polygon = (size_t)(c - cs->begin());
                
                // Loop for each point to find max min in polygon
                slice::Polygon::const_iterator l = c->begin();
                double _minX = l->x, _minY = l->y, _maxX = l->x, _maxY = l->y;
                for (l++; l != c->end(); l++) {
                    if (l->x > _maxX) {
                        _maxX = l->x;
                    }
                    if (_minX > l->x) {
                        _minX = l->x;
                    }
                    if (l->y > _maxY) {
                        _maxY = l->y;
                    }
                    if (_minY > l->y) {
                        _minY = l->y;
                    }
                }
                minX.push_back(_minX);
                maxX.push_back(_maxX);
                minY.push_back(_minY);
                maxY.push_back(_maxY);

                /*
                if (DOUBLE_GT(h, 0) && DOUBLE_GT(w, 0)) {
                    double shape[5] = { area, 4 * area / (M_PI * h * h), 2 * area / (w * h), 4 * area / (M_PI * w * h), P / l };
                }*/
                // count the outside polygon
                if (p[index_polygon] == slice::PolygonSide::OUTSIDE) {
                    // the polygon is the boundary of solid;
                    n_outside++;
                    sortedMinX.push_back(_minX);
                    sortedMaxX.push_back(_maxX);
                    sortedMinY.push_back(_minY);
                    sortedMaxY.push_back(_maxY);
                }
                else {
                    n_inside++;
                    double area = measure_polygon_area(*c);
                    FeretDiameter feret = measure_polygon_feret_diameter(*c);
                    if (!feret.empty && feret.min > 0) {
                        double w = feret.min, h = feret.perpendicularMax;
                        assert(w > 0 && h > 0 && area > 0);
                        // calculate Podczeck's shape description
                        // from: https://diplib.github.io/diplib-docs/features.html#shape_features_PodczeckShapes
                        // the polygon is the boundary of a hole; push the feret diameter to the result
                        minFeret.push_back(feret.min);
                        maxFeret.push_back(feret.perpendicularMax);
                        shapes[0].push_back(area / (w * h));
                        shapes[1].push_back(4 * area / (M_PI * h * h));
                        shapes[2].push_back(2 * area / (w * h));
                        shapes[3].push_back(4 * area / (M_PI * w * h));
                        shapes[4].push_back(feret.perimeter / feret.max);
                    }
                }
            }

            // sorting for searching performances
            std::sort(sortedMinX.begin(), sortedMinX.end());
            std::sort(sortedMaxX.begin(), sortedMaxX.end());
            std::sort(sortedMinY.begin(), sortedMinY.end());
            std::sort(sortedMaxY.begin(), sortedMaxY.end());
            
            // For each contour in slice, calculate feret
            if (n_outside > 1) {
                std::vector<double> feret;
                for (slice::Polygons::const_iterator c = cs->begin(); c != cs->end(); c++) {
                    size_t index = (size_t)(c - cs->begin());
                    if (p[index] == slice::PolygonSide::OUTSIDE) {
                        // Left direction: find previous maxX of left neighbor < minX
                        auto lower = std::lower_bound(sortedMaxX.begin(), sortedMaxX.end(), minX[index]);
                        if (lower != sortedMaxX.end() && lower != sortedMaxX.begin()) feret.push_back(minX[index] - *(lower - 1));
                        // Right direction: find next minX of right neightbor > maxX
                        auto upper = std::upper_bound(sortedMinX.begin(), sortedMinX.end(), maxX[index]);
                        if (upper != sortedMinX.end()) feret.push_back(*upper - maxX[index]);
                        // Upper direction; find previous maxY of upper neighbor < minY
                        lower = std::lower_bound(sortedMaxY.begin(), sortedMaxY.end(), minY[index]);
                        if (lower != sortedMaxY.end() && lower != sortedMaxY.begin()) feret.push_back(minY[index] - *(lower - 1));
                        // lower direction; find next minY of lower neightbor > maxY
                        upper = std::upper_bound(sortedMinY.begin(), sortedMinY.end(), maxY[index]);
                        if (upper != sortedMinY.end()) feret.push_back(*upper - maxY[index]);
                    }
                }
                std::sort(feret.begin(), feret.end());
                if (feret.size() > 0) {
                    minFeret.push_back(feret.front());
                    maxFeret.push_back(feret.back());
                }
            }
        }
    }
}