/*
 * This is the implement slice algorithm from Rodrigo, 2017 (An Optimal Algorithm for 3D Triangle Mesh Slicing)
 * which is claimed to be faster than slic3r and CGAL method.
 */
#ifndef OPTIMALSLICE_INCLUDED
#define OPTIMALSLICE_INCLUDED

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

#include <iomanip>  
#include "Mesh.h"
#include "MaxHeap.hpp"
#include "oneapi/tbb.h"

namespace optimal_slice {
    
    enum Direction {X = 0,Y,Z};
    enum PolygonSide {OUTSIDE = 0, INSIDE = 1, PAIRED = 2, UNDEFINED = 3};

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
        Point2d operator-(const Point2d &rh) const {
            return {x - rh.x, y - rh.y};
        }
        Point2d operator+(const Point2d& rh) const {
            return { x + rh.x, y + rh.y };
        }
        Point2d operator-() const {
            return { -x, -y };
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
            else if (DOUBLE_EQ(x, (ls.x)))
                if (y < ls.y) { return true; }
                else if (DOUBLE_EQ(y, (ls.y))) {
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
            vcg::Point3f point = face.P(0);
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
            if (v[0] < ls.v[0]) { return true; }
            else if (v[0] == ls.v[0]) {
                if (v[1] < ls.v[1]) { return true; }
                else if (v[1] == ls.v[1]) {
                    return (v[2] < ls.v[2]);
                }
            }
            return false;
        }
    };

    class Line {
    public:
        Point3d v[2];

        Line() {
            v[0] = Point3d{ 0, 0, 0 };
            v[1] = Point3d{ 0, 0, 0 };
        }

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

    // return the 0 < angle <= 90 degree between line and vertex P
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

    // return the 0 <= angle < 180 degree between line and vertex P
    double line_edge_angle(const SupportLine2d& line, const Point2d& delta, bool direction_ccw = false) {
        double angle;
        if (DOUBLE_EQ(delta.x, 0)) angle = 90;
        else if (DOUBLE_EQ(delta.y, 0)) angle = 0;
        else angle = atan(delta.y / delta.x) * 180 / M_PI;
        // calculate delta
        if (DOUBLE_GT_EQ(line.theta, angle)) angle = line.theta - angle;
        else angle = 180 + line.theta - angle;
        if (direction_ccw && angle > 0) angle = 180 - angle;
        // normalized angle
        if (angle < 0) angle += 180;
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
            perimeter = 0;
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
    typedef std::pair<double, size_t> DistanceIndexPair;

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

    std::ostream& operator<< (std::ostream& out, Polygon const& data) {
        out << "slice::Polygon{" << std::endl;
        for (auto d = data.begin(); d != data.end(); ++d) {
            out << " -- " << (*d) << std::endl; 
        }
        out << "}" << std::endl;
        return out;
    }
}

namespace std {
    template<> struct hash<optimal_slice::Point2d> {
        size_t operator()(const optimal_slice::Point2d& p) const noexcept {
            size_t x = hash<long long>()(llround(p.x/SLICE_PRECISION));
            size_t y = hash<long long>()(llround(p.y/SLICE_PRECISION));
            return x ^ (y << 1);
        }
    };
}

namespace optimal_slice {
    
    inline void build_triangle_list(TMesh& mesh, size_t grid_size, Plane& P, Layer& L, int direction = Direction::Z) {
        // Uniform slicing with delta > 0
        // in this case, grid_size = k from Rodrigo paper
        assert(grid_size > 1 && direction <= 2 && direction >= 0);
        vcg::tri::UpdateBounding<TMesh>::Box(mesh);
        vcg::Box3f bbox = mesh.bbox;
        double minBBox = bbox.min[direction];
        double maxBBox = bbox.max[direction];
        vcg::Point3f dim = bbox.Dim();
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
#ifndef USE_PARALLEL
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
#else
        tbb::spin_mutex writeMutex;
        static tbb::affinity_partitioner ap;
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, mesh.face.size()),
            [&](const tbb::blocked_range<size_t> r) {
                // Prepare local_L
                Layer _L(grid_size + 2);
                _L.resize(grid_size + 2);
                for (size_t i = 0; i <= grid_size + 1; i++) _L[i].clear();
                for (size_t i = r.begin(); i < r.end(); i++) {
                    if (!mesh.face[i].IsD()) {
                        Triangle triangle(mesh.face[i]);
                        size_t level = size_t(ceil((triangle.min[direction] - P[1]) / delta) + 1);
                        assert(level > 0 && level <= grid_size + 1);
                        _L[level].push_back(triangle);
                    }
                }
                {
                    tbb::spin_mutex::scoped_lock lock(writeMutex);
                    for (size_t i = 0; i <= grid_size + 1; i++) {
                        L[i].reserve(L[i].size() + _L[i].size());
                        L[i].insert(L[i].end(), _L[i].begin(), _L[i].end());
                    }
                }
            }, ap
        );
#endif
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

    Slice incremental_slicing(TMesh& mesh, size_t k_size, int direction = Direction::Z) {
        optimal_slice::Plane P;
        optimal_slice::Layer L;
        optimal_slice::build_triangle_list(mesh, k_size, P, L, direction);
        Slice S(k_size);
        
        Triangles A;
        for (size_t i = 1; i <= k_size; i++) {
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
                    Lines lines;
                    lines.reserve(local_A.size());
                    std::copy(A.begin() + r.begin(), A.begin() + r.end(), std::back_inserter(local_A));
                    for (Triangles::iterator t = local_A.begin(); t != local_A.end();) {
                        if (t->max[direction] < P[i]) {
                            t = local_A.erase(t);
                        }
                        else {
                            Line line;
                            if (t->max[direction] >= P[i] && t->min[direction] <= P[i]) {
                                if (compute_intersection(*t, P[i], line, direction)) {
                                    //tbb::spin_mutex::scoped_lock lock(writeMutex);
                                    //S[i - 1].push_back(line);
                                    lines.push_back(line);
                                }
                            }
                            t++;
                        }
                    }
                    if (local_A.size() > 0) {
                        tbb::spin_mutex::scoped_lock lock(concatMutex);
                        new_A.insert(new_A.end(), local_A.begin(), local_A.end());
                    }
                    if (lines.size() > 0) {
                        tbb::spin_mutex::scoped_lock lock(writeMutex);
                        S[i - 1].reserve(S[i - 1].size() + lines.size());
                        S[i - 1].insert(S[i - 1].end(), lines.begin(), lines.end());
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
                            // Sometimes, the program got unexpected errors
                            // We prevented it by breaking the process
                            // Perhaps, we have to investigate this incident
                            if (item == hash.end()) {
                                CS[i].clear();
                                break;
                            }
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

    PolygonSides contour_inside_test(const Polygons& C, std::vector<std::pair<size_t, size_t>>& pairList) {
        PolygonSides position(C.size());
        pairList.clear();
        std::vector<size_t> countInside(C.size(), 0);
        std::map<size_t, std::vector<size_t>> pairMapping;
        for (size_t i = 0, size = C.size(); i < size; i++) {
            std::vector<size_t> insidePairedList;
            for (size_t j = 0, size = C.size(); j < size; j++) {
                if (i == j) continue;
                if (is_point_inside_polygon(C[i][0], C[j])) {
                    countInside[i]++;
                    insidePairedList.push_back(j);
                }
            }
            // Fill the outside polygon
            if (countInside[i] == 0) position[i] = PolygonSide::OUTSIDE;
            else if ((countInside[i] % 2) == 0) {
                position[i] = PolygonSide::PAIRED;
                pairMapping.emplace(i, insidePairedList);
            }
        }
        // Match a pair polygon
        for (auto it = pairMapping.begin(); it != pairMapping.end(); ++it) {
            auto target = countInside[it->first] - 1;
            assert(target > 0);
            for (auto list_it = it->second.begin(); list_it != it->second.end(); ++list_it) {
                if (target == countInside[*list_it]) {
                    position[*list_it] = PolygonSide::PAIRED;
                    pairList.push_back(make_pair(it->first, *list_it));
                }
            }
        }
        // Finally, fill the inside polygon
        for (size_t i = 0, len = position.size(); i < len; i++) {
            if ((countInside[i] % 2) == 1 && position[i] == PolygonSide::UNDEFINED)
                position[i] = PolygonSide::INSIDE;
        }
        return position;
    }

    PolygonSides contour_inside_test(const Polygons& C) {
        std::vector<std::pair<size_t, size_t>> pairs;
        return contour_inside_test(C, pairs);
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

    double measure_point_distance(const Point2d& p1, const Point2d& p2) {
        return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
    }

    double measure_point_magnitude(const Point2d& p) {
        return sqrt(p.x * p.x + p.y * p.y);
    }

    // To find orientation of ordered triplet (p, q, r). 
    // The function returns following values 
    // 0 --> p, q and r are colinear 
    // 1 --> Clockwise (Turn Right)
    // -1 --> Counterclockwise (Turn Left)
    int orientation(const Point2d& p, const Point2d& q, const Point2d& r) {
        double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
        if (DOUBLE_EQ(val, 0)) return 0;
        if (val > 0) return 1;
        else if (val < 0) return -1;
        return 0;
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
    // the direction is counterclockwise
    // https://www.geeksforgeeks.org/convex-hull-set-2-graham-scan/
    Polygon convexhull(Polygon points) {
        if (points.size() < 3) return points;
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
        Polygon p;
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
        if (points.size() < 3) return points;
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
            while (optimal_slice::orientation(prev, current, next) > 0) {
                p.pop_back();
                if (p.size() < 2) break;
                current = p.back();
                prev = *(p.end() - 2);
            }
            p.push_back(next);
        }
        return p;
    }

    // Support function generate the farthest point in direction d
    // Complexity O(m+n)
    Point2d gjk_support(const Polygon& P, const Polygon& Q, const Point2d& d) {
        assert(P.size() > 0 && Q.size() > 0);
        Point2d p = P[0];
        double max_dot = p.x * d.x + p.y * d.y;
        for (auto it = P.begin(); it != P.end(); ++it) {
            double dot = it->x * d.x + it->y * d.y;
            if (dot > max_dot) {
                max_dot = dot;
                p = (*it);
            }
        }
        Point2d q = Q[0];
        max_dot = q.x * -d.x + q.y * -d.y;
        for (auto it = Q.begin(); it != Q.end(); ++it) {
            double dot = it->x * -d.x + it->y * -d.y;
            if (dot > max_dot) {
                max_dot = dot;
                q = (*it);
            }
        }
        return p - q;
    }

    // return the origin-closest point on simplex A,B
    // complexity: O(1)
    Point2d gjk_closest_point_to_origin(const Point2d& A, const Point2d& B) {
        Point2d AB = B - A;
        double norm_AB = (AB.x * AB.x + AB.y * AB.y);
        if (DOUBLE_EQ(norm_AB, 0)) return A;
        double lambda_2 = (-AB.x * A.x + -AB.y * A.y) / norm_AB;
        if (lambda_2 > 1) return B;
        else if (lambda_2 < -1) return A;
        double lambda_1 = 1 - lambda_2;
        return { lambda_1 * A.x + lambda_2 * B.x, lambda_1 * A.y + lambda_2 * B.y };
    }

    Point2d center_point_distance(const Polygon& P, const Polygon& Q) {
        Point2d d;
        double x = 0, y = 0;
        for (auto it = P.begin(); it != P.end(); ++it) {
            x += it->x;
            y += it->y;
        }
        d.x = x / P.size();
        d.y = y / P.size();
        x = 0;
        y = 0;
        for (auto it = Q.begin(); it != Q.end(); ++it) {
            x += it->x;
            y += it->y;
        }
        d.x -= x / Q.size();
        d.y -= y / Q.size();
        return d;
    }

    // Minkowski difference of two convex: P and Q
    // return S = P - Q, Complexity O(m+n)
    Polygon minkwoski_difference(const Polygon& P, Polygon Q) {
        Polygon S;
        if (P.size() >= 3 && Q.size() >= 3) {
            // compute -Q and find the bottom-left vertex of Q
            size_t index_q = 0;
            double min_y = -Q[0].y;
            for (Polygon::iterator it = Q.begin(); it != Q.end(); it++) {
                it->x *= -1;
                it->y *= -1;
                if (min_y > it->y || (DOUBLE_EQ(min_y, it->y) && it->x < Q[index_q].x )) {
                    min_y = it->y;
                    index_q = (size_t)(it - Q.begin());
                }
            }
            size_t index_p = 0;
            min_y = P[0].y;
            // find the bottom-left vertex of P
            for (Polygon::const_iterator it = P.begin(); it != P.end(); it++) {
                if (min_y > it->y || (DOUBLE_EQ(min_y, it->y) && it->x < P[index_p].x)) {
                    min_y = it->y;
                    index_p = (size_t)(it - P.begin());
                }
            }
            SupportLine2d line(P[index_p] + Q[index_q], 0);
            S.push_back({ line.x, line.y });
            size_t last_p = index_p;
            size_t _count = 0;
            Point2d first_point = S.front();
            //std::cout << line << std::endl;
            do {
                size_t next_p = (index_p == P.size() - 1) ? 0 : index_p + 1;
                size_t next_q = (index_q == Q.size() - 1) ? 0 : index_q + 1;
                Point2d edge_q = Q[next_q] - Q[index_q];
                Point2d edge_p = P[next_p] - P[index_p];
                double angle_p = line_edge_angle(line, edge_p, true);
                double angle_q = line_edge_angle(line, edge_q, true);
                //std::cout << "Edge P:" << edge_p << " " << angle_p << " Q:" << edge_q << " " << angle_q << std::endl;
                if (angle_p > angle_q) {
                    // insert edge q into S and update index_q
                    line.update(line.x + edge_q.x, line.y + edge_q.y);
                    line += angle_q;
                    index_q = next_q;
                }
                else {
                    line.update(line.x + edge_p.x, line.y + edge_p.y);
                    line += angle_p;
                    index_p = next_p;
                    _count++;
                }
                //std::cout << line << std::endl;
                if (index_p != last_p || _count == 0) {
                    S.push_back({ line.x,line.y });
                }
            } while (index_p != last_p || _count == 0);
        }
        return S;
    }

    // Implement from JAVA lib
    // ref: http://www.dyn4j.org/2010/04/gjk-distance-closest-points/#gjk-closest
    // return -1 if the polygon is intersecting
    double gjk_minimal_distance(const Polygon& P, const Polygon& Q, double tolerance = 1e-2) {
        if (P.size() < 3 || Q.size() < 3) return -1;
        Polygon P1 = convexhull(P);
        Polygon Q1 = convexhull(Q);
        if (P1.size() < 3 || Q1.size() < 3) return -1;
        //Polygon Diff = minkwoski_difference(P1, Q1);
        //assert(Diff.size() > 3);
        //std::copy(Diff.begin(), Diff.end(), std::ostream_iterator<slice::Point2d>(std::cout, " "));
        Point2d d = center_point_distance(P1, Q1);
        Point2d simplex_a = gjk_support(P1, Q1, d);
        Point2d simplex_b = gjk_support(P1, Q1, -d);
        d = -gjk_closest_point_to_origin(simplex_a, simplex_b);
        //std::cout << std::endl << "Pre: a:" << simplex_a << " b:" << simplex_b << " d:" << d << std::endl;
        size_t _count = P1.size() + Q1.size();
        while (_count > 0) {
            if (DOUBLE_EQ(d.x * d.x + d.y * d.y, 0)) return 0;
            Point2d c = gjk_support(P1, Q1, d);
            //std::cout << "  -- c:" << c << std::endl;
            // check new point c is better than simplex a, b
            if (( d.x * c.x + d.y * c.y ) - (simplex_a.x * d.x + simplex_a.y * d.y) < tolerance) {
                return measure_point_magnitude(d);
            } 
            Point2d p1 = gjk_closest_point_to_origin(simplex_a, c);
            Point2d p2 = gjk_closest_point_to_origin(c, simplex_b);
            if (measure_point_magnitude(p1) < measure_point_magnitude(p2)) {
                simplex_b = c;
                d = -p1;
            }
            else {
                simplex_a = c;
                d = -p2;
            }
            //std::cout << "  -- a:" << simplex_a << " b:" << simplex_b << " d:" << d << std::endl;
            _count--;
        }
        //[Warning] GJK Iteration timeout: Found the intersect polygon
        return measure_point_magnitude(d);
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

namespace optimal_slice {
    // Export the contour to SVG
    bool write_svg(std::string filename, const Polygons &C, const float width, const float height, const int min_x, const int min_y, bool show_convexhull = false) {
        if (C.empty()) return false;
        PolygonSides p = contour_inside_test(C);
        assert(p.size() == C.size());
        std::ofstream svg(filename, std::ofstream::out);
        float strokeWidth = (width + height) / 400.0f;
        svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
        svg << "<!DOCTYPE svg PUBLIC \" -//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">" << std::endl;
        svg << "<svg viewBox=\"" << min_x-strokeWidth << " " << min_y-strokeWidth << " " << width+4*strokeWidth << " " << height+4*strokeWidth << "\" xmlns=\"http://www.w3.org/2000/svg\">" << std::endl;
       
        for (Polygons::const_iterator t = C.begin(); t != C.end(); t++) {
            if (t->size() < 2) continue;
            size_t index = (size_t)(t - C.begin());
            std:string color = "eee";
            if (p[index] == PolygonSide::INSIDE) color = "f00";
            else if (p[index] == PolygonSide::OUTSIDE) color = "000";
            else if (p[index] == PolygonSide::PAIRED) color = "00f";
            // Write origin polygon
            Polygon::const_iterator c = t->begin();
            svg << "<path d=\"M" << c->x << "," << c->y << " L";
            c++;
            for (; c != t->end(); c++) {
                svg << c->x << "," << c->y << " ";
            }
            svg << "z\" fill=\"transparent\" stroke=\"#" << color << "\" stroke-width=\""<< std::fixed << std::setprecision(4) << strokeWidth << "\" stroke-linejoin=\"round\"/>" << std::endl;
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
                    svg << "z\" stroke-dasharray=\"1,1\" fill=\"transparent\" stroke=\"#d33\" stroke-width=\"" << std::fixed << std::setprecision(4) << strokeWidth << "\" stroke-linejoin=\"round\"/>" << std::endl;
                }
            }
        }
        svg << "</svg>" << std::endl;
        svg.close();
        return true;
    }

    // Measure feret diameter
    void measure_feret_and_shape(const optimal_slice::ContourSlice& CS, size_t k_polygon, std::vector<double>& minFeret, std::vector<double>& maxFeret, std::vector<double> (&shapes)[5]) {
        // Foreach slice in 3D
#ifndef USE_PARALLEL
        for (size_t cs_index = 0, cs_size = CS.size(); cs_index < cs_size; cs_index++) {
#else
        static tbb::affinity_partitioner ap;
        tbb::spin_mutex feretMutex, shapeMutex;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, CS.size()), [&](tbb::blocked_range<size_t> r) {
            for (size_t cs_index = r.begin(); cs_index < r.end(); cs_index++) {
#endif
                if (CS[cs_index].empty()) continue;
                // Inside or outside testing
                std::vector<std::pair<size_t, size_t>> pairPolygons;
                optimal_slice::PolygonSides p = contour_inside_test(CS[cs_index], pairPolygons);

                assert(p.size() == CS[cs_index].size());
                // For each contour in slice, calculate min max of contour
                std::vector< Point2d > polygonCenters;
                std::vector< std::pair<Point2d, size_t> > centerPolygonMap;
                size_t n_outside = 0, n_inside = 0;
                for (optimal_slice::Polygons::const_iterator c = CS[cs_index].begin(); c != CS[cs_index].end(); c++) {
                    // the polygon must have more-than 2 lines
                    if (c->size() <= 2) {
                        polygonCenters.push_back({ 0, 0 });
                        continue;
                    }
                    size_t index_polygon = (size_t)(c - CS[cs_index].begin());

                    // Loop for each point to find max min in polygon
                    optimal_slice::Polygon::const_iterator l = c->begin();
                    double _x = 0, _y = 0;
                    for (l++; l != c->end(); l++) {
                        _x += l->x;
                        _y += l->y;
                    }
                    _x /= c->size();
                    _y /= c->size();
                    polygonCenters.push_back({ _x, _y });

                    if (p[index_polygon] == optimal_slice::PolygonSide::OUTSIDE) {
                        // the polygon is the boundary of solid;
                        // count the outside polygon
                        n_outside++;
                        centerPolygonMap.push_back(std::make_pair(Point2d{ _x, _y }, index_polygon));
                    }
                    // 1. Rotating caliber for INSIDE polygon
                    else if (p[index_polygon] == optimal_slice::PolygonSide::INSIDE) {
                        n_inside++;
                        double area = measure_polygon_area(*c);
                        FeretDiameter feret = measure_polygon_feret_diameter(*c);
                        if (!feret.empty && feret.min > 0) {
                            double w = feret.min, h = feret.perpendicularMax;
                            assert(w > 0 && h > 0 && area > 0);
                            // calculate Podczeck's shape description
                            // from: https://diplib.github.io/diplib-docs/features.html#shape_features_PodczeckShapes
                            // the polygon is the boundary of a hole; push the feret diameter to the result
                            {
#ifdef USE_PARALLEL
                                tbb::spin_mutex::scoped_lock lock(feretMutex);
#endif
                                minFeret.push_back(feret.min);
                                maxFeret.push_back(feret.perpendicularMax);
                            }
                            {
#ifdef USE_PARALLEL
                                tbb::spin_mutex::scoped_lock lock(shapeMutex);
#endif
                                shapes[0].push_back(area / (w * h));
                                shapes[1].push_back(4 * area / (M_PI * h * h));
                                shapes[2].push_back(2 * area / (w * h));
                                shapes[3].push_back(4 * area / (M_PI * w * h));
                                shapes[4].push_back(feret.perimeter / feret.max);
                            }
                        }
                    }
                }

                // 2. GJK distance for OUTSIDE polygon
                if (n_outside > 1 && centerPolygonMap.size() > 1) {
                    // For each Polygon in slice layer, find 4 minimum-distance polygon
                    for (optimal_slice::Polygons::const_iterator c = CS[cs_index].begin(); c != CS[cs_index].end(); c++) {
                        if (c->size() <= 2) continue;
                        size_t index = (size_t)(c - CS[cs_index].begin());
                        std::vector<double> feret;
                        double _feret;
                        std::vector<DistanceIndexPair> pairs;
                        int k = std::min(centerPolygonMap.size() - 1, k_polygon);
                        if (p[index] == optimal_slice::PolygonSide::OUTSIDE) {
                            size_t i = 0;
                            // Take the first k in centerPolygonMap for initializing the MaxHeap
                            for (size_t j = 0; j < k; i++) {
                                if (i != index && i < centerPolygonMap.size()) {
                                    double distance = measure_point_distance(polygonCenters[index], centerPolygonMap[i].first);
                                    pairs.push_back(std::make_pair(distance, centerPolygonMap[i].second));
                                    j++;
                                }
                            }
                            MaxHeap<DistanceIndexPair> minimum(k, &pairs);
                            // Process the other centerPolygonMap
                            for (; i < centerPolygonMap.size(); i++) {
                                double distance = measure_point_distance(polygonCenters[index], centerPolygonMap[i].first);
                                minimum.update(std::make_pair(distance, centerPolygonMap[i].second));
                            }
                            // Use the k closest centerPolygonMap for calculating the pore size
                            for (i = 0; i < k; i++) {
                                size_t next_index = pairs[i].second;
                                if (next_index < CS[cs_index].size()) {
                                    _feret = gjk_minimal_distance(*c, CS[cs_index].at(next_index));
                                    if (_feret > 0) feret.push_back(_feret);
                                }
                            }
                        }
                        if (feret.size() > 0) {
                            std::sort(feret.begin(), feret.end());
#ifdef USE_PARALLEL
                            tbb::spin_mutex::scoped_lock lock(feretMutex);
#endif
                            minFeret.push_back(feret.front());
                            maxFeret.push_back(feret.back());
                        }
                    }
                }

                // 3. GJK distance for PAIRED polygon
                for (auto it = pairPolygons.begin(); it != pairPolygons.end(); ++it) {
                    double feret = gjk_minimal_distance(CS[cs_index].at(it->first), CS[cs_index].at(it->second));
#ifdef USE_PARALLEL
                    tbb::spin_mutex::scoped_lock lock(feretMutex);
#endif
                    minFeret.push_back(feret);
                    maxFeret.push_back(feret);
                }
            }
#ifdef USE_PARALLEL
            }, ap);
#endif
    }
}
#endif