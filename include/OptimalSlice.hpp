/*
 * This is the implement slice algorithm from Rodrigo, 2017 (An Optimal Algorithm for 3D Triangle Mesh Slicing)
 * which is claimed to be faster than slic3r and CGAL method.
 */
#ifndef OPTIMAL_SLICE_HPP
#define OPTIMAL_SLICE_HPP

#ifndef SLICE_PRECISION
// There's a problem of double hashing with the precision less than 1e-8 (e.g. 1e-10)
// when performed contour constructing
#define SLICE_PRECISION 1e-8
#endif
#define DOUBLE_EQ(x,y) (abs(x - y) < SLICE_PRECISION)
#define DOUBLE_GT(x,y) ((x - y) > SLICE_PRECISION)
#define DOUBLE_LT(x,y) ((y - x) > SLICE_PRECISION)
#define USE_PARALLEL

#include "Mesh.h"
#include "tbb/tbb.h"

namespace slice {
    
    enum Direction {X = 0,Y,Z};
    enum ContourPosition {OUTSIDE = 0, INSIDE};

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
                else if DOUBLE_EQ(y, ls.y) return (z < ls.z);
            return false;
        }
    };

    Point2d make_point2d(Point3d p) {
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
                else if (v[1] == ls.v[1]) return (v[2] < ls.v[2]);
            return false;
        }
    };

    class Line {
    public:
        size_t index;
        Point3d v[2];

        Line() {}

        Line(Point3d v0, Point3d v1, size_t index) : index(index)  {
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

    typedef std::vector<double> Plane;
    typedef std::vector<Triangle> Triangles;
    typedef std::vector<Triangles> Layer;
    typedef std::vector<Line> Lines;
    typedef std::vector<Lines> Slice;
    typedef std::vector<Point2d> Contour;
    typedef std::vector<Contour> Contours;
    typedef std::vector<Contours> ContourSlice;
    typedef std::pair<Point2d, Point2d> PairPoint2d;
    typedef std::unordered_map<Point2d, PairPoint2d> ContourHash;
    typedef std::vector<ContourPosition> ContourPositions;

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
        slice::build_triangle_list(mesh, grid_size, P, L);
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
    ContourSlice contour_construct(Slice const& S) {
        ContourSlice CS(S.size());
        ContourHash hash;
        for (size_t i = 0, len = S.size(); i < len; i++) {
            CS[i].clear();
            hash.clear();
            hash.reserve(S[i].size() + 1);
            for (Lines::const_iterator l = S[i].begin(); l != S[i].end(); l++) {
                Point2d u = make_point2d(l->v[0]);
                Point2d v = make_point2d(l->v[1]);
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
            //std::cout << " [Hash OK: " << hash.size() << "]";
            /* TODO: remove this debug message
            std::cout << "Hash Slice " << i << std::endl;
            for (ContourHash::const_iterator item = hash.begin(); item != hash.end(); item++) {
                std::cout << "  " << item->first << ": " << item->second << (item->second.first == item->second.second ? " [***]" : "") << std::endl;
                //std::cout << "  " << item->first << ": " << std::hash<Point2d>()(item->first) << std::endl;
            }
            //*/
            while (!hash.empty()) {
                ContourHash::const_iterator item = hash.begin();
                assert(item != hash.end());
                Contour C;
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
    ContourSlice contour_construct(Slice const& S) {
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
                        Point2d u = make_point2d(l->v[0]);
                        Point2d v = make_point2d(l->v[1]);
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
                        Contour C;
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
    bool is_point_inside_polygon(const Point2d &p, const Contour &C) {
        bool oddNodes = 0;
        size_t size = C.size();
        for (size_t i = 0; i < size; i++) {
            size_t j = (i == size - 1) ? 0 : i + 1;
            if ((C[i].y < p.y && C[j].y >= p.y) || (C[j].y < p.y && C[i].y >= p.y) && (C[i].x <= p.x || C[j].x <= p.x)) {
                oddNodes ^= (((p.y - C[i].y) / (C[j].y - C[i].y) * (C[j].x - C[i].x) + C[i].x) < p.x);
            }
        }
        return oddNodes;
    }

    ContourPositions contour_inside_test(const Contours& C) {
        ContourPositions position(C.size());
        for (size_t i = 0, size = C.size(); i < size; i++) {
            ContourPosition pos = ContourPosition::OUTSIDE;
            for (size_t j = 0, size = C.size(); j < size; j++) {
                if (i == j) continue;
                if (is_point_inside_polygon(C[i][0], C[j])) {
                    pos = ContourPosition::INSIDE;
                    break;
                }
            }
            position[i] = pos;
        }
        return position;
    }
}

namespace slice {
    bool write_svg(std::string filename, const Contours &C, const int width, const int height, const int min_x, const int min_y) {
        if (C.empty()) return false;
        ContourPositions p = contour_inside_test(C);
        assert(p.size() == C.size());
        std::ofstream svg(filename, 'r');
        svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
        svg << "<!DOCTYPE svg PUBLIC \" -//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">" << std::endl;
        svg << "<svg viewBox=\"" << min_x << " " << min_y << " " << width << " " << height << "\" xmlns=\"http://www.w3.org/2000/svg\">" << std::endl;
        for (Contours::const_iterator t = C.begin(); t != C.end(); t++) {
            if (t->size() < 2) continue;
            size_t index = (size_t)(t - C.begin());
            bool is_inside = (p[index] == ContourPosition::INSIDE);
            Contour::const_iterator c = t->begin();
            svg << "<path d=\"M" << c->x << "," << c->y << " L";
            c++;
            for (; c != t->end(); c++) {
                svg << c->x << "," << c->y << " ";
            }
            svg << "z\" fill=\"transparent\" stroke=\"#" << (is_inside ? "f00" : "000") << "\" stroke-width=\"0.1\" stroke-linejoin=\"round\"/>" << std::endl;
        }
        svg << "</svg>" << std::endl;
        svg.close();
        return true;
    }
}
#endif