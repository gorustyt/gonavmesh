package recast

type RcPolyMesh struct {
	Npolys       int       ///< The number of polygons.
	Verts        []int     ///< The mesh vertices.
	Polys        []int     ///< Polygon and neighbor data.
	Regs         []int     ///< The region id assigned to each polygon. [Length: #maxpolys]
	Flags        []int     ///< The user defined flags for each polygon. [Length: #maxpolys]
	Areas        []int     ///< The area id assigned to each polygon. [Length: #maxpolys]
	Bmin         []float64 ///< The minimum bounds in world space. [(x, y, z)]
	Bmax         []float64 ///< The maximum bounds in world space. [(x, y, z)]
	Cs           float64   ///< The size of each cell. (On the xz-plane.)
	Ch           float64   ///< The height of each cell. (The minimum increment along the y-axis.)
	BorderSize   int       ///< The AABB border size used to generate the source data from which the mesh was derived.
	MaxEdgeError float64   ///< The max error of the polygon edges in the mesh.
	Nvp          int       //多边形最大能生成的顶点数
	Nverts       int       ///< The number of vertices.
	Maxpolys     int       ///< The number of allocated polygons.
}

type RcContourSet struct {
	conts      []*rcContour ///< An array of the contours in the set. [Size: #nconts]
	Nconts     int          ///< The number of contours in the set.
	bmin       []float64    ///< The minimum bounds in world space. [(x, y, z)]
	bmax       []float64    ///< The maximum bounds in world space. [(x, y, z)]
	cs         float64      ///< The size of each cell. (On the xz-plane.)
	ch         float64      ///< The height of each cell. (The minimum increment along the y-axis.)
	width      int          ///< The width of the set. (Along the x-axis in cell units.)
	height     int          ///< The height of the set. (Along the z-axis in cell units.)
	borderSize int          ///< The AABB border size used to generate the source data from which the mesh was derived.
	maxError   float64      ///< The max edge error that this contour set was simplified with.
}

// 轮廓
type rcContour struct {
	verts   []int ///< The mesh vertices.
	nverts  int   ///< The number of vertices in the simplified contour.
	rverts  []int ///< Raw contour vertex and connection data. [Size: 4 * #nrverts]
	nrverts int   ///< The number of vertices in the raw contour.
	reg     int   ///< The region id of the contour.
	area    int   ///< The area id of the contour.
}
