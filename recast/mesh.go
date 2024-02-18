package recast

type RcPolyMesh struct {
	Npolys       int32     ///< The number of polygons.
	Verts        []uint16  ///< The mesh vertices.
	Polys        []uint16  ///< Polygon and neighbor data.
	Regs         []uint16  ///< The region id assigned to each polygon. [Length: #maxpolys]
	Flags        []uint16  ///< The user defined flags for each polygon. [Length: #maxpolys]
	Areas        []uint8   ///< The area id assigned to each polygon. [Length: #maxpolys]
	Bmin         []float32 ///< The minimum bounds in world space. [(x, y, z)]
	Bmax         []float32 ///< The maximum bounds in world space. [(x, y, z)]
	Cs           float32   ///< The size of each cell. (On the xz-plane.)
	Ch           float32   ///< The height of each cell. (The minimum increment along the y-axis.)
	BorderSize   int32     ///< The AABB border size used to generate the source data from which the mesh was derived.
	MaxEdgeError float32   ///< The max error of the polygon edges in the mesh.
	Nvp          int32     //多边形最大能生成的顶点数
	Nverts       int32     ///< The number of vertices.
	Maxpolys     int32     ///< The number of allocated polygons.
}

type RcContourSet struct {
	Conts      []*RcContour ///< An array of the contours in the set. [Size: #nconts]
	Nconts     int32        ///< The number of contours in the set.
	Bmin       []float32    ///< The minimum bounds in world space. [(x, y, z)]
	Bmax       []float32    ///< The maximum bounds in world space. [(x, y, z)]
	Cs         float32      ///< The size of each cell. (On the xz-plane.)
	Ch         float32      ///< The height of each cell. (The minimum increment along the y-axis.)
	Width      int32        ///< The width of the set. (Along the x-axis in cell units.)
	Height     int32        ///< The height of the set. (Along the z-axis in cell units.)
	BorderSize int32        ///< The AABB border size used to generate the source data from which the mesh was derived.
	MaxError   float32      ///< The max edge error that this contour set was simplified with.
}

// 轮廓
type RcContour struct {
	Verts   []int32 ///< The mesh vertices.
	Nverts  int32   ///< The number of vertices in the simplified contour.
	Rverts  []int32 ///< Raw contour vertex and connection data. [Size: 4 * #nrverts]
	Nrverts int     ///< The number of vertices in the raw contour.
	Reg     uint16  ///< The region id of the contour.
	Area    uint8   ///< The area id of the contour.
}
