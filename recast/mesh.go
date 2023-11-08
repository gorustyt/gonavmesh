package recast

// 多边形
type polygon struct {
	areaId string       //多边形区域id
	verts  []*Point     //顶点
	flags  map[int]bool //标志
}

type rcPolyMesh struct {
	verts        []*Point   ///< The mesh vertices.
	polys        []*polygon ///< Polygon and neighbor data.
	bmin         Point      ///< The minimum bounds in world space. [(x, y, z)]
	bmax         Point      ///< The maximum bounds in world space. [(x, y, z)]
	cs           float64    ///< The size of each cell. (On the xz-plane.)
	ch           float64    ///< The height of each cell. (The minimum increment along the y-axis.)
	borderSize   int        ///< The AABB border size used to generate the source data from which the mesh was derived.
	maxEdgeError float64    ///< The max error of the polygon edges in the mesh.
	nvp          int        //多边形最大能生成的顶点数
}

type rcContourSet struct {
	conts      []*rcContour ///< An array of the contours in the set. [Size: #nconts]
	bmin       Point        ///< The minimum bounds in world space. [(x, y, z)]
	bmax       Point        ///< The maximum bounds in world space. [(x, y, z)]
	cs         float64      ///< The size of each cell. (On the xz-plane.)
	ch         float64      ///< The height of each cell. (The minimum increment along the y-axis.)
	width      int          ///< The width of the set. (Along the x-axis in cell units.)
	height     int          ///< The height of the set. (Along the z-axis in cell units.)
	borderSize int          ///< The AABB border size used to generate the source data from which the mesh was derived.
	maxError   float64      ///< The max edge error that this contour set was simplified with.
}

// 轮廓
type rcContour struct {
	verts  []*Point ///< The mesh vertices.
	rverts []*Point ///< Raw contour vertex and connection data. [Size: 4 * #nrverts]
	reg    int      ///< The region id of the contour.
	areaId string   ///< The area id of the contour.
}
