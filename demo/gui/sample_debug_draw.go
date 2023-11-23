package gui

type duDebugDrawPrimitives int

const (
	DU_DRAW_POINTS duDebugDrawPrimitives = iota
	DU_DRAW_LINES
	DU_DRAW_TRIS
	DU_DRAW_QUADS
)

// / OpenGL debug draw implementation.
// / Abstract debug draw interface.
type DebugDrawGL interface {
	DepthMask(state bool)
	Texture(state bool)
	/// Begin drawing primitives.
	///  @param prim [in] primitive type to draw, one of rcDebugDrawPrimitives.
	///  @param size [in] size of a primitive, applies to point size and line width only.
	Begin(prim *duDebugDrawPrimitives, size float64) //size 默认等于1

	/// Submit a vertex
	///  @param pos [in] position of the verts.
	///  @param color [in] color of the verts.
	Vertex1(pos []float64, color int)

	/// Submit a vertex
	///  @param x,y,z [in] position of the verts.
	///  @param color [in] color of the verts.
	Vertex2(x, y, z float64, color int)

	/// Submit a vertex
	///  @param pos [in] position of the verts.
	///  @param color [in] color of the verts.
	///  @param uv [in] the uv coordinates of the verts.
	Vertex3(pos []float64, color int, uv []float64)

	/// Submit a vertex
	///  @param x,y,z [in] position of the verts.
	///  @param color [in] color of the verts.
	///  @param u,v [in] the uv coordinates of the verts.
	Vertex(x, y, z float64, color int, u, v float64)

	/// End drawing primitives.
	End()

	/// Compute a color for given area.
	AreaToCol(area int) int
}

type SampleTool interface {
	Type() int
	init(sample *Sample)
	reset()
	handleMenu()
	handleClick(s []float64, p []float64, shift bool)
	handleRender()
	handleRenderOverlay(proj, model []float64, view *int)
	handleToggle()
	handleStep()
	handleUpdate(dt float64)
}

type SampleToolState interface {
	init(sample *Sample)
	reset()
	handleRender()
	handleRenderOverlay(proj, model []float64, view *int)
	handleUpdate(dt float64)
}
