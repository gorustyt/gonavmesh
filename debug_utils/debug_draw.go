package debug_utils

import (
	"github.com/gorustyt/gonavmesh/common"
	"math"
)

type DuDebugDrawBase struct {
}

func NewDuDebugDraw() DuDebugDraw {
	return &DuDebugDrawBase{}
}
func (d *DuDebugDrawBase) DepthMask(state bool)                                {}
func (d *DuDebugDrawBase) Texture(state bool)                                  {}
func (d *DuDebugDrawBase) Begin(prim DuDebugDrawPrimitives, size ...float32)   {} //TODO float size = 1.0f
func (d *DuDebugDrawBase) Vertex(pos []float32, color Colorb)                  {}
func (d *DuDebugDrawBase) Vertex1(x, y, z float32, color Colorb)               {}
func (d *DuDebugDrawBase) Vertex2(pos []float32, color Colorb, uv []float32)   {}
func (d *DuDebugDrawBase) Vertex3(x, y, z float32, color Colorb, u, v float32) {}
func (d *DuDebugDrawBase) End()                                                {}

func (d *DuDebugDrawBase) AreaToCol(area int) Colorb {
	if area == 0 {
		// Treat zero area type as default.
		return DuRGBA(0, 192, 255, 255)
	} else {
		return DuIntToCol(area, 255)
	}
}

func Bit(a, b int) int {
	return (a & (1 << b)) >> b
}

func DuIntToCol(i, a int) Colorb {
	r := Bit(i, 1) + Bit(i, 3)*2 + 1
	g := Bit(i, 2) + Bit(i, 4)*2 + 1
	b := Bit(i, 0) + Bit(i, 5)*2 + 1
	return DuRGBA(r*63, g*63, b*63, a)
}

func DuIntToCol1(i int, col []float32) {
	r := Bit(i, 0) + Bit(i, 3)*2 + 1
	g := Bit(i, 1) + Bit(i, 4)*2 + 1
	b := Bit(i, 2) + Bit(i, 5)*2 + 1
	col[0] = 1 - float32(r)*63.0/255.0
	col[1] = 1 - float32(g)*63.0/255.0
	col[2] = 1 - float32(b)*63.0/255.0
}

func DuCalcBoxColors(colors []Colorb, colTop Colorb, colSide Colorb) {
	if len(colors) == 0 {
		return
	}

	colors[0] = duMultCol(colTop, 250)
	colors[1] = duMultCol(colSide, 140)
	colors[2] = duMultCol(colSide, 165)
	colors[3] = duMultCol(colSide, 217)
	colors[4] = duMultCol(colSide, 165)
	colors[5] = duMultCol(colSide, 217)
}

func DuDebugDrawCylinderWire(dd DuDebugDraw, minx float32, miny, minz,
	maxx, maxy, maxz float32, col Colorb, lineWidth float32) {
	if dd == nil {
		return
	}

	dd.Begin(DU_DRAW_LINES, lineWidth)
	DuAppendCylinderWire(dd, minx, miny, minz, maxx, maxy, maxz, col)
	dd.End()
}

func DuDebugDrawBoxWire(dd DuDebugDraw, minx, miny, minz,
	maxx, maxy, maxz float32, col Colorb, lineWidth float32) {
	if dd == nil {
		return
	}

	dd.Begin(DU_DRAW_LINES, lineWidth)
	DuAppendBoxWire(dd, minx, miny, minz, maxx, maxy, maxz, col)
	dd.End()
}

func DuDebugDrawArc(dd DuDebugDraw, x0, y0, z0,
	x1, y1, z1, h,
	as0, as1 float32, col Colorb, lineWidth float32) {
	if dd == nil {
		return
	}

	dd.Begin(DU_DRAW_LINES, lineWidth)
	DuAppendArc(dd, x0, y0, z0, x1, y1, z1, h, as0, as1, col)
	dd.End()
}

func DuDebugDrawArrow(dd DuDebugDraw, x0, y0, z0,
	x1, y1, z1,
	as0, as1 float32, col Colorb, lineWidth float32) {
	if dd == nil {
		return
	}

	dd.Begin(DU_DRAW_LINES, lineWidth)
	DuAppendArrow(dd, x0, y0, z0, x1, y1, z1, as0, as1, col)
	dd.End()
}

func DuDebugDrawCircle(dd DuDebugDraw, x, y, z,
	r float32, col Colorb, lineWidth float32) {
	if dd == nil {
		return
	}

	dd.Begin(DU_DRAW_LINES, lineWidth)
	DuAppendCircle(dd, x, y, z, r, col)
	dd.End()
}

func DuDebugDrawCross(dd DuDebugDraw, x, y, z,
	size float32, col Colorb, lineWidth float32) {
	if dd == nil {
		return
	}

	dd.Begin(DU_DRAW_LINES, lineWidth)
	DuAppendCross(dd, x, y, z, size, col)
	dd.End()
}

func DuDebugDrawBox(dd DuDebugDraw, minx, miny, minz,
	maxx, maxy, maxz float32, fcol []Colorb) {
	if dd == nil {
		return
	}

	dd.Begin(DU_DRAW_QUADS)
	DuAppendBox(dd, minx, miny, minz, maxx, maxy, maxz, fcol)
	dd.End()
}

func DuDebugDrawCylinder(dd DuDebugDraw, minx, miny, minz,
	maxx, maxy, maxz float32, col Colorb) {
	if dd == nil {
		return
	}

	dd.Begin(DU_DRAW_TRIS)
	DuAppendCylinder(dd, minx, miny, minz, maxx, maxy, maxz, col)
	dd.End()
}

func DuDebugDrawGridXZ(dd DuDebugDraw, ox, oy, oz float32,
	w, h int, size float32,
	col Colorb, lineWidth float32) {
	if dd == nil {
		return
	}

	dd.Begin(DU_DRAW_LINES, lineWidth)
	for i := 0; i <= h; i++ {
		dd.Vertex1(ox, oy, oz+float32(i)*size, col)
		dd.Vertex1(ox+float32(w)*size, oy, oz+float32(i)*size, col)
	}
	for i := 0; i <= w; i++ {
		dd.Vertex1(ox+float32(i)*size, oy, oz, col)
		dd.Vertex1(ox+float32(i)*size, oy, oz+float32(h)*size, col)
	}
	dd.End()
}

func DuAppendCylinderWire(dd DuDebugDraw, minx, miny, minz,
	maxx, maxy, maxz float32, col Colorb) {
	if dd == nil {
		return
	}

	const NUM_SEG = 16
	dir := make([]float32, NUM_SEG*2)
	init := false
	if !init {
		init = true
		for i := 0; i < NUM_SEG; i++ {
			a := float32(i) / NUM_SEG * math.Pi * 2
			dir[i*2] = float32(math.Cos(float64(a)))
			dir[i*2+1] = float32(math.Sin(float64(a)))
		}
	}

	cx := (maxx + minx) / 2
	cz := (maxz + minz) / 2
	rx := (maxx - minx) / 2
	rz := (maxz - minz) / 2

	i := 0
	j := NUM_SEG - 1
	for i < NUM_SEG {
		dd.Vertex1(cx+dir[j*2+0]*rx, miny, cz+dir[j*2+1]*rz, col)
		dd.Vertex1(cx+dir[i*2+0]*rx, miny, cz+dir[i*2+1]*rz, col)
		dd.Vertex1(cx+dir[j*2+0]*rx, maxy, cz+dir[j*2+1]*rz, col)
		dd.Vertex1(cx+dir[i*2+0]*rx, maxy, cz+dir[i*2+1]*rz, col)
		j = i
		i++
	}
	for i := 0; i < NUM_SEG; i += NUM_SEG / 4 {
		dd.Vertex1(cx+dir[i*2+0]*rx, miny, cz+dir[i*2+1]*rz, col)
		dd.Vertex1(cx+dir[i*2+0]*rx, maxy, cz+dir[i*2+1]*rz, col)
	}
}

func DuAppendBoxWire(dd DuDebugDraw, minx, miny, minz, maxx, maxy, maxz float32, col Colorb) {
	if dd == nil {
		return
	}
	// Top
	dd.Vertex1(minx, miny, minz, col)
	dd.Vertex1(maxx, miny, minz, col)
	dd.Vertex1(maxx, miny, minz, col)
	dd.Vertex1(maxx, miny, maxz, col)
	dd.Vertex1(maxx, miny, maxz, col)
	dd.Vertex1(minx, miny, maxz, col)
	dd.Vertex1(minx, miny, maxz, col)
	dd.Vertex1(minx, miny, minz, col)

	// bottom
	dd.Vertex1(minx, maxy, minz, col)
	dd.Vertex1(maxx, maxy, minz, col)
	dd.Vertex1(maxx, maxy, minz, col)
	dd.Vertex1(maxx, maxy, maxz, col)
	dd.Vertex1(maxx, maxy, maxz, col)
	dd.Vertex1(minx, maxy, maxz, col)
	dd.Vertex1(minx, maxy, maxz, col)
	dd.Vertex1(minx, maxy, minz, col)

	// Sides
	dd.Vertex1(minx, miny, minz, col)
	dd.Vertex1(minx, maxy, minz, col)
	dd.Vertex1(maxx, miny, minz, col)
	dd.Vertex1(maxx, maxy, minz, col)
	dd.Vertex1(maxx, miny, maxz, col)
	dd.Vertex1(maxx, maxy, maxz, col)
	dd.Vertex1(minx, miny, maxz, col)
	dd.Vertex1(minx, maxy, maxz, col)
}

func DuAppendBoxPoints(dd DuDebugDraw, minx, miny, minz,
	maxx, maxy, maxz float32, col Colorb) {
	if dd == nil {
		return
	}
	// Top
	dd.Vertex1(minx, miny, minz, col)
	dd.Vertex1(maxx, miny, minz, col)
	dd.Vertex1(maxx, miny, minz, col)
	dd.Vertex1(maxx, miny, maxz, col)
	dd.Vertex1(maxx, miny, maxz, col)
	dd.Vertex1(minx, miny, maxz, col)
	dd.Vertex1(minx, miny, maxz, col)
	dd.Vertex1(minx, miny, minz, col)

	// bottom
	dd.Vertex1(minx, maxy, minz, col)
	dd.Vertex1(maxx, maxy, minz, col)
	dd.Vertex1(maxx, maxy, minz, col)
	dd.Vertex1(maxx, maxy, maxz, col)
	dd.Vertex1(maxx, maxy, maxz, col)
	dd.Vertex1(minx, maxy, maxz, col)
	dd.Vertex1(minx, maxy, maxz, col)
	dd.Vertex1(minx, maxy, minz, col)
}

func DuAppendBox(dd DuDebugDraw, minx, miny, minz,
	maxx, maxy, maxz float32, fcol []Colorb) {
	if dd == nil {
		return
	}
	verts := [8 * 3]float32{
		minx, miny, minz,
		maxx, miny, minz,
		maxx, miny, maxz,
		minx, miny, maxz,
		minx, maxy, minz,
		maxx, maxy, minz,
		maxx, maxy, maxz,
		minx, maxy, maxz,
	}
	inds := [6 * 4]int{
		7, 6, 5, 4,
		0, 1, 2, 3,
		1, 5, 6, 2,
		3, 7, 4, 0,
		2, 6, 7, 3,
		0, 4, 5, 1,
	}

	in := 0
	for i := 0; i < 6; i++ {
		dd.Vertex(verts[inds[in]*3:], fcol[i])
		in++
		dd.Vertex(verts[inds[in]*3:], fcol[i])
		in++
		dd.Vertex(verts[inds[in]*3:], fcol[i])
		in++
		dd.Vertex(verts[inds[in]*3:], fcol[i])
		in++
	}
}

func DuAppendCylinder(dd DuDebugDraw, minx, miny, minz,
	maxx, maxy, maxz float32, col Colorb) {
	if dd == nil {
		return
	}

	const NUM_SEG = 16
	dir := make([]float32, NUM_SEG*2)
	init := false
	if !init {
		init = true
		for i := 0; i < NUM_SEG; i++ {
			a := float32(i) / NUM_SEG * math.Pi * 2
			dir[i*2] = float32(math.Cos(float64(a)))
			dir[i*2+1] = float32(math.Sin(float64(a)))
		}
	}

	col2 := duMultCol(col, 160)

	cx := (maxx + minx) / 2
	cz := (maxz + minz) / 2
	rx := (maxx - minx) / 2
	rz := (maxz - minz) / 2

	for i := 2; i < NUM_SEG; i++ {
		var a = 0
		var b = i - 1
		var c = i
		dd.Vertex1(cx+dir[a*2+0]*rx, miny, cz+dir[a*2+1]*rz, col2)
		dd.Vertex1(cx+dir[b*2+0]*rx, miny, cz+dir[b*2+1]*rz, col2)
		dd.Vertex1(cx+dir[c*2+0]*rx, miny, cz+dir[c*2+1]*rz, col2)
	}
	for i := 2; i < NUM_SEG; i++ {
		var a = 0
		var b = i
		var c = i - 1
		dd.Vertex1(cx+dir[a*2+0]*rx, maxy, cz+dir[a*2+1]*rz, col)
		dd.Vertex1(cx+dir[b*2+0]*rx, maxy, cz+dir[b*2+1]*rz, col)
		dd.Vertex1(cx+dir[c*2+0]*rx, maxy, cz+dir[c*2+1]*rz, col)
	}
	i := 0
	j := NUM_SEG - 1
	for i < NUM_SEG {
		dd.Vertex1(cx+dir[i*2+0]*rx, miny, cz+dir[i*2+1]*rz, col2)
		dd.Vertex1(cx+dir[j*2+0]*rx, miny, cz+dir[j*2+1]*rz, col2)
		dd.Vertex1(cx+dir[j*2+0]*rx, maxy, cz+dir[j*2+1]*rz, col)

		dd.Vertex1(cx+dir[i*2+0]*rx, miny, cz+dir[i*2+1]*rz, col2)
		dd.Vertex1(cx+dir[j*2+0]*rx, maxy, cz+dir[j*2+1]*rz, col)
		dd.Vertex1(cx+dir[i*2+0]*rx, maxy, cz+dir[i*2+1]*rz, col)
		j = i
		i++
	}
}
func evalArc(x0, y0, z0,
	dx, dy, dz,
	h, u float32, res []float32) {
	res[0] = x0 + dx*u
	res[1] = y0 + dy*u + h*(1-(u*2-1)*(u*2-1))
	res[2] = z0 + dz*u
}

func AppendArrowHead(dd DuDebugDraw, p, q []float32,
	s float32, col Colorb) {
	eps := 0.001
	if dd == nil {
		return
	}
	if common.VdistSqr(p, q) < eps*eps {
		return
	}
	ax := []float32{0, 1, 0}
	ay := []float32{0, 1, 0}
	az := make([]float32, 3)
	common.Vsub(az, q, p)
	common.Vnormalize(az)
	common.Vcross(ax, ay, az)
	common.Vcross(ay, az, ax)
	common.Vnormalize(ay)

	dd.Vertex(p, col)
	//	dd->vertex(p[0]+az[0]*s+ay[0]*s/2, p[1]+az[1]*s+ay[1]*s/2, p[2]+az[2]*s+ay[2]*s/2, col);
	dd.Vertex1(p[0]+az[0]*s+ax[0]*s/3, p[1]+az[1]*s+ax[1]*s/3, p[2]+az[2]*s+ax[2]*s/3, col)

	dd.Vertex(p, col)
	dd.Vertex1(p[0]+az[0]*s-ax[0]*s/3, p[1]+az[1]*s-ax[1]*s/3, p[2]+az[2]*s-ax[2]*s/3, col)

}

func DuAppendArc(dd DuDebugDraw, x0, y0, z0,
	x1, y1, z1, h,
	as0, as1 float32, col Colorb) {
	if dd == nil {
		return
	}
	const NUM_ARC_PTS = 8
	PAD := float32(0.05)
	ARC_PTS_SCALE := (1.0 - PAD*2) / NUM_ARC_PTS
	dx := x1 - x0
	dy := y1 - y0
	dz := z1 - z0
	length := float32(math.Sqrt(float64(dx*dx + dy*dy + dz*dz)))
	prev := make([]float32, 3)
	evalArc(x0, y0, z0, dx, dy, dz, length*h, PAD, prev)
	for i := 1; i <= NUM_ARC_PTS; i++ {
		u := PAD + float32(i)*ARC_PTS_SCALE
		pt := make([]float32, 3)
		evalArc(x0, y0, z0, dx, dy, dz, length*h, u, pt)
		dd.Vertex1(prev[0], prev[1], prev[2], col)
		dd.Vertex1(pt[0], pt[1], pt[2], col)
		prev[0] = pt[0]
		prev[1] = pt[1]
		prev[2] = pt[2]
	}

	// End arrows
	if as0 > 0.001 {
		p := make([]float32, 3)
		q := make([]float32, 3)
		evalArc(x0, y0, z0, dx, dy, dz, length*h, PAD, p)
		evalArc(x0, y0, z0, dx, dy, dz, length*h, PAD+0.05, q)
		AppendArrowHead(dd, p, q, as0, col)
	}

	if as1 > 0.001 {
		p := make([]float32, 3)
		q := make([]float32, 3)
		evalArc(x0, y0, z0, dx, dy, dz, length*h, 1-PAD, p)
		evalArc(x0, y0, z0, dx, dy, dz, length*h, 1-(PAD+0.05), q)
		AppendArrowHead(dd, p, q, as1, col)
	}
}

func DuAppendArrow(dd DuDebugDraw, x0, y0, z0,
	x1, y1, z1,
	as0, as1 float32, col Colorb) {
	if dd == nil {
		return
	}

	dd.Vertex1(x0, y0, z0, col)
	dd.Vertex1(x1, y1, z1, col)

	// End arrows
	p := []float32{x0, y0, z0}
	q := []float32{x1, y1, z1}
	if as0 > 0.001 {
		AppendArrowHead(dd, p, q, as0, col)
	}

	if as1 > 0.001 {
		AppendArrowHead(dd, q, p, as1, col)
	}

}

func DuAppendCircle(dd DuDebugDraw, x, y, z,
	r float32, col Colorb) {
	if dd == nil {
		return
	}
	const NUM_SEG = 40
	dir := make([]float32, 40*2)
	init := false
	if !init {
		init = true
		for i := 0; i < NUM_SEG; i++ {
			a := float32(i) / NUM_SEG * math.Pi * 2
			dir[i*2] = float32(math.Cos(float64(a)))
			dir[i*2+1] = float32(math.Sin(float64(a)))
		}
	}
	i := 0
	j := NUM_SEG - 1
	for i < NUM_SEG {
		dd.Vertex1(x+dir[j*2+0]*r, y, z+dir[j*2+1]*r, col)
		dd.Vertex1(x+dir[i*2+0]*r, y, z+dir[i*2+1]*r, col)
		j = i
		i++
	}
}

func DuAppendCross(dd DuDebugDraw, x, y, z,
	s float32, col Colorb) {
	if dd == nil {
		return
	}
	dd.Vertex1(x-s, y, z, col)
	dd.Vertex1(x+s, y, z, col)
	dd.Vertex1(x, y-s, z, col)
	dd.Vertex1(x, y+s, z, col)
	dd.Vertex1(x, y, z-s, col)
	dd.Vertex1(x, y, z+s, col)
}

func DuRGBA[T int | int32 | int32 | uint8](r, g, b, a T) Colorb {
	return Colorb{uint8(r), uint8(g), uint8(b), uint8(a)}
}

func DuRGBAf(fr, fg, fb, fa float32) Colorb {
	r := int(fr * 255.0)
	g := int(fg * 255.0)
	b := int(fb * 255.0)
	a := int(fa * 255.0)
	return DuRGBA(r, g, b, a)
}

func duMultCol(col Colorb, d uint8) Colorb {
	r := int(col.R() & 0xff)
	g := int((col.G()) & 0xff)
	b := int((col.B()) & 0xff)
	a := int(col.A() & 0xff)
	di := int(d)
	return DuRGBA(r*di>>8, g*di>>8, b*di>>8, a)
}

func DuDarkenCol(col Colorb) (res Colorb) {
	i := col.Int()
	res.FromInt(((i >> 1) & 0x007f7f7f) | (i & 0xff000000))
	return res
}

func DuLerpCol(ca, cb Colorb, u uint8) Colorb {
	ra := ca.R() & 0xff
	ga := (ca.G()) & 0xff
	ba := (ca.B()) & 0xff
	aa := (ca.A()) & 0xff
	rb := cb.R() & 0xff
	gb := (cb.G()) & 0xff
	bb := (cb.B()) & 0xff
	ab := (cb.A()) & 0xff

	r := (ra*(255-u) + rb*u) / 255
	g := (ga*(255-u) + gb*u) / 255
	b := (ba*(255-u) + bb*u) / 255
	a := (aa*(255-u) + ab*u) / 255
	return DuRGBA(r, g, b, a)
}

func DuTransCol(c Colorb, a uint8) Colorb {
	return Colorb{c.R() & 0xff, c.G() & 0xff, c.B() & 0xff, a}

}

type DuDebugDraw interface {
	DepthMask(state bool)

	Texture(state bool)

	/// Begin drawing primitives.
	///  @param prim [in] primitive type to draw, one of rcDebugDrawPrimitives.
	///  @param size [in] size of a primitive, applies to point size and line width only.
	Begin(prim DuDebugDrawPrimitives, size ...float32) //TODO float size = 1.0f

	/// Submit a vertex
	///  @param pos [in] position of the verts.
	///  @param color [in] color of the verts.
	Vertex(pos []float32, color Colorb)

	/// Submit a vertex
	///  @param x,y,z [in] position of the verts.
	///  @param color [in] color of the verts.
	Vertex1(x, y, z float32, color Colorb)

	/// Submit a vertex
	///  @param pos [in] position of the verts.
	///  @param color [in] color of the verts.
	///  @param uv [in] the uv coordinates of the verts.
	Vertex2(pos []float32, color Colorb, uv []float32)

	/// Submit a vertex
	///  @param x,y,z [in] position of the verts.
	///  @param color [in] color of the verts.
	///  @param u,v [in] the uv coordinates of the verts.
	Vertex3(x, y, z float32, color Colorb, u, v float32)

	/// End drawing primitives.
	End()

	/// Compute a color for given area.
	AreaToCol(area int) Colorb
}
type DuDisplayList struct {
	m_pos   []float32
	m_color []Colorb
	m_size  int
	m_cap   int

	m_prim      DuDebugDrawPrimitives
	m_primSize  float32
	m_depthMask bool
}

func NewDuDisplayList(cap int) *DuDisplayList {
	if cap == 0 {
		cap = 512
	}
	d := &DuDisplayList{
		m_primSize:  1.0,
		m_prim:      DU_DRAW_LINES,
		m_depthMask: true,
	}
	if cap < 8 {
		cap = 8
	}

	d.resize(cap)
	return d
}

type DuDebugDrawPrimitives int

const (
	DU_DRAW_POINTS = iota
	DU_DRAW_LINES
	DU_DRAW_TRIS
	DU_DRAW_QUADS
)

func (d *DuDisplayList) resize(cap int) {
	newPos := make([]float32, cap*3)
	if d.m_size != 0 {
		copy(newPos, d.m_pos[:3*d.m_size])
	}
	d.m_pos = newPos

	newColor := make([]Colorb, cap)
	if d.m_size != 0 {
		copy(newColor, d.m_color[:d.m_size])
	}

	d.m_color = newColor

	d.m_cap = cap
}

func (d *DuDisplayList) clear() {
	d.m_size = 0
}

func (d *DuDisplayList) depthMask(state bool) {
	d.m_depthMask = state
}

func (d *DuDisplayList) begin(prim DuDebugDrawPrimitives, size float32) {
	d.clear()
	d.m_prim = prim
	d.m_primSize = size
}

func (d *DuDisplayList) vertex(x, y, z float32, color Colorb) {
	if d.m_size+1 >= d.m_cap {
		d.resize(d.m_cap * 2)
	}

	p := d.m_pos[d.m_size*3:]
	p[0] = x
	p[1] = y
	p[2] = z
	d.m_color[d.m_size] = color
	d.m_size++
}

func (d *DuDisplayList) Vertex(pos []float32, color Colorb) {
	d.vertex(pos[0], pos[1], pos[2], color)
}
func (d *DuDisplayList) end() {
}

func (d *DuDisplayList) draw(dd DuDebugDraw) {
	if dd == nil {
		return
	}
	if d.m_size == 0 {
		return
	}
	dd.DepthMask(d.m_depthMask)
	dd.Begin(d.m_prim, d.m_primSize)
	for i := 0; i < d.m_size; i++ {
		dd.Vertex(d.m_pos[i*3:], d.m_color[i])
	}
	dd.End()
}
