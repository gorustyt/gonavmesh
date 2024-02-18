package recast

import (
	"fmt"
	"strconv"
	"testing"
)

func assertTrue(t *testing.T, value bool, msg string) {
	if !value {
		t.Errorf(msg)
	}
}

func TestRcClamp(t *testing.T) {
	assertTrue(t, rcClamp(2, 0, 1) == 1, "Higher than range error")
	assertTrue(t, rcClamp(1, 0, 2) == 1, "Within range error")
	assertTrue(t, rcClamp(0, 1, 2) == 1, "Lower than range error")
}

func TestSqr(t *testing.T) {
	if rcSqr(2) != 4 {
		t.Errorf("Sqr squares a number")
	}
	if rcSqr(-4) != 16 {
		t.Errorf("Sqr squares a number")
	}
	if rcSqr(0) != 0 {
		t.Errorf("Sqr squares a number")
	}
}

func TestVcross(t *testing.T) {
	v1 := []float32{3, -3, 1}
	v2 := []float32{4, 9, 2}
	result := make([]float32, 3)
	rcVcross(result, v1, v2)
	if result[0] != -15 {
		t.Errorf("Computes cross product")
	}
	if result[1] != -2 {
		t.Errorf("Computes cross product")
	}
	if result[2] != 39 {
		t.Errorf("Computes cross product")
	}
	result = make([]float32, 3)
	v1 = []float32{3, -3, 1}
	rcVcross(result, v1, v1)
	if result[0] != 0 {
		t.Errorf("Cross product with itself is zero")
	}
	if result[1] != 0 {
		t.Errorf("Cross product with itself is zero")
	}
	if result[2] != 0 {
		t.Errorf("Cross product with itself is zero")
	}
}

func TestVdot(t *testing.T) {
	v1 := []float32{1, 0, 0}
	result := rcVdot(v1, v1)
	if result != 1 {
		t.Errorf("Dot normalized vector with itself")
	}

	v1 = []float32{1, 2, 3}
	v2 := []float32{0, 0, 0}
	result = rcVdot(v1, v2)
	if result != 0 {
		t.Errorf("Dot zero vector with anything is zero")
	}
}
func TestVdist(t *testing.T) {
	v1 := []float32{3, 1, 3}
	v2 := []float32{1, 3, 1}
	result := rcVdist(v1, v2)
	value, _ := strconv.ParseFloat(fmt.Sprintf("%.4f", result), 64)
	if value != 3.4641 {
		t.Errorf("distance between two vectors")
	}

	v1 = []float32{3, 1, 3}
	v2 = []float32{0, 0, 0}
	result = rcVdist(v1, v2)
	magnitude := rcSqrt(rcSqr(v1[0]) + rcSqr(v1[1]) + rcSqr(v1[2]))
	if result != magnitude {
		t.Errorf("Distance from zero is magnitude")
	}
}
func TestVdistSqr(t *testing.T) {
	v1 := []float32{3, 1, 3}
	v2 := []float32{1, 3, 1}
	result := rcVdistSqr(v1, v2)
	if result != 12 {
		t.Errorf("squared distance between two vectors")
	}

	v1 = []float32{3, 1, 3}
	v2 = []float32{0, 0, 0}
	result = rcVdistSqr(v1, v2)
	magnitude := rcSqr(v1[0]) + rcSqr(v1[1]) + rcSqr(v1[2])
	if result != magnitude {
		t.Errorf("squared distance from zero is squared magnitude")
	}

}

func TestVnormalize(t *testing.T) {
	v := []float32{3, 3, 3}
	rcVnormalize(v)
	if v[0] != rcSqrt(1.0/3.0) {
		t.Errorf("normalizing reduces magnitude to 1")
	}
	if v[1] != rcSqrt(1.0/3.0) {
		t.Errorf("normalizing reduces magnitude to 1")
	}
	if v[2] != rcSqrt(1.0/3.0) {
		t.Errorf("normalizing reduces magnitude to 1")
	}
	magnitude := rcSqrt(rcSqr(v[0]) + rcSqr(v[1]) + rcSqr(v[2]))
	if magnitude != 1 {
		t.Errorf("normalizing reduces magnitude to 1")
	}
}

func TestCalcBounds(t *testing.T) {
	verts := []float32{1, 2, 3}
	bmin := make([]float32, 3)
	bmax := make([]float32, 3)
	rcCalcBounds(verts, 1, bmin, bmax)

	if bmin[0] != verts[0] {
		t.Errorf("bounds of one vector")
	}
	if bmin[1] != verts[1] {
		t.Errorf("bounds of one vector")
	}
	if bmin[2] != verts[2] {
		t.Errorf("bounds of one vector")
	}

	if bmax[0] != verts[0] {
		t.Errorf("bounds of one vector")
	}
	if bmax[1] != verts[1] {
		t.Errorf("bounds of one vector")
	}
	if bmax[2] != verts[2] {
		t.Errorf("bounds of one vector")
	}

	verts = []float32{1, 2, 3,
		0, 2, 5}
	bmin = make([]float32, 3)
	bmax = make([]float32, 3)
	rcCalcBounds(verts, 2, bmin, bmax)

	if bmin[0] != 0 {
		t.Errorf("bounds of one vector")
	}
	if bmin[1] != 2 {
		t.Errorf("bounds of one vector")
	}
	if bmin[2] != 3 {
		t.Errorf("bounds of one vector")
	}

	if bmax[0] != 1 {
		t.Errorf("bounds of one vector")
	}
	if bmax[1] != 2 {
		t.Errorf("bounds of one vector")
	}
	if bmax[2] != 5 {
		t.Errorf("bounds of one vector")
	}
}
func TestCalcGridSize(t *testing.T) {
	verts := []float32{1, 2, 3,
		0, 2, 6}
	bmin := make([]float32, 3)
	bmax := make([]float32, 3)
	rcCalcBounds(verts, 2, bmin, bmax)
	cellSize := 1.5

	var width int
	var height int

	rcCalcGridSize(bmin, bmax, cellSize, &width, &height)

	if width != 1 {
		t.Errorf("computes the size of an x & z axis grid")
	}
	if height != 2 {
		t.Errorf("computes the size of an x & z axis grid")
	}
}

func TestCreateHeightfield(t *testing.T) {
	verts := []float32{1, 2, 3,
		0, 2, 6}
	bmin := make([]float32, 3)
	bmax := make([]float32, 3)
	rcCalcBounds(verts, 2, bmin, bmax)
	cellSize := 1.5
	cellHeight := 2.0
	var width int
	var height int

	rcCalcGridSize(bmin, bmax, cellSize, &width, &height)
	heightfield := &RcHeightfield{}

	result := rcCreateHeightfield(heightfield, width, height, bmin, bmax, cellSize, cellHeight)
	msg := "create a heightfield"
	assertTrue(t, result, msg)

	assertTrue(t, heightfield.width == width, msg)
	assertTrue(t, heightfield.height == height, msg)

	assertTrue(t, heightfield.bmin[0] == bmin[0], msg)
	assertTrue(t, heightfield.bmin[1] == bmin[1], msg)
	assertTrue(t, heightfield.bmin[2] == bmin[2], msg)

	assertTrue(t, heightfield.bmax[0] == bmax[0], msg)
	assertTrue(t, heightfield.bmax[1] == bmax[1], msg)
	assertTrue(t, heightfield.bmax[2] == bmax[2], msg)

	assertTrue(t, heightfield.cs == cellSize, msg)
	assertTrue(t, heightfield.ch == cellHeight, msg)

	assertTrue(t, heightfield.spans != nil, msg)
	assertTrue(t, heightfield.pools == nil, msg)
	assertTrue(t, heightfield.freelist == nil, msg)
}

func TestMarkWalkableTriangles(t *testing.T) {
	walkableSlopeAngle := 45.0
	nt := 1
	nv := 3
	var (
		verts          []float32
		walkable_tri   []int
		unwalkable_tri []int
		areas          []int
	)
	reset := func() {
		verts = []float32{
			0, 0, 0,
			1, 0, 0,
			0, 0, -1,
		}
		walkable_tri = []int{0, 1, 2}
		unwalkable_tri = []int{0, 2, 1}
		areas = []int{RC_NULL_AREA}
	}
	reset()
	rcMarkWalkableTriangles(walkableSlopeAngle, verts, nv, walkable_tri, nt, areas)
	if areas[0] != RC_WALKABLE_AREA {
		t.Errorf("One walkable triangle")
	}
	reset()
	rcMarkWalkableTriangles(walkableSlopeAngle, verts, nv, unwalkable_tri, nt, areas)
	if areas[0] != RC_NULL_AREA {
		t.Errorf("One non-walkable triangle")
	}
	reset()
	areas[0] = 42
	rcMarkWalkableTriangles(walkableSlopeAngle, verts, nv, unwalkable_tri, nt, areas)
	if areas[0] != 42 {
		t.Errorf("Non-walkable triangle area id's are not modified")
	}
	reset()
	walkableSlopeAngle = 0
	rcMarkWalkableTriangles(walkableSlopeAngle, verts, nv, walkable_tri, nt, areas)
	if areas[0] != RC_NULL_AREA {
		t.Errorf("Slopes equal to the max slope are considered unwalkable.")
	}
}

func TestClearUnwalkableTriangles(t *testing.T) {
	walkableSlopeAngle := 45.0
	nt := 1
	nv := 3
	var (
		verts          []float32
		walkable_tri   []int
		unwalkable_tri []int
		areas          []int
	)
	reset := func() {
		verts = []float32{
			0, 0, 0,
			1, 0, 0,
			0, 0, -1,
		}
		walkable_tri = []int{0, 1, 2}
		unwalkable_tri = []int{0, 2, 1}
		areas = []int{42}
	}
	reset()
	rcClearUnwalkableTriangles(walkableSlopeAngle, verts, nv, unwalkable_tri, nt, areas)
	if areas[0] != RC_NULL_AREA {
		t.Errorf("Sets area ID of unwalkable triangle to RC_NULL_AREA")
	}
	reset()
	rcClearUnwalkableTriangles(walkableSlopeAngle, verts, nv, walkable_tri, nt, areas)
	if areas[0] != 42 {
		t.Errorf("Does not modify walkable triangle aread ID's")
	}
	reset()
	walkableSlopeAngle = 0
	rcClearUnwalkableTriangles(walkableSlopeAngle, verts, nv, walkable_tri, nt, areas)
	if areas[0] != RC_NULL_AREA {
		t.Errorf("Does not modify walkable triangle aread ID's")
	}
}

func TestAddSpan(t *testing.T) {
	var (
		verts        []float32
		bmin         []float32
		bmax         []float32
		cellSize     float32
		cellHeight   float32
		width        int
		height       int
		hf           RcHeightfield
		x            int
		y            int
		smin         int
		smax         int
		area         int
		flagMergeThr int
	)
	reset := func() {
		verts = []float32{
			1, 2, 3,
			0, 2, 6,
		}
		bmin = make([]float32, 3)
		bmax = make([]float32, 3)
		rcCalcBounds(verts, 2, bmin, bmax)

		cellSize = 1.5
		cellHeight = 2.0
		width = 0
		height = 0

		rcCalcGridSize(bmin, bmax, cellSize, &width, &height)

		hf = RcHeightfield{}
		if !rcCreateHeightfield(&hf, width, height, bmin, bmax, cellSize, cellHeight) {
			t.Errorf("rcAddSpan")
		}

		x = 0
		y = 0
		smin = 0
		smax = 1
		area = 42
		flagMergeThr = 1

	}
	reset()
	result := rcAddSpan(&hf, x, y, smin, smax, area, flagMergeThr)
	if !result {
		t.Errorf("Add a span to an empty heightfield.")
	}
	if hf.spans[0] == nil {
		t.Errorf("Add a span to an empty heightfield.")
	}
	if hf.spans[0].smin != smin {
		t.Errorf("Add a span to an empty heightfield.")
	}
	if hf.spans[0].smax != smax {
		t.Errorf("Add a span to an empty heightfield.")
	}
	if hf.spans[0].area != area {
		t.Errorf("Add a span to an empty heightfield.")
	}
	reset()
	result = rcAddSpan(&hf, x, y, smin, smax, area, flagMergeThr)
	if !result {
		t.Errorf("Add a span that gets merged with an existing span.")
	}
	if hf.spans[0] == nil {
		t.Errorf("Add a span that gets merged with an existing span.")
	}
	if hf.spans[0].smin != smin {
		t.Errorf("Add a span that gets merged with an existing span.")
	}
	if hf.spans[0].smax != smax {
		t.Errorf("Add a span that gets merged with an existing span.")
	}
	if hf.spans[0].area != area {
		t.Errorf("Add a span that gets merged with an existing span.")
	}

	smin = 1
	smax = 2
	result = rcAddSpan(&hf, x, y, smin, smax, area, flagMergeThr)
	if !result {
		t.Errorf("Add a span that gets merged with an existing span.")
	}
	if hf.spans[0] == nil {
		t.Errorf("Add a span that gets merged with an existing span.")
	}
	if hf.spans[0].smin != 0 {
		t.Errorf("Add a span that gets merged with an existing span.")
	}
	if hf.spans[0].smax != 2 {
		t.Errorf("Add a span that gets merged with an existing span.")
	}
	if hf.spans[0].area != area {
		t.Errorf("Add a span that gets merged with an existing span.")
	}
	reset()
	smin = 0
	smax = 1
	if !rcAddSpan(&hf, x, y, smin, smax, area, flagMergeThr) {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0] == nil {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0].smin != smin {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0].smax != smax {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0].area != area {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0].next != nil {
		t.Errorf("Add a span that merges with two spans above and below.")
	}

	smin = 2
	smax = 3
	if !rcAddSpan(&hf, x, y, smin, smax, area, flagMergeThr) {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0].next == nil {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0].next.smin != smin {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0].next.smax != smax {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0].next.area != area {
		t.Errorf("Add a span that merges with two spans above and below.")
	}

	smin = 1
	smax = 2
	if !rcAddSpan(&hf, x, y, smin, smax, area, flagMergeThr) {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0] == nil {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0].smin != 0 {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0].smax != 3 {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0].area != area {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
	if hf.spans[0].next != nil {
		t.Errorf("Add a span that merges with two spans above and below.")
	}
}
func TestRasterizeTriangle(t *testing.T) {
	verts := []float32{
		0, 0, 0,
		1, 0, 0,
		0, 0, -1,
	}
	bmin := make([]float32, 3)
	bmax := make([]float32, 3)
	rcCalcBounds(verts, 3, bmin, bmax)

	cellSize := .5
	cellHeight := .5

	var width int
	var height int

	rcCalcGridSize(bmin, bmax, cellSize, &width, &height)

	var solid RcHeightfield
	assertTrue(t, rcCreateHeightfield(&solid, width, height, bmin, bmax, cellSize, cellHeight), "Rasterize a triangle")

	area := 42
	flagMergeThr := 1

	assertTrue(t, rcRasterizeTriangle(verts[:], verts[3:], verts[6:], area, &solid, flagMergeThr), "Rasterize a triangle")

	assertTrue(t, solid.spans[0+0*width] != nil, "Rasterize a triangle")
	assertTrue(t, solid.spans[1+0*width] == nil, "Rasterize a triangle")
	assertTrue(t, solid.spans[0+1*width] != nil, "Rasterize a triangle")
	assertTrue(t, solid.spans[1+1*width] != nil, "Rasterize a triangle")

	assertTrue(t, solid.spans[0+0*width].smin == 0, "Rasterize a triangle")
	assertTrue(t, solid.spans[0+0*width].smax == 1, "Rasterize a triangle")
	assertTrue(t, solid.spans[0+0*width].area == area, "Rasterize a triangle")
	assertTrue(t, solid.spans[0+0*width].next == nil, "Rasterize a triangle")

	assertTrue(t, solid.spans[0+1*width].smin == 0, "Rasterize a triangle")
	assertTrue(t, solid.spans[0+1*width].smax == 1, "Rasterize a triangle")
	assertTrue(t, solid.spans[0+1*width].area == area, "Rasterize a triangle")
	assertTrue(t, solid.spans[0+1*width].next == nil, "Rasterize a triangle")

	assertTrue(t, solid.spans[1+1*width].smin == 0, "Rasterize a triangle")
	assertTrue(t, solid.spans[1+1*width].smax == 1, "Rasterize a triangle")
	assertTrue(t, solid.spans[1+1*width].area == area, "Rasterize a triangle")
	assertTrue(t, solid.spans[1+1*width].next == nil, "Rasterize a triangle")

}
func TestRasterizeTriangle1(t *testing.T) {
	// This is a minimal repro case for the issue fixed in PR #476 (https://github.com/recastnavigation/recastnavigation/pull/476)

	// create a heightfield
	cellSize := 1.0
	cellHeight := 1.0
	width := 10
	height := 10
	bmin := []float32{0, 0, 0}
	bmax := []float32{10, 10, 10}
	var heightfield RcHeightfield
	assertTrue(t, rcCreateHeightfield(&heightfield, width, height, bmin, bmax, cellSize, cellHeight), "rcRasterizeTriangle overlapping bb but non-overlapping triangle")

	// rasterize a triangle outside of the heightfield.
	area := 42
	flagMergeThr := 1
	verts :=
		[]float32{
			-10.0, 5.5, -10.0,
			-10.0, 5.5, 3,
			3.0, 5.5, -10.0,
		}
	assertTrue(t, rcRasterizeTriangle(verts[0:], verts[3:], verts[6:], area, &heightfield, flagMergeThr), "rcRasterizeTriangle overlapping bb but non-overlapping triangle")

	// ensure that no spans were created
	for x := 0; x < width; x++ {
		for z := 0; z < height; z++ {
			span := heightfield.spans[x+z*heightfield.width]
			assertTrue(t, span == nil, "rcRasterizeTriangle overlapping bb but non-overlapping triangle")
		}
	}
}

func TestRasterizeTriangle2(t *testing.T) {

	verts := []float32{
		5, 0, 0.005,
		5, 0, -0.005,
		-5, 0, 0.005,

		-5, 0, 0.005,
		5, 0, -0.005,
		-5, 0, -0.005,
	}
	bmin := make([]float32, 3)
	bmax := make([]float32, 3)
	rcCalcBounds(verts, 3, bmin, bmax)

	cellSize := 1.0
	cellHeight := 1.0

	var width = 0
	var height = 0

	rcCalcGridSize(bmin, bmax, cellSize, &width, &height)

	var solid RcHeightfield
	assertTrue(t, rcCreateHeightfield(&solid, width, height, bmin, bmax, cellSize, cellHeight), "Skinny triangle along x axis")

	areas := []int{42, 42}
	flagMergeThr := 1
	assertTrue(t, rcRasterizeTriangles1(verts, areas, 2, &solid, flagMergeThr), "Skinny triangle along x axis")

	verts = []float32{
		0.005, 0, 5,
		-0.005, 0, 5,
		0.005, 0, -5,

		0.005, 0, -5,
		-0.005, 0, 5,
		-0.005, 0, -5,
	}
	bmin = make([]float32, 3)
	bmax = make([]float32, 3)
	rcCalcBounds(verts, 3, bmin, bmax)

	cellSize = 1.0
	cellHeight = 1.0

	width = 0
	height = 0

	rcCalcGridSize(bmin, bmax, cellSize, &width, &height)

	solid = RcHeightfield{}
	assertTrue(t, rcCreateHeightfield(&solid, width, height, bmin, bmax, cellSize, cellHeight), "Skinny triangle along z axis")

	areas = []int{42, 42}
	flagMergeThr = 1
	assertTrue(t, rcRasterizeTriangles1(verts, areas, 2, &solid, flagMergeThr), "Skinny triangle along z axis")

}

func TestRasterizeTriangles1(t *testing.T) {

	verts := []float32{
		0, 0, 0,
		1, 0, 0,
		0, 0, -1,
		0, 0, 1,
	}
	tris := []int{
		0, 1, 2,
		0, 3, 1,
	}
	areas := []int{
		1,
		2,
	}
	bmin := make([]float32, 3)
	bmax := make([]float32, 3)
	rcCalcBounds(verts, 4, bmin, bmax)

	cellSize := .5
	cellHeight := .5

	var width int
	var height int

	rcCalcGridSize(bmin, bmax, cellSize, &width, &height)

	var solid RcHeightfield
	assertTrue(t, rcCreateHeightfield(&solid, width, height, bmin, bmax, cellSize, cellHeight), "rcRasterizeTriangles")
	msg := "Rasterize some triangles"
	flagMergeThr := 1
	assertTrue(t, rcRasterizeTriangles(verts, 4, tris, areas, 2, &solid, flagMergeThr), msg)

	assertTrue(t, solid.spans[0+0*width] != nil, msg)
	assertTrue(t, solid.spans[0+1*width] != nil, msg)
	assertTrue(t, solid.spans[0+2*width] != nil, msg)
	assertTrue(t, solid.spans[0+3*width] != nil, msg)
	assertTrue(t, solid.spans[1+0*width] == nil, msg)
	assertTrue(t, solid.spans[1+1*width] != nil, msg)
	assertTrue(t, solid.spans[1+2*width] != nil, msg)
	assertTrue(t, solid.spans[1+3*width] == nil, msg)

	assertTrue(t, solid.spans[0+0*width].smin == 0, msg)
	assertTrue(t, solid.spans[0+0*width].smax == 1, msg)
	assertTrue(t, solid.spans[0+0*width].area == 1, msg)
	assertTrue(t, solid.spans[0+0*width].next == nil, msg)

	assertTrue(t, solid.spans[0+1*width].smin == 0, msg)
	assertTrue(t, solid.spans[0+1*width].smax == 1, msg)
	assertTrue(t, solid.spans[0+1*width].area == 1, msg)
	assertTrue(t, solid.spans[0+1*width].next == nil, msg)

	assertTrue(t, solid.spans[0+2*width].smin == 0, msg)
	assertTrue(t, solid.spans[0+2*width].smax == 1, msg)
	assertTrue(t, solid.spans[0+2*width].area == 2, msg)
	assertTrue(t, solid.spans[0+2*width].next == nil, msg)

	assertTrue(t, solid.spans[0+3*width].smin == 0, msg)
	assertTrue(t, solid.spans[0+3*width].smax == 1, msg)
	assertTrue(t, solid.spans[0+3*width].area == 2, msg)
	assertTrue(t, solid.spans[0+3*width].next == nil, msg)

	assertTrue(t, solid.spans[1+1*width].smin == 0, msg)
	assertTrue(t, solid.spans[1+1*width].smax == 1, msg)
	assertTrue(t, solid.spans[1+1*width].area == 1, msg)
	assertTrue(t, solid.spans[1+1*width].next == nil, msg)

	assertTrue(t, solid.spans[1+2*width].smin == 0, msg)
	assertTrue(t, solid.spans[1+2*width].smax == 1, msg)
	assertTrue(t, solid.spans[1+2*width].area == 2, msg)
	assertTrue(t, solid.spans[1+2*width].next == nil, msg)
}

func TestRasterizeTriangles2(t *testing.T) {

	verts := []float32{
		0, 0, 0,
		1, 0, 0,
		0, 0, -1,
		0, 0, 1,
	}
	areas := []int{
		1,
		2,
	}
	bmin := make([]float32, 3)
	bmax := make([]float32, 3)
	rcCalcBounds(verts, 4, bmin, bmax)

	cellSize := .5
	cellHeight := .5

	var width int
	var height int

	rcCalcGridSize(bmin, bmax, cellSize, &width, &height)

	var solid RcHeightfield
	assertTrue(t, rcCreateHeightfield(&solid, width, height, bmin, bmax, cellSize, cellHeight), "rcRasterizeTriangles")

	flagMergeThr := 1
	msg := "Unsigned short overload"
	utris := []int{
		0, 1, 2,
		0, 3, 1,
	}
	assertTrue(t, rcRasterizeTriangles(verts, 4, utris, areas, 2, &solid, flagMergeThr), msg)

	assertTrue(t, solid.spans[0+0*width] != nil, msg)
	assertTrue(t, solid.spans[0+1*width] != nil, msg)
	assertTrue(t, solid.spans[0+2*width] != nil, msg)
	assertTrue(t, solid.spans[0+3*width] != nil, msg)
	assertTrue(t, solid.spans[1+0*width] == nil, msg)
	assertTrue(t, solid.spans[1+1*width] != nil, msg)
	assertTrue(t, solid.spans[1+2*width] != nil, msg)
	assertTrue(t, solid.spans[1+3*width] == nil, msg)

	assertTrue(t, solid.spans[0+0*width].smin == 0, msg)
	assertTrue(t, solid.spans[0+0*width].smax == 1, msg)
	assertTrue(t, solid.spans[0+0*width].area == 1, msg)
	assertTrue(t, solid.spans[0+0*width].next == nil, msg)

	assertTrue(t, solid.spans[0+1*width].smin == 0, msg)
	assertTrue(t, solid.spans[0+1*width].smax == 1, msg)
	assertTrue(t, solid.spans[0+1*width].area == 1, msg)
	assertTrue(t, solid.spans[0+1*width].next == nil, msg)

	assertTrue(t, solid.spans[0+2*width].smin == 0, msg)
	assertTrue(t, solid.spans[0+2*width].smax == 1, msg)
	assertTrue(t, solid.spans[0+2*width].area == 2, msg)
	assertTrue(t, solid.spans[0+2*width].next == nil, msg)

	assertTrue(t, solid.spans[0+3*width].smin == 0, msg)
	assertTrue(t, solid.spans[0+3*width].smax == 1, msg)
	assertTrue(t, solid.spans[0+3*width].area == 2, msg)
	assertTrue(t, solid.spans[0+3*width].next == nil, msg)

	assertTrue(t, solid.spans[1+1*width].smin == 0, msg)
	assertTrue(t, solid.spans[1+1*width].smax == 1, msg)
	assertTrue(t, solid.spans[1+1*width].area == 1, msg)
	assertTrue(t, solid.spans[1+1*width].next == nil, msg)

	assertTrue(t, solid.spans[1+2*width].smin == 0, msg)
	assertTrue(t, solid.spans[1+2*width].smax == 1, msg)
	assertTrue(t, solid.spans[1+2*width].area == 2, msg)
	assertTrue(t, solid.spans[1+2*width].next == nil, msg)
}
func TestRasterizeTriangles3(t *testing.T) {

	verts := []float32{
		0, 0, 0,
		1, 0, 0,
		0, 0, -1,
		0, 0, 1,
	}
	areas := []int{
		1,
		2,
	}
	bmin := make([]float32, 3)
	bmax := make([]float32, 3)
	rcCalcBounds(verts, 4, bmin, bmax)

	cellSize := .5
	cellHeight := .5

	var width int
	var height int

	rcCalcGridSize(bmin, bmax, cellSize, &width, &height)

	var solid RcHeightfield
	assertTrue(t, rcCreateHeightfield(&solid, width, height, bmin, bmax, cellSize, cellHeight), "rcRasterizeTriangles")

	flagMergeThr := 1
	vertsList := []float32{
		0, 0, 0,
		1, 0, 0,
		0, 0, -1,
		0, 0, 0,
		0, 0, 1,
		1, 0, 0,
	}
	msg := "Triangle list overload"
	assertTrue(t, rcRasterizeTriangles1(vertsList, areas, 2, &solid, flagMergeThr), msg)

	assertTrue(t, solid.spans[0+0*width] != nil, msg)
	assertTrue(t, solid.spans[0+1*width] != nil, msg)
	assertTrue(t, solid.spans[0+2*width] != nil, msg)
	assertTrue(t, solid.spans[0+3*width] != nil, msg)
	assertTrue(t, solid.spans[1+0*width] == nil, msg)
	assertTrue(t, solid.spans[1+1*width] != nil, msg)
	assertTrue(t, solid.spans[1+2*width] != nil, msg)
	assertTrue(t, solid.spans[1+3*width] == nil, msg)

	assertTrue(t, solid.spans[0+0*width].smin == 0, msg)
	assertTrue(t, solid.spans[0+0*width].smax == 1, msg)
	assertTrue(t, solid.spans[0+0*width].area == 1, msg)
	assertTrue(t, solid.spans[0+0*width].next == nil, msg)

	assertTrue(t, solid.spans[0+1*width].smin == 0, msg)
	assertTrue(t, solid.spans[0+1*width].smax == 1, msg)
	assertTrue(t, solid.spans[0+1*width].area == 1, msg)
	assertTrue(t, solid.spans[0+1*width].next == nil, msg)

	assertTrue(t, solid.spans[0+2*width].smin == 0, msg)
	assertTrue(t, solid.spans[0+2*width].smax == 1, msg)
	assertTrue(t, solid.spans[0+2*width].area == 2, msg)
	assertTrue(t, solid.spans[0+2*width].next == nil, msg)

	assertTrue(t, solid.spans[0+3*width].smin == 0, msg)
	assertTrue(t, solid.spans[0+3*width].smax == 1, msg)
	assertTrue(t, solid.spans[0+3*width].area == 2, msg)
	assertTrue(t, solid.spans[0+3*width].next == nil, msg)

	assertTrue(t, solid.spans[1+1*width].smin == 0, msg)
	assertTrue(t, solid.spans[1+1*width].smax == 1, msg)
	assertTrue(t, solid.spans[1+1*width].area == 1, msg)
	assertTrue(t, solid.spans[1+1*width].next == nil, msg)

	assertTrue(t, solid.spans[1+2*width].smin == 0, msg)
	assertTrue(t, solid.spans[1+2*width].smax == 1, msg)
	assertTrue(t, solid.spans[1+2*width].area == 2, msg)
	assertTrue(t, solid.spans[1+2*width].next == nil, msg)
}
