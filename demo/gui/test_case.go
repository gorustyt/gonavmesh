package gui

import "gonavamesh/recast"

type TestType int

const (
	TEST_PATHFIND TestType = iota
	TEST_RAYCAST
)

type Test struct {
	Type         TestType
	spos         [3]float64
	epos         [3]float64
	nspos        [3]float64
	nepos        [3]float64
	radius       float64
	includeFlags int
	excludeFlags int
	expand       bool

	straight  []float64
	nstraight int
	polys     []recast.DtPolyRef
	npolys    int

	findNearestPolyTime  int
	findPathTime         int
	findStraightPathTime int

	next *Test
}
type TestCase struct {
	m_sampleName   string
	m_geomFileName string
	m_tests        *Test
}

func newTestCase() *TestCase {
	return &TestCase{}
}
func (t *TestCase) load(filePath string) bool {
	return true
}

func (t *TestCase) getSampleName() string   { return t.m_sampleName }
func (t *TestCase) getGeomFileName() string { return t.m_geomFileName }

func (t *TestCase) doTests(navmesh *recast.DtNavMesh, navquery recast.DtNavMeshQuery) {

}

func (t *TestCase) handleRender() {

}
func (t *TestCase) handleRenderOverlay(proj, model []float64, view []int) bool {
	return true
}
