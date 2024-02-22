package mesh

import (
	"github.com/gorustyt/gonavmesh/common"
	"github.com/gorustyt/gonavmesh/detour"
	"github.com/gorustyt/gonavmesh/recast"
)

type SampleDebug struct {
	*Sample
	m_chf         *recast.RcCompactHeightfield
	m_cset        *recast.RcContourSet
	m_pmesh       *recast.RcPolyMesh
	m_halfExtents [3]float32
	m_center      [3]float32
	m_bmin        [3]float32
	m_bmax        [3]float32
	m_ref         detour.DtPolyRef
}

func newSampleDebug() *SampleDebug {
	return &SampleDebug{
		Sample: &Sample{},
	}
}

func (s *SampleDebug) handleSettings() {}

func (s *SampleDebug) handleDebugMode() {

}
func (s *SampleDebug) handleClick(ss, p []float32, shift bool) {
	if s.m_tool != nil {
		s.m_tool.handleClick(ss, p, shift)
	}
}
func (s *SampleDebug) handleToggle() {
	if s.m_tool != nil {
		s.m_tool.handleToggle()
	}

}
func (s *SampleDebug) handleRender() {

}
func (s *SampleDebug) handleRenderOverlay(proj, model []float32, view []int) {

}
func (s *SampleDebug) handleMeshChanged(geom *InputGeom) {
	s.m_geom = geom
}
func (s *SampleDebug) handleBuild() {

}

func (s *SampleDebug) getBoundsMin() []float32 {
	if s.m_cset != nil {
		return common.SliceTToSlice[float32, float32](s.m_cset.Bmin)
	}

	if s.m_chf != nil {
		return common.SliceTToSlice[float32, float32](s.m_chf.Bmin[:])
	}

	if s.m_navMesh != nil {
		return s.m_bmin[:]
	}

	return nil
}
func (s *SampleDebug) getBoundsMax() []float32 {
	if s.m_cset != nil {
		return common.SliceTToSlice[float32, float32](s.m_cset.Bmax)
	}

	if s.m_chf != nil {
		return common.SliceTToSlice[float32, float32](s.m_chf.Bmax[:])
	}

	if s.m_navMesh != nil {
		return s.m_bmax[:]
	}

	return nil
}
