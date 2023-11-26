package gui

import (
	"gonavamesh/debug_utils"
	"gonavamesh/recast"
)

type TileFlags struct {
	flags  []int
	nflags int
	base   recast.DtPolyRef
}

func (t *TileFlags) purge() {
	t.flags = nil
}

type NavmeshFlags struct {
	m_nav    recast.IDtNavMesh
	m_tiles  []*TileFlags
	m_ntiles int
}

func newNavmeshFlags(nav recast.IDtNavMesh) *NavmeshFlags {
	n := &NavmeshFlags{m_nav: nav}
	n.m_ntiles = nav.GetMaxTiles()
	if n.m_ntiles == 0 {
		return n
	}

	n.m_tiles = make([]*TileFlags, n.m_ntiles)
	for i := range n.m_tiles {
		n.m_tiles[i] = &TileFlags{}
	}
	// Alloc flags for each tile.
	for i := 0; i < nav.GetMaxTiles(); i++ {
		tile := nav.GetTile(i)
		if tile.Header == nil {
			continue
		}
		tf := n.m_tiles[i]
		tf.nflags = tile.Header.PolyCount
		tf.base = nav.GetPolyRefBase(tile)
		if tf.nflags > 0 {
			tf.flags = make([]int, tf.nflags)
		}
	}

	return n
}

func (n *NavmeshFlags) clearAllFlags() {
	for i := 0; i < n.m_ntiles; i++ {
		tf := n.m_tiles[i]
		tf.flags = make([]int, tf.nflags)
	}
}

func (n *NavmeshFlags) getFlags(ref recast.DtPolyRef) int {

	// Assume the ref is valid, no bounds checks.
	_, it, ip := n.m_nav.DecodePolyId(ref)
	return n.m_tiles[it].flags[ip]
}

func (n *NavmeshFlags) setFlags(ref recast.DtPolyRef, flags int) {
	// Assume the ref is valid, no bounds checks.
	_, it, ip := n.m_nav.DecodePolyId(ref)
	n.m_tiles[it].flags[ip] = flags
}

func FloodNavmesh(nav recast.IDtNavMesh, flags *NavmeshFlags, start recast.DtPolyRef, flag int) {
	// If already visited, skip.
	if flags.getFlags(start) != 0 {
		return
	}

	flags.setFlags(start, flag)

	openList := recast.NewStack(func() recast.DtPolyRef {
		return 0
	})
	openList.Push(start)

	for openList.Len() != 0 {
		ref := openList.Pop()

		// Get current poly and tile.
		// The API input has been checked already, skip checking internal data.
		tile, poly := nav.GetTileAndPolyByRefUnsafe(ref)

		// Visit linked polygons.
		for i := poly.FirstLink; i != recast.DT_NULL_LINK; i = tile.Links[i].Next {
			neiRef := tile.Links[i].Ref
			// Skip invalid and already visited.
			if neiRef == 0 || flags.getFlags(neiRef) != 0 {
				continue
			}

			// Mark as visited
			flags.setFlags(neiRef, flag)
			// Visit neighbours
			openList.Push(neiRef)
		}
	}
}

func disableUnvisitedPolys(nav recast.IDtNavMesh, flags *NavmeshFlags) {
	for i := 0; i < nav.GetMaxTiles(); i++ {
		tile := nav.GetTile(i)
		if tile.Header == nil {
			continue
		}
		base := nav.GetPolyRefBase(tile)
		for j := 0; j < tile.Header.PolyCount; j++ {
			ref := base | recast.DtPolyRef(j)
			if flags.getFlags(ref) == 0 {
				f, _ := nav.GetPolyFlags(ref)
				nav.SetPolyFlags(ref, f|SAMPLE_POLYFLAGS_DISABLED)
			}
		}
	}
}

type NavMeshPruneTool struct {
	m_sample *Sample

	m_flags     *NavmeshFlags
	m_hitPos    []float64
	m_hitPosSet bool
	gs          *guiState
}

func newNavMeshPruneTool(gs *guiState) *NavMeshPruneTool {
	return &NavMeshPruneTool{
		m_hitPos: make([]float64, 3),
		gs:       gs,
	}
}
func (n *NavMeshPruneTool) init(sample *Sample) {
	n.m_sample = sample
}
func (n *NavMeshPruneTool) Type() int { return TOOL_NAVMESH_PRUNE }
func (n *NavMeshPruneTool) reset() {
	n.m_hitPosSet = false
	n.m_flags = nil
}

func (n *NavMeshPruneTool) handleMenu() {
	nav := n.m_sample.getNavMesh()
	if nav == nil {
		return
	}
	if n.m_flags == nil {
		return
	}

	if n.gs.imguiButton("Clear Selection") {
		n.m_flags.clearAllFlags()
	}

	if n.gs.imguiButton("Prune Unselected") {
		disableUnvisitedPolys(nav, n.m_flags)
		n.m_flags = nil
	}
}

func (n *NavMeshPruneTool) handleClick(s, p []float64, shift bool) {
	if n.m_sample == nil {
		return
	}
	geom := n.m_sample.getInputGeom()
	if geom == nil {
		return
	}
	nav := n.m_sample.getNavMesh()
	if nav == nil {
		return
	}
	query := n.m_sample.getNavMeshQuery()
	if query == nil {
		return
	}

	copy(n.m_hitPos, p)
	n.m_hitPosSet = true

	if n.m_flags == nil {
		n.m_flags = newNavmeshFlags(nav)
	}

	halfExtents := []float64{2, 4, 2}
	ref, _ := query.FindNearestPoly(p, halfExtents, &recast.DtQueryFilter{}, []float64{})

	FloodNavmesh(nav, n.m_flags, ref, 1)
}

func (n *NavMeshPruneTool) handleToggle() {
}

func (n *NavMeshPruneTool) handleStep() {
}

func (n *NavMeshPruneTool) handleUpdate(dt float64) {
}

func (n *NavMeshPruneTool) handleRender() {
	dd := n.m_sample.getDebugDraw()

	if n.m_hitPosSet {
		s := n.m_sample.getAgentRadius()
		col := debug_utils.DuRGBA(255, 255, 255, 255)
		dd.Begin(debug_utils.DU_DRAW_LINES)
		dd.Vertex1(n.m_hitPos[0]-s, n.m_hitPos[1], n.m_hitPos[2], col)
		dd.Vertex1(n.m_hitPos[0]+s, n.m_hitPos[1], n.m_hitPos[2], col)
		dd.Vertex1(n.m_hitPos[0], n.m_hitPos[1]-s, n.m_hitPos[2], col)
		dd.Vertex1(n.m_hitPos[0], n.m_hitPos[1]+s, n.m_hitPos[2], col)
		dd.Vertex1(n.m_hitPos[0], n.m_hitPos[1], n.m_hitPos[2]-s, col)
		dd.Vertex1(n.m_hitPos[0], n.m_hitPos[1], n.m_hitPos[2]+s, col)
		dd.End()
	}

	nav := n.m_sample.getNavMesh()
	if n.m_flags != nil && nav != nil {
		for i := 0; i < nav.GetMaxTiles(); i++ {
			tile := nav.GetTile(i)
			if tile.Header == nil {
				continue
			}
			base := nav.GetPolyRefBase(tile)
			for j := 0; j < tile.Header.PolyCount; j++ {
				ref := base | recast.DtPolyRef(j)
				if n.m_flags.getFlags(ref) != 0 {
					debug_utils.DuDebugDrawNavMeshPoly(dd, nav, ref, debug_utils.DuRGBA(255, 255, 255, 128))
				}
			}
		}
	}

}

func (n *NavMeshPruneTool) handleRenderOverlay(proj, model []float64, view []int) {
	// Tool help
	h := view[3]
	n.gs.imguiDrawText(280, h-40, IMGUI_ALIGN_LEFT, "LMB: Click fill area.", imguiRGBA(255, 255, 255, 192))
}
