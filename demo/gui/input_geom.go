package gui

import "gonavamesh/debug_utils"

const MAX_CONVEXVOL_PTS = 12

const MAX_OFFMESH_CONNECTIONS = 256

const MAX_VOLUMES = 256

type ConvexVolume struct {
	verts      []float64
	hmin, hmax float64
	nverts     int
	area       int
}

func newConvexVolume() *ConvexVolume {
	return &ConvexVolume{
		verts: make([]float64, MAX_CONVEXVOL_PTS*3),
	}
}

type InputGeom struct {
	m_chunkyMesh       *rcChunkyTriMesh
	m_mesh             *mesh
	m_meshBMin         []float64
	m_meshBMax         []float64
	m_buildSettings    *BuildSettings
	m_hasBuildSettings bool

	/// @name Off-Mesh connections.
	///@{
	m_offMeshConVerts []float64
	m_offMeshConRads  []float64
	m_offMeshConDirs  []int
	m_offMeshConAreas []int
	m_offMeshConFlags []int
	m_offMeshConId    []int
	m_offMeshConCount int
	///@}

	/// @name Convex Volumes.
	///@{
	m_volumes     []*ConvexVolume
	m_volumeCount int
	///@}
}

func newInputGeom() *InputGeom {
	return &InputGeom{
		m_meshBMin:        make([]float64, 3),
		m_meshBMax:        make([]float64, 3),
		m_offMeshConDirs:  make([]int, MAX_OFFMESH_CONNECTIONS),
		m_offMeshConAreas: make([]int, MAX_OFFMESH_CONNECTIONS),
		m_offMeshConFlags: make([]int, MAX_OFFMESH_CONNECTIONS),
		m_offMeshConId:    make([]int, MAX_OFFMESH_CONNECTIONS),
		m_offMeshConRads:  make([]float64, MAX_OFFMESH_CONNECTIONS),
		m_offMeshConVerts: make([]float64, MAX_OFFMESH_CONNECTIONS*3*2),
		m_volumes:         make([]*ConvexVolume, MAX_VOLUMES),
	}
}
func (g *InputGeom) loadMesh(l *logger, filepath string) {
	panic("impl")
}
func (g *InputGeom) loadGeomSet(l *logger, filepath string) bool {
	panic("impl")
}

// / Method to return static mesh data.
func (g *InputGeom) getMesh() *mesh              { return g.m_mesh }
func (g *InputGeom) getMeshBoundsMin() []float64 { return g.m_meshBMin }
func (g *InputGeom) getMeshBoundsMax() []float64 { return g.m_meshBMax }
func (g *InputGeom) getNavMeshBoundsMin() []float64 {
	if g.m_hasBuildSettings {
		return g.m_buildSettings.navMeshBMin[:]
	}
	return g.m_meshBMin
}
func (g *InputGeom) getNavMeshBoundsMax() []float64 {
	if g.m_hasBuildSettings {
		return g.m_buildSettings.navMeshBMax[:]
	}
	return g.m_meshBMax
}
func (g *InputGeom) getChunkyMesh() *rcChunkyTriMesh { return g.m_chunkyMesh }
func (g *InputGeom) getBuildSettings() *BuildSettings {
	if g.m_hasBuildSettings {
		return g.m_buildSettings
	}
	return nil
}

func (g *InputGeom) addConvexVolume(verts []float64, nverts int, minh, maxh float64, area int) {
	if g.m_volumeCount >= MAX_VOLUMES {
		return
	}
	vol := g.m_volumes[g.m_volumeCount]
	if vol == nil {
		vol = newConvexVolume()
		g.m_volumes[g.m_volumeCount] = vol
	}
	g.m_volumeCount++
	copy(vol.verts, verts[:3*nverts])
	vol.hmin = minh
	vol.hmax = maxh
	vol.nverts = nverts
	vol.area = area
}

func (g *InputGeom) drawOffMeshConnections(dd debug_utils.DuDebugDraw, hilight bool) {
	conColor := debug_utils.DuRGBA(192, 0, 128, 192)
	baseColor := debug_utils.DuRGBA(0, 0, 0, 64)
	dd.DepthMask(false)

	dd.Begin(debug_utils.DU_DRAW_LINES, 2.0)
	for i := 0; i < g.m_offMeshConCount; i++ {
		v := g.m_offMeshConVerts[i*3*2:]

		dd.Vertex1(v[0], v[1], v[2], baseColor)
		dd.Vertex1(v[0], v[1]+0.2, v[2], baseColor)

		dd.Vertex1(v[3], v[4], v[5], baseColor)
		dd.Vertex1(v[3], v[4]+0.2, v[5], baseColor)

		debug_utils.DuAppendCircle(dd, v[0], v[1]+0.1, v[2], g.m_offMeshConRads[i], baseColor)
		debug_utils.DuAppendCircle(dd, v[3], v[4]+0.1, v[5], g.m_offMeshConRads[i], baseColor)

		if hilight {
			tmp := 0.0
			if g.m_offMeshConDirs[i]&1 > 0 {
				tmp = 0.6
			}

			debug_utils.DuAppendArc(dd, v[0], v[1], v[2], v[3], v[4], v[5], 0.25, tmp, 0.6, conColor)
		}
	}
	dd.End()

	dd.DepthMask(true)
}

// / @name Off-Mesh connections.
// /@{
func (g *InputGeom) getOffMeshConnectionCount() int       { return g.m_offMeshConCount }
func (g *InputGeom) getOffMeshConnectionVerts() []float64 { return g.m_offMeshConVerts }
func (g *InputGeom) getOffMeshConnectionRads() []float64  { return g.m_offMeshConRads }
func (g *InputGeom) getOffMeshConnectionDirs() []int      { return g.m_offMeshConDirs }
func (g *InputGeom) getOffMeshConnectionAreas() []int     { return g.m_offMeshConAreas }
func (g *InputGeom) getOffMeshConnectionFlags() []int     { return g.m_offMeshConFlags }
func (g *InputGeom) getOffMeshConnectionId() []int        { return g.m_offMeshConId }
func (g *InputGeom) addOffMeshConnection(spos, epos []float64, rad float64,
	bidir int, area int, flags int) {
	if g.m_offMeshConCount >= MAX_OFFMESH_CONNECTIONS {
		return
	}
	v := g.m_offMeshConVerts[g.m_offMeshConCount*3*2:]
	g.m_offMeshConRads[g.m_offMeshConCount] = rad
	g.m_offMeshConDirs[g.m_offMeshConCount] = bidir
	g.m_offMeshConAreas[g.m_offMeshConCount] = area
	g.m_offMeshConFlags[g.m_offMeshConCount] = flags
	g.m_offMeshConId[g.m_offMeshConCount] = 1000 + g.m_offMeshConCount
	copy(v[0:3], spos)
	copy(v[3:6], epos)
	g.m_offMeshConCount++
}

func (g *InputGeom) deleteOffMeshConnection(i int) {
	g.m_offMeshConCount--
	src := g.m_offMeshConVerts[g.m_offMeshConCount*3*2:]
	dst := g.m_offMeshConVerts[i*3*2:]
	copy(dst[0:3], src[0:3])
	copy(dst[3:6], src[3:5])
	g.m_offMeshConRads[i] = g.m_offMeshConRads[g.m_offMeshConCount]
	g.m_offMeshConDirs[i] = g.m_offMeshConDirs[g.m_offMeshConCount]
	g.m_offMeshConAreas[i] = g.m_offMeshConAreas[g.m_offMeshConCount]
	g.m_offMeshConFlags[i] = g.m_offMeshConFlags[g.m_offMeshConCount]
}

// / @name Box Volumes.
// /@{
func (g *InputGeom) getConvexVolumeCount() int         { return g.m_volumeCount }
func (g *InputGeom) getConvexVolumes() []*ConvexVolume { return g.m_volumes }
func (g *InputGeom) deleteConvexVolume(index int) {
	g.m_volumeCount--
	g.m_volumes[index] = g.m_volumes[g.m_volumeCount]
}

type BuildSettings struct {
	// Cell size in world units
	cellSize float64
	// Cell height in world units
	cellHeight float64
	// Agent height in world units
	agentHeight float64
	// Agent radius in world units
	agentRadius float64
	// Agent max climb in world units
	agentMaxClimb float64
	// Agent max slope in degrees
	agentMaxSlope float64
	// Region minimum size in voxels.
	// regionMinSize = sqrt(regionMinArea)
	regionMinSize float64
	// Region merge size in voxels.
	// regionMergeSize = sqrt(regionMergeArea)
	regionMergeSize float64
	// Edge max length in world units
	edgeMaxLen float64
	// Edge max error in voxels
	edgeMaxError float64
	vertsPerPoly float64
	// Detail sample distance in voxels
	detailSampleDist float64
	// Detail sample max error in voxel heights.
	detailSampleMaxError float64
	// Partition type, see SamplePartitionType
	partitionType int
	// Bounds of the area to mesh
	navMeshBMin [3]float64
	navMeshBMax [3]float64
	// Size of the tiles in voxels
	tileSize float64
}
