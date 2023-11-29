package gui

import (
	"bufio"
	"fmt"
	"gonavamesh/common"
	"gonavamesh/debug_utils"
	"gonavamesh/recast"
	"io"
	"log"
	"os"
	"path/filepath"
	"strings"
)

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
	m_mesh             *rcMeshLoaderObj
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
func (g *InputGeom) loadMesh(path string) bool {
	if g.m_mesh != nil {
		g.m_chunkyMesh = nil
		g.m_mesh = nil
	}
	g.m_offMeshConCount = 0
	g.m_volumeCount = 0

	g.m_mesh = newRcMeshLoaderObj()
	if g.m_mesh == nil {
		log.Printf("loadMesh: Out of memory 'm_mesh'.")
		return false
	}
	if err := g.m_mesh.Load(path); err != nil {
		log.Printf("buildTiledNavigation: Could not load '%s' err:%v", path, err)
		return false
	}

	recast.RcCalcBounds(g.m_mesh.getVerts(), g.m_mesh.getVertCount(), g.m_meshBMin, g.m_meshBMax)

	g.m_chunkyMesh = newRcChunkyTriMesh()
	if !rcCreateChunkyTriMesh(g.m_mesh.getVerts(), g.m_mesh.getTris(), g.m_mesh.getTriCount(), 256, g.m_chunkyMesh) {
		log.Printf("buildTiledNavigation: Failed to build chunky rcMeshLoaderObj.")
		return false
	}

	return true
}
func (g *InputGeom) parseRow(ss []string) {
	switch ss[0] {
	case "f":
		// File name.
		if ss[1] != "" {
			if !g.loadMesh(ss[1]) {
				return
			}
		}
	case "c":
		// Off-mesh connection
		if g.m_offMeshConCount < MAX_OFFMESH_CONNECTIONS {
			v := g.m_offMeshConVerts[g.m_offMeshConCount*3*2:]
			bidir := 0
			area := 0
			flags := 0
			var rad float64
			_, err := fmt.Sscanf(strings.Join(ss[1:], " "), "%f %f %f  %f %f %f %f %d %d %d",
				&v[0], &v[1], &v[2], &v[3], &v[4], &v[5], &rad, &bidir, &area, &flags)
			g.m_offMeshConRads[g.m_offMeshConCount] = rad
			g.m_offMeshConDirs[g.m_offMeshConCount] = bidir
			g.m_offMeshConAreas[g.m_offMeshConCount] = area
			g.m_offMeshConFlags[g.m_offMeshConCount] = flags
			g.m_offMeshConCount++
			if err != nil {
				log.Println(err)
			}
		}
	case "v":

		// Convex volumes
		if g.m_volumeCount < MAX_VOLUMES {
			vol := g.m_volumes[g.m_volumeCount]
			g.m_volumeCount++
			_, err := fmt.Sscanf(strings.Join(ss[1:], " "), "%d %d %f %f", &vol.nverts, vol.area, vol.hmin, vol.hmax)
			for i := 0; i < vol.nverts; i++ {
				_, err = fmt.Sscanf(strings.Join(ss[1:], " "), "%f %f %f", &vol.verts[i*3+0], &vol.verts[i*3+1], &vol.verts[i*3+2])
			}
			if err != nil {
				log.Println(err)
			}
		}
	case "s":
		// Settings
		g.m_hasBuildSettings = true
		_, err := fmt.Sscanf(strings.Join(ss[1:], " "), "%f %f %f %f %f %f %f %f %f %f %f %f %f %d %f %f %f %f %f %f %f",
			g.m_buildSettings.cellSize,
			g.m_buildSettings.cellHeight,
			g.m_buildSettings.agentHeight,
			g.m_buildSettings.agentRadius,
			g.m_buildSettings.agentMaxClimb,
			g.m_buildSettings.agentMaxSlope,
			g.m_buildSettings.regionMinSize,
			g.m_buildSettings.regionMergeSize,
			g.m_buildSettings.edgeMaxLen,
			g.m_buildSettings.edgeMaxError,
			g.m_buildSettings.vertsPerPoly,
			g.m_buildSettings.detailSampleDist,
			g.m_buildSettings.detailSampleMaxError,
			g.m_buildSettings.partitionType,
			g.m_buildSettings.navMeshBMin[0],
			g.m_buildSettings.navMeshBMin[1],
			g.m_buildSettings.navMeshBMin[2],
			g.m_buildSettings.navMeshBMax[0],
			g.m_buildSettings.navMeshBMax[1],
			g.m_buildSettings.navMeshBMax[2],
			g.m_buildSettings.tileSize)
		if err != nil {
			log.Println(err)
		}
	}
}
func (g *InputGeom) loadGeomSet(path string) bool {
	f, err := os.Open(path)
	if err != nil {
		log.Println(err)
		panic("")
	}
	defer f.Close()
	reader := bufio.NewReader(f)
	line := ""
	for {
		lines, prefix, err := reader.ReadLine()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Println(err)
			panic("")
		}
		if prefix {
			line += string(lines)
		} else {
			g.parseRow(strings.Split(line, " "))
		}
	}

	return true
}
func (g *InputGeom) load(path string) bool {
	extension := filepath.Ext(path)
	if extension == ".gset" {
		return g.loadGeomSet(path)
	}

	if extension == ".obj" {
		return g.loadMesh(path)
	}
	return false
}

func (g *InputGeom) saveGeomSet(settings *BuildSettings) bool {
	if g.m_mesh == nil {
		return false
	}

	// Change extension
	path := g.m_mesh.getFileName() + ".gset"
	fp, err := os.Open(path)
	defer fp.Close()
	common.AssertTrue(err == nil)
	// Store rcMeshLoaderObj filename.
	_, err = fmt.Fscanln(fp, "f %s", g.m_mesh.getFileName())
	common.AssertTrue(err == nil)
	// Store settings if any
	if settings != nil {
		_, err := fmt.Fscanln(fp,
			"s %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f %f %f %f %f %f %f",
			settings.cellSize,
			settings.cellHeight,
			settings.agentHeight,
			settings.agentRadius,
			settings.agentMaxClimb,
			settings.agentMaxSlope,
			settings.regionMinSize,
			settings.regionMergeSize,
			settings.edgeMaxLen,
			settings.edgeMaxError,
			settings.vertsPerPoly,
			settings.detailSampleDist,
			settings.detailSampleMaxError,
			settings.partitionType,
			settings.navMeshBMin[0],
			settings.navMeshBMin[1],
			settings.navMeshBMin[2],
			settings.navMeshBMax[0],
			settings.navMeshBMax[1],
			settings.navMeshBMax[2],
			settings.tileSize)
		if err != nil {
			log.Println(err)
			return false
		}
	}

	// Store off-rcMeshLoaderObj links.
	for i := 0; i < g.m_offMeshConCount; i++ {
		v := g.m_offMeshConVerts[i*3*2:]
		rad := g.m_offMeshConRads[i]
		bidir := g.m_offMeshConDirs[i]
		area := g.m_offMeshConAreas[i]
		flags := g.m_offMeshConFlags[i]
		_, err := fmt.Fscanln(fp, "c %f %f %f  %f %f %f  %f %d %d %d",
			v[0], v[1], v[2], v[3], v[4], v[5], rad, bidir, area, flags)
		if err != nil {
			log.Println(err)
			return false
		}
	}

	// Convex volumes
	for i := 0; i < g.m_volumeCount; i++ {
		vol := g.m_volumes[i]
		_, err := fmt.Fscanln(fp, "v %d %d %f %f", vol.nverts, vol.area, vol.hmin, vol.hmax)
		if err != nil {
			log.Println(err)
			return false
		}
		for j := 0; j < vol.nverts; j++ {
			_, err := fmt.Fscanln(fp, "%f %f %f", vol.verts[j*3+0], vol.verts[j*3+1], vol.verts[j*3+2])
			if err != nil {
				log.Println(err)
				return false
			}
		}

	}

	return true
}

// / Method to return static rcMeshLoaderObj data.
func (g *InputGeom) getMesh() *rcMeshLoaderObj   { return g.m_mesh }
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

func (g *InputGeom) drawConvexVolumes(dd debug_utils.DuDebugDraw, hilights ...bool) {
	var hilight bool
	if len(hilights) > 0 {
		hilight = hilights[0]
	}
	_ = hilight
	dd.DepthMask(false)

	dd.Begin(debug_utils.DU_DRAW_TRIS)

	for i := 0; i < g.m_volumeCount; i++ {
		vol := g.m_volumes[i]
		col := debug_utils.DuTransCol(dd.AreaToCol(vol.area), 32)
		j := 0
		k := vol.nverts - 1
		for j < vol.nverts {
			va := common.GetVs3(vol.verts, k)
			vb := common.GetVs3(vol.verts, j)

			dd.Vertex1(vol.verts[0], vol.hmax, vol.verts[2], col)
			dd.Vertex1(vb[0], vol.hmax, vb[2], col)
			dd.Vertex1(va[0], vol.hmax, va[2], col)

			dd.Vertex1(va[0], vol.hmin, va[2], debug_utils.DuDarkenCol(col))
			dd.Vertex1(va[0], vol.hmax, va[2], col)
			dd.Vertex1(vb[0], vol.hmax, vb[2], col)

			dd.Vertex1(va[0], vol.hmin, va[2], debug_utils.DuDarkenCol(col))
			dd.Vertex1(vb[0], vol.hmax, vb[2], col)
			dd.Vertex1(vb[0], vol.hmin, vb[2], debug_utils.DuDarkenCol(col))
			k = j
			j++
		}
	}

	dd.End()

	dd.Begin(debug_utils.DU_DRAW_LINES, 2.0)
	for i := 0; i < g.m_volumeCount; i++ {
		vol := g.m_volumes[i]
		col := debug_utils.DuTransCol(dd.AreaToCol(vol.area), 220)
		j := 0
		k := vol.nverts - 1
		for j < vol.nverts {
			va := common.GetVs3(vol.verts, k)
			vb := common.GetVs3(vol.verts, j)
			dd.Vertex1(va[0], vol.hmin, va[2], debug_utils.DuDarkenCol(col))
			dd.Vertex1(vb[0], vol.hmin, vb[2], debug_utils.DuDarkenCol(col))
			dd.Vertex1(va[0], vol.hmax, va[2], col)
			dd.Vertex1(vb[0], vol.hmax, vb[2], col)
			dd.Vertex1(va[0], vol.hmin, va[2], debug_utils.DuDarkenCol(col))
			dd.Vertex1(va[0], vol.hmax, va[2], col)
			k = j
			j++
		}
	}
	dd.End()

	dd.Begin(debug_utils.DU_DRAW_POINTS, 3.0)
	for i := 0; i < g.m_volumeCount; i++ {
		vol := g.m_volumes[i]
		col := debug_utils.DuDarkenCol(debug_utils.DuTransCol(dd.AreaToCol(vol.area), 220))
		for j := 0; j < vol.nverts; j++ {
			dd.Vertex1(vol.verts[j*3+0], vol.verts[j*3+1]+0.1, vol.verts[j*3+2], col)
			dd.Vertex1(vol.verts[j*3+0], vol.hmin, vol.verts[j*3+2], col)
			dd.Vertex1(vol.verts[j*3+0], vol.hmax, vol.verts[j*3+2], col)
		}
	}
	dd.End()

	dd.DepthMask(true)
}
func (g *InputGeom) raycastMesh(src, dst []float64, tmin *float64) bool {
	// Prune hit ray.
	var btmin, btmax float64
	if !isectSegAABB(src, dst, g.m_meshBMin, g.m_meshBMax, btmin, btmax) {
		return false
	}

	var p, q [2]float64
	p[0] = src[0] + (dst[0]-src[0])*btmin
	p[1] = src[2] + (dst[2]-src[2])*btmin
	q[0] = src[0] + (dst[0]-src[0])*btmax
	q[1] = src[2] + (dst[2]-src[2])*btmax

	cid := make([]int, 512)
	ncid := rcGetChunksOverlappingSegment(g.m_chunkyMesh, p, q, cid, 512)
	if ncid == 0 {
		return false
	}

	*tmin = 1.0
	hit := false
	verts := g.m_mesh.getVerts()

	for i := 0; i < ncid; i++ {
		node := g.m_chunkyMesh.nodes[cid[i]]
		tris := g.m_chunkyMesh.tris[node.i*3:]
		ntris := node.n

		for j := 0; j < ntris*3; j += 3 {
			t := 1.0
			if intersectSegmentTriangle(src, dst,
				common.GetVs3(verts, tris[j]),
				common.GetVs3(verts, tris[j+1]),
				common.GetVs3(verts, tris[j+2]), t) {
				if t < *tmin {
					*tmin = t
				}
				hit = true
			}
		}
	}

	return hit
}
func intersectSegmentTriangle(sp, sq, a, b, c []float64, t float64) bool {
	var v, w float64
	ab := make([]float64, 3)
	ac := make([]float64, 3)
	qp := make([]float64, 3)
	ap := make([]float64, 3)
	norm := make([]float64, 3)
	e := make([]float64, 3)
	common.Vsub(ab, b, a)
	common.Vsub(ac, c, a)
	common.Vsub(qp, sp, sq)

	// Compute triangle normal. Can be precalculated or cached if
	// intersecting multiple segments against the same triangle
	common.Vcross(norm, ab, ac)

	// Compute denominator d. If d <= 0, segment is parallel to or points
	// away from triangle, so exit early
	d := common.Vdot(qp, norm)
	if d <= 0.0 {
		return false
	}

	// Compute intersection t value of pq with plane of triangle. A ray
	// intersects iff 0 <= t. Segment intersects iff 0 <= t <= 1. Delay
	// dividing by d until intersection has been found to pierce triangle
	common.Vsub(ap, sp, a)
	t = common.Vdot(ap, norm)
	if t < 0.0 {
		return false
	}
	if t > d {
		return false
	} // For segment; exclude this code line for a ray test

	// Compute barycentric coordinate components and test if within bounds
	common.Vcross(e, qp, ap)
	v = common.Vdot(ac, e)
	if v < 0.0 || v > d {
		return false
	}
	w = -common.Vdot(ab, e)
	if w < 0.0 || v+w > d {
		return false
	}

	// Segment/ray intersects triangle. Perform delayed division
	t /= d

	return true
}

func (g *InputGeom) drawOffMeshConnections(dd debug_utils.DuDebugDraw, hilights ...bool) {
	var hilight bool
	if len(hilights) > 0 {
		hilight = hilights[0]
	}
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
	// Bounds of the area to rcMeshLoaderObj
	navMeshBMin [3]float64
	navMeshBMax [3]float64
	// Size of the tiles in voxels
	tileSize float64
}
