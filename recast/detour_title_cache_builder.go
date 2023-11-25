package recast

import (
	"gonavamesh/common"
	"math"
)

const DT_TILECACHE_MAGIC = 'D'<<24 | 'T'<<16 | 'L'<<8 | 'R' ///< 'DTLR';
const DT_TILECACHE_VERSION = 1

const DT_TILECACHE_NULL_AREA = 0

const DT_TILECACHE_WALKABLE_AREA = 63

const DT_TILECACHE_NULL_IDX = 0xffff

const MAX_VERTS_PER_POLY = 6 // TODO: use the DT_VERTS_PER_POLYGON
const MAX_REM_EDGES = 48     // TODO: make this an expression.

func getDirOffsetX(dir int) int {
	offset := [4]int{-1, 0, 1, 0}
	return offset[dir&0x03]
}

func getDirOffsetY(dir int) int {
	offset := [4]int{0, 1, 0, -1}
	return offset[dir&0x03]
}

type DtTileCacheLayerHeader struct {
	Magic                  int ///< Data magic
	Version                int ///< Data version
	Tx, Ty, Tlayer         int
	Bmin, Bmax             [3]float64
	Hmin, Hmax             int ///< Height min/max range
	Width, Height          int ///< Dimension of the layer.
	Minx, Maxx, Miny, Maxy int ///< Usable sub-region.
}

type DtTileCacheLayer struct {
	Header   *DtTileCacheLayerHeader
	RegCount int ///< Region count.
	Heights  []int
	Areas    []int
	Cons     []int
	Regs     []int
}

type DtTileCacheContour struct {
	Nverts int
	Verts  []int
	Reg    int
	Area   int
}

type DtTileCacheContourSet struct {
	Nconts int
	Conts  []*DtTileCacheContour
}

type DtTileCachePolyMesh struct {
	Nvp    int
	Nverts int   ///< Number of vertices.
	Npolys int   ///< Number of polygons.
	Verts  []int ///< Vertices of the mesh, 3 elements per vertex.
	Polys  []int ///< Polygons of the mesh, nvp*2 elements per polygon.
	Flags  []int ///< Per polygon flags.
	Areas  []int ///< Area ID of polygons.
}

type dtLayerSweepSpan struct {
	ns  int // number samples
	id  int // region id
	nei int // neighbour id
}

const DT_LAYER_MAX_NEIS = 16

type dtLayerMonotoneRegion struct {
	area   int
	neis   [DT_LAYER_MAX_NEIS]int
	nneis  int
	regId  int
	areaId int
}

type dtTempContour struct {
	verts  []int
	nverts int
	cverts int
	poly   []int
	npoly  int
	cpoly  int
}

func newDtTempContour(vbuf []int, nvbuf int, pbuf []int, npbuf int) *dtTempContour {
	d := &dtTempContour{
		verts: vbuf, cverts: nvbuf,
		poly: pbuf, cpoly: npbuf,
	}
	return d
}

type dtTileCacheCompressor interface {
	maxCompressedSize(bufferSize int) int
	compress(buffer []byte, bufferSize int,
		compressed, maxCompressedSize int, compressedSize *int) DtStatus
	decompress(compressed []byte, compressedSize int,
		buffer []byte, maxBufferSizeint, bufferSize *int) DtStatus
}

type dtTileCacheMeshProcess interface {
	process(params *DtNavMeshCreateParams, polyAreas []int, polyFlags []int)
}

func overlapRangeExl(amin, amax,
	bmin, bmax int) bool {
	if amin >= bmax || amax <= bmin {
		return false
	}
	return true
}

func addUniqueLast(a []int, an *int, v int) {
	n := *an
	if n > 0 && a[n-1] == v {
		return
	}
	a[*an] = v
	*an++
}

func isConnected(layer *DtTileCacheLayer,
	ia, ib int, walkableClimb int) bool {
	if layer.Areas[ia] != layer.Areas[ib] {
		return false
	}
	if common.Abs(layer.Heights[ia]-layer.Heights[ib]) > walkableClimb {
		return false
	}
	return true
}

func canMerge(oldRegId, newRegId int, regs []*dtLayerMonotoneRegion, nregs int) bool {
	count := 0
	for i := 0; i < nregs; i++ {
		reg := regs[i]
		if reg.regId != oldRegId {
			continue
		}
		nnei := reg.nneis
		for j := 0; j < nnei; j++ {
			if regs[reg.neis[j]].regId == newRegId {
				count++
			}

		}
	}
	return count == 1
}
func dtBuildTileCacheRegions(
	layer *DtTileCacheLayer,
	walkableClimb int) DtStatus {

	w := layer.Header.Width
	h := layer.Header.Height
	layer.Regs = make([]int, w*h)
	for i := range layer.Regs {
		layer.Regs[i] = 0xff
	}

	nsweeps := w
	sweeps := make([]*dtLayerSweepSpan, nsweeps)
	for i := range sweeps {
		sweeps[i] = &dtLayerSweepSpan{}
	}
	// Partition walkable area into monotone regions.
	prevCount := make([]int, 256)
	regId := 0

	for y := 0; y < h; y++ {
		if regId > 0 {
			for i := 0; i < regId; i++ {
				prevCount[i] = 0
			}
		}

		sweepId := 0

		for x := 0; x < w; x++ {
			idx := x + y*w
			if layer.Areas[idx] == DT_TILECACHE_NULL_AREA {
				continue
			}

			sid := 0xff

			// -x
			xidx := (x - 1) + y*w
			if x > 0 && isConnected(layer, idx, xidx, walkableClimb) {
				if layer.Regs[xidx] != 0xff {
					sid = layer.Regs[xidx]
				}

			}

			if sid == 0xff {
				sid = sweepId
				sweepId++
				sweeps[sid].nei = 0xff
				sweeps[sid].ns = 0
			}

			// -y
			yidx := x + (y-1)*w
			if y > 0 && isConnected(layer, idx, yidx, walkableClimb) {
				nr := layer.Regs[yidx]
				if nr != 0xff {
					// Set neighbour when first valid neighbour is encoutered.
					if sweeps[sid].ns == 0 {
						sweeps[sid].nei = nr
					}

					if sweeps[sid].nei == nr {
						// Update existing neighbour
						sweeps[sid].ns++
						prevCount[nr]++
					} else {
						// This is hit if there is nore than one neighbour.
						// Invalidate the neighbour.
						sweeps[sid].nei = 0xff
					}
				}
			}

			layer.Regs[idx] = sid
		}

		// Create unique ID.
		for i := 0; i < sweepId; i++ {
			// If the neighbour is set and there is only one continuous connection to it,
			// the sweep will be merged with the previous one, else new region is created.
			if sweeps[i].nei != 0xff && prevCount[sweeps[i].nei] == sweeps[i].ns {
				sweeps[i].id = sweeps[i].nei
			} else {
				if regId == 255 {
					// Region ID's overflow.
					return DT_FAILURE | DT_BUFFER_TOO_SMALL
				}
				sweeps[i].id = regId
				regId++
			}
		}

		// Remap local sweep ids to region ids.
		for x := 0; x < w; x++ {
			idx := x + y*w
			if layer.Regs[idx] != 0xff {
				layer.Regs[idx] = sweeps[layer.Regs[idx]].id
			}

		}
	}

	// Allocate and init layer regions.
	nregs := regId
	regs := make([]*dtLayerMonotoneRegion, nregs)
	for i := range regs {
		regs[i] = &dtLayerMonotoneRegion{}
	}
	for i := 0; i < nregs; i++ {
		regs[i].regId = 0xff
	}

	// Find region neighbours.
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			idx := x + y*w
			ri := layer.Regs[idx]
			if ri == 0xff {
				continue
			}

			// Update area.
			regs[ri].area++
			regs[ri].areaId = layer.Areas[idx]

			// Update neighbours
			ymi := x + (y-1)*w
			if y > 0 && isConnected(layer, idx, ymi, walkableClimb) {
				rai := layer.Regs[ymi]
				if rai != 0xff && rai != ri {
					addUniqueLast(regs[ri].neis[:], &regs[ri].nneis, rai)
					addUniqueLast(regs[rai].neis[:], &regs[rai].nneis, ri)
				}
			}
		}
	}

	for i := 0; i < nregs; i++ {
		regs[i].regId = i
	}

	for i := 0; i < nregs; i++ {
		reg := regs[i]

		merge := -1
		mergea := 0
		for j := 0; j < reg.nneis; j++ {
			nei := reg.neis[j]
			regn := regs[nei]
			if reg.regId == regn.regId {
				continue
			}

			if reg.areaId != regn.areaId {
				continue
			}

			if regn.area > mergea {
				if canMerge(reg.regId, regn.regId, regs, nregs) {
					mergea = regn.area
					merge = nei
				}
			}
		}
		if merge != -1 {
			oldId := reg.regId
			newId := regs[merge].regId
			for j := 0; j < nregs; j++ {
				if regs[j].regId == oldId {
					regs[j].regId = newId
				}
			}

		}
	}

	// Compact ids.
	remap := make([]int, 256)

	// Find number of unique regions.
	regId = 0
	for i := 0; i < nregs; i++ {
		remap[regs[i].regId] = 1
	}

	for i := 0; i < 256; i++ {
		if remap[i] != 0 {
			remap[i] = regId
			regId++
		}
	}

	// Remap ids.
	for i := 0; i < nregs; i++ {
		regs[i].regId = remap[regs[i].regId]
	}

	layer.RegCount = regId

	for i := 0; i < w*h; i++ {
		if layer.Regs[i] != 0xff {
			layer.Regs[i] = regs[layer.Regs[i]].regId
		}

	}

	return DT_SUCCESS
}

func appendVertex(cont *dtTempContour, x, y, z, r int) bool {
	// Try to merge with existing segments.
	if cont.nverts > 1 {
		pa := rcGetVert4(cont.verts, (cont.nverts - 2))
		pb := rcGetVert4(cont.verts, (cont.nverts - 1))
		if pb[3] == r {
			if pa[0] == pb[0] && pb[0] == x {
				// The verts are aligned aling x-axis, update z.
				pb[1] = y
				pb[2] = z
				return true
			} else if pa[2] == pb[2] && pb[2] == z {
				// The verts are aligned aling z-axis, update x.
				pb[0] = x
				pb[1] = y
				return true
			}
		}
	}

	// Add new point.
	if cont.nverts+1 > cont.cverts {
		return false
	}

	v := rcGetVert4(cont.verts, cont.nverts)
	v[0] = x
	v[1] = y
	v[2] = z
	v[3] = r
	cont.nverts++

	return true
}
func getNeighbourReg(layer *DtTileCacheLayer,
	ax, ay, dir int) int {
	w := layer.Header.Width
	ia := ax + ay*w

	con := layer.Cons[ia] & 0xf
	portal := layer.Cons[ia] >> 4
	mask := (1 << dir)

	if (con & mask) == 0 {
		// No connection, return portal or hard edge.
		if portal&mask > 0 {
			return 0xf8 + dir
		}

		return 0xff
	}

	bx := ax + getDirOffsetX(dir)
	by := ay + getDirOffsetY(dir)
	ib := bx + by*w

	return layer.Regs[ib]
}

func titleCacheWalkContour(layer *DtTileCacheLayer, x, y int, cont *dtTempContour) bool {
	w := layer.Header.Width
	h := layer.Header.Height

	cont.nverts = 0

	startX := x
	startY := y
	startDir := -1

	for i := 0; i < 4; i++ {
		dir := (i + 3) & 3
		rn := getNeighbourReg(layer, x, y, dir)
		if rn != layer.Regs[x+y*w] {
			startDir = dir
			break
		}
	}
	if startDir == -1 {
		return true
	}

	dir := startDir
	maxIter := w * h

	iter := 0
	for iter < maxIter {
		rn := getNeighbourReg(layer, x, y, dir)

		nx := x
		ny := y
		ndir := dir

		if rn != layer.Regs[x+y*w] {
			// Solid edge.
			px := x
			pz := y
			switch dir {
			case 0:
				pz++
				break
			case 1:
				px++
				pz++
				break
			case 2:
				px++
				break
			}

			// Try to merge with previous vertex.
			if !appendVertex(cont, px, layer.Heights[x+y*w], pz, rn) {
				return false
			}

			ndir = (dir + 1) & 0x3 // Rotate CW
		} else {
			// Move to next.
			nx = x + getDirOffsetX(dir)
			ny = y + getDirOffsetY(dir)
			ndir = (dir + 3) & 0x3 // Rotate CCW
		}

		if iter > 0 && x == startX && y == startY && dir == startDir {
			break
		}

		x = nx
		y = ny
		dir = ndir

		iter++
	}

	// Remove last vertex if it is duplicate of the first one.
	pa := rcGetVert4(cont.verts, (cont.nverts - 1))
	pb := cont.verts[0:]
	if pa[0] == pb[0] && pa[2] == pb[2] {
		cont.nverts--
	}

	return true
}

func titleCacheDistancePtSeg(x, z, px, pz, qx, qz int) float64 {
	pqx := float64(qx - px)
	pqz := float64(qz - pz)
	dx := float64(x - px)
	dz := float64(z - pz)
	d := pqx*pqx + pqz*pqz
	t := pqx*dx + pqz*dz
	if d > 0 {
		t /= d
	}

	if t < 0 {
		t = 0
	} else if t > 1 {
		t = 1
	}

	dx = float64(px) + t*pqx - float64(x)
	dz = float64(pz) + t*pqz - float64(z)

	return dx*dx + dz*dz
}

func titleCacheSimplifyContour(cont *dtTempContour, maxError float64) {
	cont.npoly = 0

	for i := 0; i < cont.nverts; i++ {
		j := (i + 1) % cont.nverts
		// Check for start of a wall segment.
		ra := cont.verts[j*4+3]
		rb := cont.verts[i*4+3]
		if ra != rb {
			cont.poly[cont.npoly] = i
			cont.npoly++
		}

	}
	if cont.npoly < 2 {
		// If there is no transitions at all,
		// create some initial points for the simplification process.
		// Find lower-left and upper-right vertices of the contour.
		llx := cont.verts[0]
		llz := cont.verts[2]
		lli := 0
		urx := cont.verts[0]
		urz := cont.verts[2]
		uri := 0
		for i := 1; i < cont.nverts; i++ {
			x := cont.verts[i*4+0]
			z := cont.verts[i*4+2]
			if x < llx || (x == llx && z < llz) {
				llx = x
				llz = z
				lli = i
			}
			if x > urx || (x == urx && z > urz) {
				urx = x
				urz = z
				uri = i
			}
		}
		cont.npoly = 0
		cont.poly[cont.npoly] = lli
		cont.npoly++
		cont.poly[cont.npoly] = uri
		cont.npoly++
	}

	// Add points until all raw points are within
	// error tolerance to the simplified shape.
	for i := 0; i < cont.npoly; {
		ii := (i + 1) % cont.npoly

		ai := cont.poly[i]
		ax := cont.verts[ai*4+0]
		az := cont.verts[ai*4+2]

		bi := cont.poly[ii]
		bx := cont.verts[bi*4+0]
		bz := cont.verts[bi*4+2]

		// Find maximum deviation from the segment.
		maxd := float64(0)
		maxi := -1
		var ci, cinc, endi int

		// Traverse the segment in lexilogical order so that the
		// max deviation is calculated similarly when traversing
		// opposite segments.
		if bx > ax || (bx == ax && bz > az) {
			cinc = 1
			ci = (ai + cinc) % cont.nverts
			endi = bi
		} else {
			cinc = cont.nverts - 1
			ci = (bi + cinc) % cont.nverts
			endi = ai
		}

		// Tessellate only outer edges or edges between areas.
		for ci != endi {
			d := titleCacheDistancePtSeg(cont.verts[ci*4+0], cont.verts[ci*4+2], ax, az, bx, bz)
			if d > maxd {
				maxd = d
				maxi = ci
			}
			ci = (ci + cinc) % cont.nverts
		}

		// If the max deviation is larger than accepted error,
		// add new point, else continue to next segment.
		if maxi != -1 && maxd > (maxError*maxError) {
			cont.npoly++
			for j := cont.npoly - 1; j > i; j-- {
				cont.poly[j] = cont.poly[j-1]
			}

			cont.poly[i+1] = maxi
		} else {
			i++
		}
	}

	// Remap vertices
	start := 0
	for i := 1; i < cont.npoly; i++ {
		if cont.poly[i] < cont.poly[start] {
			start = i
		}
	}

	cont.nverts = 0
	for i := 0; i < cont.npoly; i++ {
		j := (start + i) % cont.npoly
		src := rcGetVert4(cont.verts, cont.poly[j])
		dst := rcGetVert4(cont.verts, cont.nverts)
		dst[0] = src[0]
		dst[1] = src[1]
		dst[2] = src[2]
		dst[3] = src[3]
		cont.nverts++
	}
}

func titleCacheGetCornerHeight(layer *DtTileCacheLayer,
	x, y, z,
	walkableClimb int,
	shouldRemove *bool) int {
	w := layer.Header.Width
	h := layer.Header.Height

	n := 0

	portal := 0xf
	height := 0
	preg := 0xff
	allSameReg := true

	for dz := -1; dz <= 0; dz++ {
		for dx := -1; dx <= 0; dx++ {
			px := x + dx
			pz := z + dz
			if px >= 0 && pz >= 0 && px < w && pz < h {
				idx := px + pz*w
				lh := layer.Heights[idx]
				if common.Abs(lh-y) <= walkableClimb && layer.Areas[idx] != DT_TILECACHE_NULL_AREA {
					height = common.Max(height, lh)
					portal &= (layer.Cons[idx] >> 4)
					if preg != 0xff && preg != layer.Regs[idx] {
						allSameReg = false
					}

					preg = layer.Regs[idx]
					n++
				}
			}
		}
	}

	portalCount := 0
	for dir := 0; dir < 4; dir++ {
		if portal&(1<<dir) > 0 {
			portalCount++
		}
	}

	*shouldRemove = false
	if n > 1 && portalCount == 1 && allSameReg {
		*shouldRemove = true
	}

	return height
}
func titleCacheBuildMeshAdjacency(
	polys []int, npolys int,
	verts []int, nverts int,
	lcset *DtTileCacheContourSet) bool {
	// Based on code by Eric Lengyel from:
	// https://web.archive.org/web/20080704083314/http://www.terathon.com/code/edges.php

	maxEdgeCount := npolys * MAX_VERTS_PER_POLY
	firstEdge := make([]int, nverts+maxEdgeCount)

	nextEdge := firstEdge[nverts:]
	edgeCount := 0

	edges := make([]*rcEdge, maxEdgeCount)
	for i := 0; i < nverts; i++ {
		firstEdge[i] = DT_TILECACHE_NULL_IDX
	}

	for i := 0; i < npolys; i++ {
		t := polys[i*MAX_VERTS_PER_POLY*2:]
		for j := 0; j < MAX_VERTS_PER_POLY; j++ {
			if t[j] == DT_TILECACHE_NULL_IDX {
				break
			}
			v0 := t[j]
			v1 := t[j+1]
			if j+1 >= MAX_VERTS_PER_POLY || t[j+1] == DT_TILECACHE_NULL_IDX {
				v1 = t[0]
			}
			if v0 < v1 {
				edge := edges[edgeCount]
				edge.vert[0] = v0
				edge.vert[1] = v1
				edge.poly[0] = i
				edge.polyEdge[0] = j
				edge.poly[1] = i
				edge.polyEdge[1] = 0xff
				// Insert edge
				nextEdge[edgeCount] = firstEdge[v0]
				firstEdge[v0] = edgeCount
				edgeCount++
			}
		}
	}

	for i := 0; i < npolys; i++ {
		t := polys[i*MAX_VERTS_PER_POLY*2:]
		for j := 0; j < MAX_VERTS_PER_POLY; j++ {
			if t[j] == DT_TILECACHE_NULL_IDX {
				break
			}
			v0 := t[j]
			v1 := t[j+1]
			if j+1 >= MAX_VERTS_PER_POLY || t[j+1] == DT_TILECACHE_NULL_IDX {
				v1 = t[0]
			}
			if v0 > v1 {
				found := false
				for e := firstEdge[v1]; e != DT_TILECACHE_NULL_IDX; e = nextEdge[e] {
					edge := edges[e]
					if edge.vert[1] == v0 && edge.poly[0] == edge.poly[1] {
						edge.poly[1] = i
						edge.polyEdge[1] = j
						found = true
						break
					}
				}
				if !found {
					// Matching edge not found, it is an open edge, add it.
					edge := edges[edgeCount]
					edge.vert[0] = v1
					edge.vert[1] = v0
					edge.poly[0] = i
					edge.polyEdge[0] = j
					edge.poly[1] = i
					edge.polyEdge[1] = 0xff
					// Insert edge
					nextEdge[edgeCount] = firstEdge[v1]
					firstEdge[v1] = edgeCount
					edgeCount++
				}
			}
		}
	}

	// Mark portal edges.
	for i := 0; i < lcset.Nconts; i++ {
		cont := lcset.Conts[i]
		if cont.Nverts < 3 {
			continue
		}

		j := 0
		k := cont.Nverts - 1
		for j < cont.Nverts {
			va := rcGetVert4(cont.Verts, k)
			vb := rcGetVert4(cont.Verts, j)
			dir := va[3] & 0xf
			if dir == 0xf {
				k = j
				j++
				continue
			}
			if dir == 0 || dir == 2 {
				// Find matching vertical edge
				x := va[0]
				zmin := va[2]
				zmax := vb[2]
				if zmin > zmax {
					zmin, zmax = zmax, zmin
				}

				for m := 0; m < edgeCount; m++ {
					e := edges[m]
					// Skip connected edges.
					if e.poly[0] != e.poly[1] {
						k = j
						j++
						continue
					}

					eva := rcGetVert(verts, e.vert[0])
					evb := rcGetVert(verts, e.vert[1])
					if eva[0] == x && evb[0] == x {
						ezmin := eva[2]
						ezmax := evb[2]
						if ezmin > ezmax {
							ezmin, ezmax = ezmax, ezmin
						}

						if overlapRangeExl(zmin, zmax, ezmin, ezmax) {
							// Reuse the other polyedge to store dir.
							e.polyEdge[1] = dir
						}
					}
				}
			} else {
				// Find matching vertical edge
				z := va[2]
				xmin := va[0]
				xmax := vb[0]
				if xmin > xmax {
					xmin, xmax = xmax, xmin
				}

				for m := 0; m < edgeCount; m++ {
					e := edges[m]
					// Skip connected edges.
					if e.poly[0] != e.poly[1] {
						k = j
						j++
						continue
					}

					eva := rcGetVert(verts, e.vert[0])
					evb := rcGetVert(verts, e.vert[1])
					if eva[2] == z && evb[2] == z {
						exmin := eva[0]
						exmax := evb[0]
						if exmin > exmax {
							exmin, exmax = exmax, exmin
						}

						if overlapRangeExl(xmin, xmax, exmin, exmax) {
							// Reuse the other polyedge to store dir.
							e.polyEdge[1] = dir
						}
					}
				}
			}
			k = j
			j++
		}
	}

	// Store adjacency
	for i := 0; i < edgeCount; i++ {
		e := edges[i]
		if e.poly[0] != e.poly[1] {
			p0 := polys[e.poly[0]*MAX_VERTS_PER_POLY*2:]
			p1 := polys[e.poly[1]*MAX_VERTS_PER_POLY*2:]
			p0[MAX_VERTS_PER_POLY+e.polyEdge[0]] = e.poly[1]
			p1[MAX_VERTS_PER_POLY+e.polyEdge[1]] = e.poly[0]
		} else if e.polyEdge[1] != 0xff {
			p0 := polys[e.poly[0]*MAX_VERTS_PER_POLY*2:]
			p0[MAX_VERTS_PER_POLY+e.polyEdge[0]] = 0x8000 | e.polyEdge[1]
		}

	}

	return true
}

// TODO: move this somewhere else, once the layer meshing is done.
func dtBuildTileCacheContours(
	layer *DtTileCacheLayer,
	walkableClimb int, maxError float64,
	lcset *DtTileCacheContourSet) DtStatus {

	w := layer.Header.Width
	h := layer.Header.Height

	lcset.Nconts = layer.RegCount
	lcset.Conts = make([]*DtTileCacheContour, lcset.Nconts)
	for i := range lcset.Conts {
		lcset.Conts[i] = &DtTileCacheContour{}
	}

	// Allocate temp buffer for contour tracing.
	maxTempVerts := (w + h) * 2 * 2 // Twice around the layer.
	tempVerts := make([]int, maxTempVerts*4)
	tempPoly := make([]int, maxTempVerts)
	temp := newDtTempContour(tempVerts, maxTempVerts, tempPoly, maxTempVerts)

	// Find contours.
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			idx := x + y*w
			ri := layer.Regs[idx]
			if ri == 0xff {
				continue
			}

			cont := lcset.Conts[ri]

			if cont.Nverts > 0 {
				continue
			}

			cont.Reg = ri
			cont.Area = layer.Areas[idx]

			if !titleCacheWalkContour(layer, x, y, temp) {
				// Too complex contour.
				// Note: If you hit here ofte, try increasing 'maxTempVerts'.
				return DT_FAILURE | DT_BUFFER_TOO_SMALL
			}

			titleCacheSimplifyContour(temp, maxError)

			// Store contour.
			cont.Nverts = temp.nverts
			if cont.Nverts > 0 {
				cont.Verts = make([]int, 4*temp.nverts)
				i := 0
				j := temp.nverts - 1
				for i < temp.nverts {
					dst := rcGetVert4(cont.Verts, j)
					v := rcGetVert4(temp.verts, j)
					vn := rcGetVert4(temp.verts, i)
					nei := vn[3] // The neighbour reg is stored at segment vertex of a segment.
					shouldRemove := false
					lh := titleCacheGetCornerHeight(layer, v[0], v[1], v[2],
						walkableClimb, &shouldRemove)

					dst[0] = v[0]
					dst[1] = lh
					dst[2] = v[2]

					// Store portal direction and remove status to the fourth component.
					dst[3] = 0x0f
					if nei != 0xff && nei >= 0xf8 {
						dst[3] = nei - 0xf8
					}

					if shouldRemove {
						dst[3] |= 0x80
					}

				}
				j = i
				i++

			}
		}
	}

	return DT_SUCCESS
}

const VERTEX_BUCKET_COUNT2 = (1 << 8)

func computeVertexHash2(x, y, z int) int {
	h1 := 0x8da6b343 // Large multiplicative constants;
	h2 := 0xd8163841 // here arbitrarily chosen primes
	h3 := 0xcb1ab31f
	n := h1*x + h2*y + h3*z
	return (n & (VERTEX_BUCKET_COUNT2 - 1))
}

func titleCacheAddVertex(x, y, z int,
	verts []int, firstVert []int, nextVert []int, nv *int) int {
	bucket := computeVertexHash2(x, 0, z)
	i := firstVert[bucket]

	for i != DT_TILECACHE_NULL_IDX {
		v := rcGetVert(verts, i)
		if v[0] == x && v[2] == z && (common.Abs(v[1]-y) <= 2) {
			return i
		}

		i = nextVert[i] // next
	}

	// Could not find, create new.
	i = *nv
	*nv++
	v := rcGetVert(verts, i)
	v[0] = x
	v[1] = y
	v[2] = z
	nextVert[i] = firstVert[bucket]
	firstVert[bucket] = i

	return i
}

type titleCacheRcEdge struct {
	vert     [2]int
	polyEdge [2]int
	poly     [2]int
}

func dtDecompressTileCacheLayer(comp *dtTileCacheCompressor,
	data []byte, compressedSize int,
) (layerOut *DtTileCacheLayer, status DtStatus) {
	//if compressedHeader.magic != DT_TILECACHE_MAGIC {
	//	return nil, DT_FAILURE | DT_WRONG_MAGIC
	//}
	//
	//if compressedHeader.version != DT_TILECACHE_VERSION {
	//	return nil, DT_FAILURE | DT_WRONG_VERSION
	//}

	//const int layerSize = dtAlign4(sizeof(DtTileCacheLayer));
	//const int headerSize = dtAlign4(sizeof(DtTileCacheLayerHeader));
	//gridSize := compressedHeader.width * compressedHeader.height;
	//const int bufferSize = layerSize + headerSize + gridSize*4;
	//
	//buffer = (unsigned char*)alloc->alloc(bufferSize);
	//
	//layer := &DtTileCacheLayer{}
	//header = &DtTileCacheLayerHeader{}
	//unsigned char* grids = buffer + layerSize + headerSize;
	//const int gridsSize = bufferSize - (layerSize + headerSize);
	//
	//// Copy header
	//memcpy(header, compressedHeader, headerSize);
	//// Decompress grid.
	//int size = 0;
	//DtStatus status = comp.decompress(compressed+headerSize, compressedSize-headerSize,
	//	grids, gridsSize, &size);
	//if (status.DtStatusFailed()) {
	//	alloc->free(buffer);
	//	return status;
	//}
	//
	//layer.header = header;
	//layer.heights = grids;
	//layer.areas = grids + gridSize;
	//layer.cons = grids + gridSize*2;
	//layer.regs = grids + gridSize*3;
	return layerOut, DT_SUCCESS

}

func dtMarkCylinderArea(layer *DtTileCacheLayer, orig []float64, cs float64, ch float64,
	pos []float64, radius float64, height float64, areaId int) DtStatus {
	bmin := make([]float64, 3)
	bmax := make([]float64, 3)
	bmin[0] = pos[0] - radius
	bmin[1] = pos[1]
	bmin[2] = pos[2] - radius
	bmax[0] = pos[0] + radius
	bmax[1] = pos[1] + height
	bmax[2] = pos[2] + radius
	r2 := common.Sqr(radius/cs + 0.5)

	w := layer.Header.Width
	h := layer.Header.Height
	ics := 1.0 / cs
	ich := 1.0 / ch

	px := (pos[0] - orig[0]) * ics
	pz := (pos[2] - orig[2]) * ics

	minx := int(math.Floor((bmin[0] - orig[0]) * ics))
	miny := int(math.Floor((bmin[1] - orig[1]) * ich))
	minz := int(math.Floor((bmin[2] - orig[2]) * ics))
	maxx := int(math.Floor((bmax[0] - orig[0]) * ics))
	maxy := int(math.Floor((bmax[1] - orig[1]) * ich))
	maxz := int(math.Floor((bmax[2] - orig[2]) * ics))

	if maxx < 0 {
		return DT_SUCCESS
	}
	if minx >= w {
		return DT_SUCCESS
	}
	if maxz < 0 {
		return DT_SUCCESS
	}
	if minz >= h {
		return DT_SUCCESS
	}

	if minx < 0 {
		minx = 0
	}
	if maxx >= w {
		maxx = w - 1
	}
	if minz < 0 {
		minz = 0
	}
	if maxz >= h {
		maxz = h - 1
	}

	for z := minz; z <= maxz; z++ {
		for x := minx; x <= maxx; x++ {
			dx := (float64(x) + 0.5) - float64(px)
			dz := (float64(z) + 0.5) - float64(pz)
			if dx*dx+dz*dz > r2 {
				continue
			}

			y := layer.Heights[x+z*w]
			if y < miny || y > maxy {
				continue
			}

			layer.Areas[x+z*w] = areaId
		}
	}

	return DT_SUCCESS
}

func dtMarkBoxArea(layer *DtTileCacheLayer, orig []float64, cs, ch float64,
	bmin, bmax []float64, areaId int) DtStatus {
	w := layer.Header.Width
	h := layer.Header.Height
	ics := 1.0 / cs
	ich := 1.0 / ch

	minx := int(math.Floor((bmin[0] - orig[0]) * ics))
	miny := int(math.Floor((bmin[1] - orig[1]) * ich))
	minz := int(math.Floor((bmin[2] - orig[2]) * ics))
	maxx := int(math.Floor((bmax[0] - orig[0]) * ics))
	maxy := int(math.Floor((bmax[1] - orig[1]) * ich))
	maxz := int(math.Floor((bmax[2] - orig[2]) * ics))

	if maxx < 0 {
		return DT_SUCCESS
	}
	if minx >= w {
		return DT_SUCCESS
	}
	if maxz < 0 {
		return DT_SUCCESS
	}
	if minz >= h {
		return DT_SUCCESS
	}

	if minx < 0 {
		minx = 0
	}
	if maxx >= w {
		maxx = w - 1
	}
	if minz < 0 {
		minz = 0
	}
	if maxz >= h {
		maxz = h - 1
	}

	for z := minz; z <= maxz; z++ {
		for x := minx; x <= maxx; x++ {
			y := layer.Heights[x+z*w]
			if y < miny || y > maxy {
				continue
			}

			layer.Areas[x+z*w] = areaId
		}
	}

	return DT_SUCCESS
}

func dtMarkBoxArea1(layer *DtTileCacheLayer, orig []float64, cs float64, ch float64,
	center []float64, halfExtents []float64, rotAux []float64, areaId int) DtStatus {
	w := layer.Header.Width
	h := layer.Header.Height
	ics := 1.0 / cs
	ich := 1.0 / ch

	cx := (center[0] - orig[0]) * ics
	cz := (center[2] - orig[2]) * ics

	maxr := 1.41 * common.Max(halfExtents[0], halfExtents[2])
	minx := int(math.Floor(cx - maxr*ics))
	maxx := int(math.Floor(cx + maxr*ics))
	minz := int(math.Floor(cz - maxr*ics))
	maxz := int(math.Floor(cz + maxr*ics))
	miny := int(math.Floor((center[1] - halfExtents[1] - orig[1]) * ich))
	maxy := int(math.Floor((center[1] + halfExtents[1] - orig[1]) * ich))

	if maxx < 0 {
		return DT_SUCCESS
	}
	if minx >= w {
		return DT_SUCCESS
	}
	if maxz < 0 {
		return DT_SUCCESS
	}
	if minz >= h {
		return DT_SUCCESS
	}

	if minx < 0 {
		minx = 0
	}
	if maxx >= w {
		maxx = w - 1
	}
	if minz < 0 {
		minz = 0
	}
	if maxz >= h {
		maxz = h - 1
	}

	xhalf := halfExtents[0]*ics + 0.5
	zhalf := halfExtents[2]*ics + 0.5

	for z := minz; z <= maxz; z++ {
		for x := minx; x <= maxx; x++ {
			x2 := 2.0 * (float64(x) - cx)
			z2 := 2.0 * (float64(z) - cz)
			xrot := rotAux[1]*x2 + rotAux[0]*z2
			if xrot > xhalf || xrot < -xhalf {
				continue
			}

			zrot := rotAux[1]*z2 - rotAux[0]*x2
			if zrot > zhalf || zrot < -zhalf {
				continue
			}

			y := layer.Heights[x+z*w]
			if y < miny || y > maxy {
				continue
			}

			layer.Areas[x+z*w] = areaId
		}
	}

	return DT_SUCCESS
}

func dtBuildTileCachePolyMesh(
	lcset *DtTileCacheContourSet,
	mesh *DtTileCachePolyMesh) DtStatus {
	maxVertices := 0
	maxTris := 0
	maxVertsPerCont := 0
	for i := 0; i < lcset.Nconts; i++ {
		// Skip null contours.
		if lcset.Conts[i].Nverts < 3 {
			continue
		}
		maxVertices += lcset.Conts[i].Nverts
		maxTris += lcset.Conts[i].Nverts - 2
		maxVertsPerCont = common.Max(maxVertsPerCont, lcset.Conts[i].Nverts)
	}

	// TODO: warn about too many vertices?

	mesh.Nvp = MAX_VERTS_PER_POLY

	vflags := make([]int, maxVertices)
	mesh.Verts = make([]int, maxVertices*3)

	mesh.Polys = make([]int, maxTris*MAX_VERTS_PER_POLY*2)
	mesh.Areas = make([]int, maxTris)
	mesh.Flags = make([]int, maxTris)
	// Just allocate and clean the mesh flags array. The user is resposible for filling it.

	mesh.Nverts = 0
	mesh.Npolys = 0
	for i := range mesh.Polys {
		mesh.Polys[i] = 0xff
	}
	firstVert := make([]int, VERTEX_BUCKET_COUNT2)
	for i := 0; i < VERTEX_BUCKET_COUNT2; i++ {
		firstVert[i] = DT_TILECACHE_NULL_IDX
	}

	nextVert := make([]int, maxVertices)

	indices := make([]int, maxVertsPerCont)
	tris := make([]int, maxVertsPerCont*3)

	polys := make([]int, maxVertsPerCont*MAX_VERTS_PER_POLY)

	for i := 0; i < lcset.Nconts; i++ {
		cont := lcset.Conts[i]

		// Skip null contours.
		if cont.Nverts < 3 {
			continue
		}

		// Triangulate contour
		for j := 0; j < cont.Nverts; j++ {
			indices[j] = j
		}

		ntris := titleCacheTriangulate(cont.Nverts, cont.Verts, indices, tris)
		if ntris <= 0 {
			// TODO: issue warning!
			ntris = -ntris
		}

		// Add and merge vertices.
		for j := 0; j < cont.Nverts; j++ {
			v := rcGetVert4(cont.Verts, j)
			indices[j] = titleCacheAddVertex(v[0], v[1], v[2],
				mesh.Verts, firstVert, nextVert, &mesh.Nverts)
			if v[3]&0x80 > 0 {
				// This vertex should be removed.
				vflags[indices[j]] = 1
			}
		}

		// Build initial polygons.
		npolys := 0
		for i := 0; i < maxVertsPerCont*MAX_VERTS_PER_POLY; i++ {
			polys[i] = 0xff
		}
		for j := 0; j < ntris; j++ {
			t := rcGetVert(tris, j)
			if t[0] != t[1] && t[0] != t[2] && t[1] != t[2] {
				polys[npolys*MAX_VERTS_PER_POLY+0] = indices[t[0]]
				polys[npolys*MAX_VERTS_PER_POLY+1] = indices[t[1]]
				polys[npolys*MAX_VERTS_PER_POLY+2] = indices[t[2]]
				npolys++
			}
		}
		if npolys == 0 {
			continue
		}

		// Merge polygons.
		maxVertsPerPoly := MAX_VERTS_PER_POLY
		if maxVertsPerPoly > 3 {
			for {
				// Find best polygons to merge.
				bestMergeVal := 0
				var bestPa = 0
				var bestPb = 0
				var bestEa = 0
				var bestEb = 0

				for j := 0; j < npolys-1; j++ {
					pj := polys[j*MAX_VERTS_PER_POLY:]
					for k := j + 1; k < npolys; k++ {
						pk := polys[k*MAX_VERTS_PER_POLY:]
						var ea, eb int
						v := titleCacheGetPolyMergeValue(pj, pk, mesh.Verts, &ea, &eb)
						if v > bestMergeVal {
							bestMergeVal = v
							bestPa = j
							bestPb = k
							bestEa = ea
							bestEb = eb
						}
					}
				}

				if bestMergeVal > 0 {
					// Found best, merge.
					pa := polys[bestPa*MAX_VERTS_PER_POLY:]
					pb := polys[bestPb*MAX_VERTS_PER_POLY:]
					mergePolys(pa, pb, bestEa, bestEb)
					copy(pb, polys[(npolys-1)*MAX_VERTS_PER_POLY:(npolys-1)*MAX_VERTS_PER_POLY+MAX_VERTS_PER_POLY])
					npolys--
				} else {
					// Could not merge any polygons, stop.
					break
				}
			}
		}

		// Store polygons.
		for j := 0; j < npolys; j++ {
			p := mesh.Polys[mesh.Npolys*MAX_VERTS_PER_POLY*2:]
			q := polys[j*MAX_VERTS_PER_POLY:]
			for k := 0; k < MAX_VERTS_PER_POLY; k++ {
				p[k] = q[k]
			}

			mesh.Areas[mesh.Npolys] = cont.Area
			mesh.Npolys++
			if mesh.Npolys > maxTris {
				return DT_FAILURE | DT_BUFFER_TOO_SMALL
			}

		}
	}

	// Remove edge vertices.
	for i := 0; i < mesh.Nverts; i++ {
		if vflags[i] != 0 {
			if !titleCacheCanRemoveVertex(mesh, i) {
				continue
			}

			status := titleCacheRemoveVertex(mesh, i, maxTris)
			if status.DtStatusFailed() {
				return status
			}

			// Remove vertex
			// Note: mesh.nverts is already decremented inside removeVertex()!
			for j := i; j < mesh.Nverts; j++ {
				vflags[j] = vflags[j+1]
			}

			i--
		}
	}

	// Calculate adjacency.
	if !titleCacheBuildMeshAdjacency(mesh.Polys, mesh.Npolys, mesh.Verts, mesh.Nverts, lcset) {
		return DT_FAILURE | DT_OUT_OF_MEMORY
	}

	return DT_SUCCESS
}

func titleCacheTriangulate(n int, verts []int, indices []int, tris []int) int {
	ntris := 0
	dst := 0

	// The last bit of the index is used to indicate if the vertex can be removed.
	for i := 0; i < n; i++ {
		i1 := next(i, n)
		i2 := next(i1, n)
		if diagonal(i, i2, n, verts, indices) {
			indices[i1] |= 0x8000
		}

	}

	for n > 3 {
		minLen := -1
		mini := -1
		for i := 0; i < n; i++ {
			i1 := next(i, n)
			if indices[i1]&0x8000 > 0 {
				p0 := rcGetVert4(verts, (indices[i] & 0x7fff))
				p2 := rcGetVert4(verts, (indices[next(i1, n)] & 0x7fff))

				dx := p2[0] - p0[0]
				dz := p2[2] - p0[2]
				length := dx*dx + dz*dz
				if minLen < 0 || length < minLen {
					minLen = length
					mini = i
				}
			}
		}

		if mini == -1 {
			// Should not happen.
			/*			printf("mini == -1 ntris=%d n=%d\n", ntris, n);
						for (int i = 0; i < n; i++)
						{
						printf("%d ", indices[i] & 0x0fffffff);
						}
						printf("\n");*/
			return -ntris
		}

		i := mini
		i1 := next(i, n)
		i2 := next(i1, n)

		tris[dst] = indices[i] & 0x7fff
		dst++
		tris[dst] = indices[i1] & 0x7fff
		dst++
		tris[dst] = indices[i2] & 0x7fff
		dst++
		ntris++

		// Removes P[i1] by copying P[i+1]...P[n-1] left one index.
		n--
		for k := i1; k < n; k++ {
			indices[k] = indices[k+1]
		}

		if i1 >= n {
			i1 = 0
		}
		i = prev(i1, n)
		// Update diagonal flags.
		if diagonal(prev(i, n), i1, n, verts, indices) {
			indices[i] |= 0x8000
		} else {
			indices[i] &= 0x7fff
		}

		if diagonal(i, next(i1, n), n, verts, indices) {
			indices[i1] |= 0x8000
		} else {
			indices[i1] &= 0x7fff
		}
	}

	// Append the remaining triangle.
	tris[dst] = indices[0] & 0x7fff
	dst++
	tris[dst] = indices[1] & 0x7fff
	dst++
	tris[dst] = indices[2] & 0x7fff
	dst++
	ntris++

	return ntris
}

func titleCacheCountPolyVerts(p []int) int {
	for i := 0; i < MAX_VERTS_PER_POLY; i++ {
		if p[i] == DT_TILECACHE_NULL_IDX {
			return i
		}
	}

	return MAX_VERTS_PER_POLY
}

func mergePolys(pa, pb []int, ea, eb int) {
	tmp := make([]int, MAX_VERTS_PER_POLY*2)

	na := titleCacheCountPolyVerts(pa)
	nb := titleCacheCountPolyVerts(pb)

	// Merge polygons.
	for i := range tmp {
		tmp[i] = 0xff
	}
	n := 0
	// Add pa
	for i := 0; i < na-1; i++ {
		tmp[n] = pa[(ea+1+i)%na]
		n++
	}

	// Add pb
	for i := 0; i < nb-1; i++ {
		tmp[n] = pb[(eb+1+i)%nb]
		n++
	}

	copy(pa, tmp[:MAX_VERTS_PER_POLY])
}

func titleCacheGetPolyMergeValue(pa, pb []int,
	verts []int, ea, eb *int) int {
	na := titleCacheCountPolyVerts(pa)
	nb := titleCacheCountPolyVerts(pb)

	// If the merged polygon would be too big, do not merge.
	if na+nb-2 > MAX_VERTS_PER_POLY {
		return -1
	}

	// Check if the polygons share an edge.
	*ea = -1
	*eb = -1

	for i := 0; i < na; i++ {
		va0 := pa[i]
		va1 := pa[(i+1)%na]
		if va0 > va1 {
			va0, va1 = va1, va0
		}

		for j := 0; j < nb; j++ {
			vb0 := pb[j]
			vb1 := pb[(j+1)%nb]
			if vb0 > vb1 {
				vb0, vb1 = vb1, vb0
			}

			if va0 == vb0 && va1 == vb1 {
				*ea = i
				*eb = j
				break
			}
		}
	}

	// No common edge, cannot merge.
	if *ea == -1 || *eb == -1 {
		return -1
	}

	// Check to see if the merged polygon would be convex.
	var va, vb, vc int

	va = pa[(*ea+na-1)%na]
	vb = pa[*ea]
	vc = pb[(*eb+2)%nb]
	if !uleft(rcGetVert(verts, va), rcGetVert(verts, vb), rcGetVert(verts, vc)) {
		return -1
	}

	va = pb[(*eb+nb-1)%nb]
	vb = pb[*eb]
	vc = pa[(*ea+2)%na]
	if !uleft(rcGetVert(verts, va), rcGetVert(verts, vb), rcGetVert(verts, vc)) {
		return -1
	}

	va = pa[*ea]
	vb = pa[(*ea+1)%na]

	dx := verts[va*3+0] - verts[vb*3+0]
	dy := verts[va*3+2] - verts[vb*3+2]

	return dx*dx + dy*dy
}

func titleCacheCanRemoveVertex(mesh *DtTileCachePolyMesh, rem int) bool {
	// Count number of polygons to remove.
	numTouchedVerts := 0
	numRemainingEdges := 0
	for i := 0; i < mesh.Npolys; i++ {
		p := mesh.Polys[i*MAX_VERTS_PER_POLY*2:]
		nv := titleCacheCountPolyVerts(p)
		numRemoved := 0
		numVerts := 0
		for j := 0; j < nv; j++ {
			if p[j] == rem {
				numTouchedVerts++
				numRemoved++
			}
			numVerts++
		}
		if numRemoved != 0 {
			numRemainingEdges += numVerts - (numRemoved + 1)
		}
	}

	// There would be too few edges remaining to create a polygon.
	// This can happen for example when a tip of a triangle is marked
	// as deletion, but there are no other polys that share the vertex.
	// In this case, the vertex should not be removed.
	if numRemainingEdges <= 2 {
		return false
	}

	// Check that there is enough memory for the test.
	maxEdges := numTouchedVerts * 2
	if maxEdges > MAX_REM_EDGES {
		return false
	}

	// Find edges which share the removed vertex.
	edges := make([]int, MAX_REM_EDGES)
	nedges := 0

	for i := 0; i < mesh.Npolys; i++ {
		p := mesh.Polys[i*MAX_VERTS_PER_POLY*2:]
		nv := titleCacheCountPolyVerts(p)
		j := 0
		k := nv - 1
		// Collect edges which touches the removed vertex.
		for j < nv {
			if p[j] == rem || p[k] == rem {
				// Arrange edge so that a=rem.
				a := p[j]
				b := p[k]
				if b == rem {
					a, b = b, a
				}

				// Check if the edge exists
				exists := false
				for m := 0; m < nedges; m++ {
					e := edges[m*3:]
					if e[1] == b {
						// Exists, increment vertex share count.
						e[2]++
						exists = true
					}
				}
				// Add new edge.
				if !exists {
					e := edges[nedges*3:]
					e[0] = a
					e[1] = b
					e[2] = 1
					nedges++
				}
			}
			k = j
			j++
		}
	}

	// There should be no more than 2 open edges.
	// This catches the case that two non-adjacent polygons
	// share the removed vertex. In that case, do not remove the vertex.
	numOpenEdges := 0
	for i := 0; i < nedges; i++ {
		if edges[i*3+2] < 2 {
			numOpenEdges++
		}

	}
	if numOpenEdges > 2 {
		return false
	}

	return true
}

func titleCacheRemoveVertex(mesh *DtTileCachePolyMesh, rem int, maxTris int) DtStatus {
	// Count number of polygons to remove.
	numRemovedVerts := 0
	for i := 0; i < mesh.Npolys; i++ {
		p := mesh.Polys[i*MAX_VERTS_PER_POLY*2:]
		nv := titleCacheCountPolyVerts(p)
		for j := 0; j < nv; j++ {
			if p[j] == rem {
				numRemovedVerts++
			}

		}
	}

	nedges := 0
	edges := make([]int, MAX_REM_EDGES*3)
	nhole := 0
	hole := make([]int, MAX_REM_EDGES)
	nharea := 0
	harea := make([]int, MAX_REM_EDGES)

	for i := 0; i < mesh.Npolys; i++ {
		p := mesh.Polys[i*MAX_VERTS_PER_POLY*2:]
		nv := titleCacheCountPolyVerts(p)
		hasRem := false
		for j := 0; j < nv; j++ {
			if p[j] == rem {
				hasRem = true
			}
		}

		if hasRem {
			// Collect edges which does not touch the removed vertex.
			j := 0
			k := nv - 1
			for j < nv {
				if p[j] != rem && p[k] != rem {
					if nedges >= MAX_REM_EDGES {
						return DT_FAILURE | DT_BUFFER_TOO_SMALL
					}

					e := edges[nedges*3:]
					e[0] = p[k]
					e[1] = p[j]
					e[2] = mesh.Areas[i]
					nedges++
				}
				k = j
				j++
			}
			// Remove the polygon.
			p2 := mesh.Polys[(mesh.Npolys-1)*MAX_VERTS_PER_POLY*2:]
			copy(p, p2[:MAX_VERTS_PER_POLY])
			for i := MAX_VERTS_PER_POLY; i < MAX_VERTS_PER_POLY+MAX_VERTS_PER_POLY; i++ {
				p[i] = 0xff
			}
			mesh.Areas[i] = mesh.Areas[mesh.Npolys-1]
			mesh.Npolys--
			i--
		}
	}

	// Remove vertex.
	for i := rem; i < mesh.Nverts-1; i++ {
		mesh.Verts[i*3+0] = mesh.Verts[(i+1)*3+0]
		mesh.Verts[i*3+1] = mesh.Verts[(i+1)*3+1]
		mesh.Verts[i*3+2] = mesh.Verts[(i+1)*3+2]
	}
	mesh.Nverts--

	// Adjust indices to match the removed vertex layout.
	for i := 0; i < mesh.Npolys; i++ {
		p := mesh.Polys[i*MAX_VERTS_PER_POLY*2:]
		nv := titleCacheCountPolyVerts(p)
		for j := 0; j < nv; j++ {
			if p[j] > rem {
				p[j]--
			}
		}

	}
	for i := 0; i < nedges; i++ {
		if edges[i*3+0] > rem {
			edges[i*3+0]--
		}
		if edges[i*3+1] > rem {
			edges[i*3+1]--
		}
	}

	if nedges == 0 {
		return DT_SUCCESS
	}

	// Start with one vertex, keep appending connected
	// segments to the start and end of the hole.
	pushBack(edges[0], hole, &nhole)
	pushBack(edges[2], harea, &nharea)

	for nedges > 0 {
		match := false
		for i := 0; i < nedges; i++ {
			ea := edges[i*3+0]
			eb := edges[i*3+1]
			a := edges[i*3+2]
			add := false
			if hole[0] == eb {
				// The segment matches the beginning of the hole boundary.
				if nhole >= MAX_REM_EDGES {
					return DT_FAILURE | DT_BUFFER_TOO_SMALL
				}

				pushFront(ea, hole, &nhole)
				pushFront(a, harea, &nharea)
				add = true
			} else if hole[nhole-1] == ea {
				// The segment matches the end of the hole boundary.
				if nhole >= MAX_REM_EDGES {
					return DT_FAILURE | DT_BUFFER_TOO_SMALL
				}

				pushBack(eb, hole, &nhole)
				pushBack(a, harea, &nharea)
				add = true
			}
			if add {
				// The edge segment was added, remove it.
				edges[i*3+0] = edges[(nedges-1)*3+0]
				edges[i*3+1] = edges[(nedges-1)*3+1]
				edges[i*3+2] = edges[(nedges-1)*3+2]
				nedges--
				match = true
				i--
			}
		}

		if !match {
			break
		}

	}

	tris := make([]int, MAX_REM_EDGES*3)
	tverts := make([]int, MAX_REM_EDGES*3)
	tpoly := make([]int, MAX_REM_EDGES*3)

	// Generate temp vertex array for triangulation.
	for i := 0; i < nhole; i++ {
		pi := hole[i]
		tverts[i*4+0] = mesh.Verts[pi*3+0]
		tverts[i*4+1] = mesh.Verts[pi*3+1]
		tverts[i*4+2] = mesh.Verts[pi*3+2]
		tverts[i*4+3] = 0
		tpoly[i] = i
	}

	// Triangulate the hole.
	ntris := triangulate(nhole, tverts, tpoly, tris)
	if ntris < 0 {
		// TODO: issue warning!
		ntris = -ntris
	}

	if ntris > MAX_REM_EDGES {
		return DT_FAILURE | DT_BUFFER_TOO_SMALL
	}

	polys := make([]int, MAX_REM_EDGES*MAX_VERTS_PER_POLY)
	pareas := make([]int, MAX_REM_EDGES)

	// Build initial polygons.
	npolys := 0
	for i := 0; i < ntris*MAX_VERTS_PER_POLY; i++ {
		polys[i] = 0xff
	}
	for j := 0; j < ntris; j++ {
		t := tris[j*3:]
		if t[0] != t[1] && t[0] != t[2] && t[1] != t[2] {
			polys[npolys*MAX_VERTS_PER_POLY+0] = hole[t[0]]
			polys[npolys*MAX_VERTS_PER_POLY+1] = hole[t[1]]
			polys[npolys*MAX_VERTS_PER_POLY+2] = hole[t[2]]
			pareas[npolys] = harea[t[0]]
			npolys++
		}
	}
	if npolys == 0 {
		return DT_SUCCESS
	}

	// Merge polygons.
	maxVertsPerPoly := MAX_VERTS_PER_POLY
	if maxVertsPerPoly > 3 {
		for {
			// Find best polygons to merge.
			bestMergeVal := 0
			var bestPa = 0
			var bestPb = 0
			var bestEa = 0
			var bestEb = 0

			for j := 0; j < npolys-1; j++ {
				pj := polys[j*MAX_VERTS_PER_POLY:]
				for k := j + 1; k < npolys; k++ {
					pk := polys[k*MAX_VERTS_PER_POLY:]
					var ea, eb int
					v := titleCacheGetPolyMergeValue(pj, pk, mesh.Verts, &ea, &eb)
					if v > bestMergeVal {
						bestMergeVal = v
						bestPa = j
						bestPb = k
						bestEa = ea
						bestEb = eb
					}
				}
			}

			if bestMergeVal > 0 {
				// Found best, merge.
				pa := polys[bestPa*MAX_VERTS_PER_POLY:]
				pb := polys[bestPb*MAX_VERTS_PER_POLY:]
				mergePolys(pa, pb, bestEa, bestEb)
				copy(pb, polys[(npolys-1)*MAX_VERTS_PER_POLY:(npolys-1)*MAX_VERTS_PER_POLY+MAX_VERTS_PER_POLY])
				pareas[bestPb] = pareas[npolys-1]
				npolys--
			} else {
				// Could not merge any polygons, stop.
				break
			}
		}
	}

	// Store polygons.
	for i := 0; i < npolys; i++ {
		if mesh.Npolys >= maxTris {
			break
		}
		p := mesh.Polys[mesh.Npolys*MAX_VERTS_PER_POLY*2:]
		for i := 0; i < MAX_VERTS_PER_POLY*2; i++ {
			p[i] = 0xff
		}
		for j := 0; j < MAX_VERTS_PER_POLY; j++ {
			p[j] = polys[i*MAX_VERTS_PER_POLY+j]
		}

		mesh.Areas[mesh.Npolys] = pareas[i]
		mesh.Npolys++
		if mesh.Npolys > maxTris {
			return DT_FAILURE | DT_BUFFER_TOO_SMALL
		}

	}

	return DT_SUCCESS
}
