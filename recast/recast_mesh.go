package recast

import (
	"github.com/gorustyt/gonavmesh/common"
	"log"
	"math"
	"reflect"
)

const (
	VERTEX_BUCKET_COUNT = (1 << 12)
	RC_MESH_NULL_IDX    = 0xffff
)

type RcEdge struct {
	Vert     [2]uint16
	PolyEdge [2]uint16
	Poly     [2]uint16
}

func buildMeshAdjacency(polys []uint16, npolys, nverts, vertsPerPoly int32) bool {
	// Based on code by Eric Lengyel from:
	// https://web.archive.org/web/20080704083314/http://www.terathon.com/code/edges.php

	maxEdgeCount := npolys * vertsPerPoly
	firstEdge := make([]uint16, nverts+maxEdgeCount)
	nextEdge := firstEdge[nverts:]
	edgeCount := int32(0)
	edges := make([]*RcEdge, maxEdgeCount)
	for i := range edges {
		edges[i] = &RcEdge{}
	}
	for i := int32(0); i < nverts; i++ {
		firstEdge[i] = RC_MESH_NULL_IDX
	}

	for i := int32(0); i < npolys; i++ {
		t := common.GetVert2(polys, i*vertsPerPoly)
		for j := int32(0); j < vertsPerPoly; j++ {
			if t[j] == RC_MESH_NULL_IDX {
				break
			}
			v0 := t[j]
			v1 := t[j+1]
			if j+1 >= vertsPerPoly || t[j+1] == RC_MESH_NULL_IDX {
				v1 = t[0]
			}
			if v0 < v1 {
				edge := edges[edgeCount]
				edge.Vert[0] = v0
				edge.Vert[1] = v1
				edge.Poly[0] = uint16(i)
				edge.PolyEdge[0] = uint16(j)
				edge.Poly[1] = uint16(i)
				edge.PolyEdge[1] = 0
				// Insert edge
				nextEdge[edgeCount] = firstEdge[v0]
				firstEdge[v0] = uint16(edgeCount)
				edgeCount++
			}
		}
	}

	for i := int32(0); i < npolys; i++ {
		t := common.GetVert2(polys, i*vertsPerPoly)
		for j := int32(0); j < vertsPerPoly; j++ {
			if t[j] == RC_MESH_NULL_IDX {
				break
			}
			v0 := t[j]
			v1 := t[j+1]
			if j+1 >= vertsPerPoly || t[j+1] == RC_MESH_NULL_IDX {
				v1 = t[0]
			}
			if v0 > v1 {
				for e := firstEdge[v1]; e != RC_MESH_NULL_IDX; e = nextEdge[e] {
					edge := edges[e]
					if edge.Vert[1] == v0 && edge.Poly[0] == edge.Poly[1] {
						edge.Poly[1] = uint16(i)
						edge.PolyEdge[1] = uint16(j)
						break
					}
				}
			}
		}
	}

	// Store adjacency
	for i := int32(0); i < edgeCount; i++ {
		e := edges[i]
		if e.Poly[0] != e.Poly[1] {
			p0 := common.GetVert2(polys, int32(e.Poly[0])*vertsPerPoly)
			p1 := common.GetVert2(polys, int32(e.Poly[1])*vertsPerPoly)
			p0[vertsPerPoly+int32(e.PolyEdge[0])] = e.Poly[1]
			p1[vertsPerPoly+int32(e.PolyEdge[1])] = e.Poly[0]
		}
	}
	return true
}

func computeVertexHash(x, y, z int32) int32 {
	h1 := uint32(0x8da6b343) // Large multiplicative constants;
	h2 := uint32(0xd8163841) // here arbitrarily chosen primes
	h3 := uint32(0xcb1ab31f)
	n := h1*uint32(x) + h2*uint32(y) + h3*uint32(z)
	return int32(n & (VERTEX_BUCKET_COUNT - 1))
}

func addVertex(x, y, z uint16,
	verts []uint16, firstVert []int32, nextVert []int32) (nv int32, t uint16) {
	bucket := computeVertexHash(int32(x), 0, int32(z))
	i := firstVert[bucket]

	for i != -1 {
		v := common.GetVert3(verts, i)
		if v[0] == x && (common.Abs(v[1]-y) <= 2) && v[2] == z {
			return nv, uint16(i)
		}

		i = nextVert[i] // next
	}

	// Could not find, create new.
	i = nv
	nv++
	v := common.GetVert3(verts, i)
	v[0] = x
	v[1] = y
	v[2] = z
	nextVert[i] = firstVert[bucket]
	firstVert[bucket] = i

	return nv, uint16(i)
}

func countPolyVerts(p []uint16, nvp int32) int32 {
	for i := int32(0); i < nvp; i++ {
		if p[i] == RC_MESH_NULL_IDX {
			return i
		}
	}

	return nvp
}

func getPolyMergeValue(pa, pb []uint16, verts []uint16, ea, eb int32, nvp int32) int32 {
	na := countPolyVerts(pa, nvp)
	nb := countPolyVerts(pb, nvp)

	// If the merged polygon would be too big, do not merge.
	if na+nb-2 > nvp {
		return -1
	}

	// Check if the polygons share an edge.
	ea = -1
	eb = -1

	for i := int32(0); i < na; i++ {
		va0 := pa[i]
		va1 := pa[(i+1)%na]
		if va0 > va1 {
			va0, va1 = va1, va0
		}

		for j := int32(0); j < nb; j++ {
			vb0 := pb[j]
			vb1 := pb[(j+1)%nb]
			if vb0 > vb1 {
				vb0, vb1 = vb1, vb0
			}

			if va0 == vb0 && va1 == vb1 {
				ea = i
				eb = j
				break
			}
		}
	}

	// No common edge, cannot merge.
	if ea == -1 || eb == -1 {
		return -1
	}

	// Check to see if the merged polygon would be convex.

	va := pa[(ea+na-1)%na]
	vb := pa[ea]
	vc := pb[(eb+2)%nb]

	if !common.Uleft(common.GetVert3(verts, va), common.GetVert3(verts, vb), common.GetVert3(verts, vc)) {
		return -1
	}

	va = pb[(eb+nb-1)%nb]
	vb = pb[eb]
	vc = pa[(ea+2)%na]
	if !common.Uleft(common.GetVert3(verts, va), common.GetVert3(verts, vb), common.GetVert3(verts, vc)) {
		return -1
	}

	va = pa[ea]
	vb = pa[(ea+1)%na]

	dx := int32(verts[va*3+0] - verts[vb*3+0])
	dy := int32(verts[va*3+2] - verts[vb*3+2])

	return dx*dx + dy*dy
}

func mergePolyVerts(pa, pb []uint16, ea, eb int32, tmp []uint16, nvp int32) {
	na := countPolyVerts(pa, nvp)
	nb := countPolyVerts(pb, nvp)

	for i := range tmp {
		if i < int(nvp) {
			tmp[i] = 0xff
		}
	}
	// Merge polygons.
	n := 0
	// Add pa
	for i := int32(0); i < na-1; i++ {
		tmp[n] = pa[(ea+1+i)%na]
		n++
	}

	// Add pb
	for i := int32(0); i < nb-1; i++ {
		tmp[n] = pb[(eb+1+i)%nb]
		n++
	}

	copy(pa, tmp[:nvp])
}

func canRemoveVertex(mesh *RcPolyMesh, rem uint16) bool {
	nvp := mesh.Nvp

	// Count number of polygons to remove.
	numTouchedVerts := 0
	numRemainingEdges := 0
	for i := nvp; i < mesh.Npolys; i++ {
		p := mesh.Polys[i*nvp*2 : i*nvp*2+2]
		nv := countPolyVerts(p, nvp)
		numRemoved := 0
		numVerts := 0
		for j := nvp; j < nv; j++ {
			if p[j] == rem {
				numTouchedVerts++
				numRemoved++
			}
			numVerts++
		}
		if numRemoved > 0 {
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

	// Find edges which share the removed vertex.
	maxEdges := numTouchedVerts * 2
	nedges := 0
	edges := make([]int32, maxEdges*3)
	for i := int32(0); i < mesh.Npolys; i++ {
		p := mesh.Polys[i*nvp*2 : i*nvp*2+2]
		nv := countPolyVerts(p, nvp)

		// Collect edges which touches the removed vertex.
		j := int32(0)
		k := nv - 1
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
					e := edges[m*3 : m*3+3]
					if e[1] == int32(b) {
						// Exists, increment vertex share count.
						e[2]++
						exists = true
					}
				}
				// Add new edge.
				if !exists {
					e := edges[nedges*3 : nedges*3+3]
					e[0] = int32(a)
					e[1] = int32(b)
					e[2] = 1
					nedges++
				}
			}
		}
		k = j
		j++
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

func removeVertex(mesh *RcPolyMesh, rem uint16, maxTris int32) bool {
	nvp := mesh.Nvp

	// Count number of polygons to remove.
	numRemovedVerts := int32(0)
	for i := int32(0); i < mesh.Npolys; i++ {
		p := mesh.Polys[i*nvp*2:]
		nv := countPolyVerts(p, nvp)
		for j := int32(0); j < nv; j++ {
			if p[j] == rem {
				numRemovedVerts++
			}

		}
	}

	nedges := int32(0)
	edges := make([]int32, numRemovedVerts*nvp*4)

	nhole := int32(0)
	hole := make([]int32, numRemovedVerts*nvp)
	nhreg := int32(0)
	hreg := make([]int32, numRemovedVerts*nvp)
	nharea := int32(0)
	harea := make([]int32, numRemovedVerts*nvp)

	for i := int32(0); i < mesh.Npolys; i++ {
		p := mesh.Polys[i*nvp*2:]
		nv := countPolyVerts(p, nvp)
		hasRem := false
		for j := int32(0); j < nv; j++ {
			if p[j] == rem {
				hasRem = true
			}
		}

		if hasRem {
			// Collect edges which does not touch the removed vertex.
			j := int32(0)
			k := nv - 1
			for j < nv {
				if p[j] != rem && p[k] != rem {
					e := common.GetVert4(edges, nedges)
					e[0] = int32(p[k])
					e[1] = int32(p[j])
					e[2] = int32(mesh.Regs[i])
					e[3] = int32(mesh.Areas[i])
					nedges++
				}
				k = j
				j++
			}
			// Remove the polygon.
			p2 := mesh.Polys[(mesh.Npolys-1)*nvp*2:]
			if !reflect.DeepEqual(p, p2) {
				copy(p, p2[:nvp])
			}
			for k := nvp; i < nvp+nvp; i++ {
				p[k] = 0xff
			}
			mesh.Regs[i] = mesh.Regs[mesh.Npolys-1]
			mesh.Areas[i] = mesh.Areas[mesh.Npolys-1]
			mesh.Npolys--
			i--
		}
	}

	// Remove vertex.
	for i := int32(rem); i < mesh.Nverts-1; i++ {
		mesh.Verts[i*3+0] = mesh.Verts[(i+1)*3+0]
		mesh.Verts[i*3+1] = mesh.Verts[(i+1)*3+1]
		mesh.Verts[i*3+2] = mesh.Verts[(i+1)*3+2]
	}
	mesh.Nverts--

	// Adjust indices to match the removed vertex layout.
	for i := int32(0); i < mesh.Npolys; i++ {
		p := mesh.Polys[i*nvp*2:]
		nv := countPolyVerts(p, nvp)
		for j := int32(0); j < nv; j++ {
			if p[j] > rem {
				p[j]--
			}
		}

	}
	for i := int32(0); i < nedges; i++ {
		if edges[i*4+0] > int32(rem) {
			edges[i*4+0]--
		}
		if edges[i*4+1] > int32(rem) {
			edges[i*4+1]--
		}
	}

	if nedges == 0 {
		return true
	}

	// Start with one vertex, keep appending connected
	// segments to the start and end of the hole.
	common.PushBack(edges[0], hole, &nhole)
	common.PushBack(edges[2], hreg, &nhreg)
	common.PushBack(edges[3], harea, &nharea)

	for nedges > 0 {
		match := false

		for i := int32(0); i < nedges; i++ {
			ea := edges[i*4+0]
			eb := edges[i*4+1]
			r := edges[i*4+2]
			a := edges[i*4+3]
			add := false
			if hole[0] == eb {
				// The segment matches the beginning of the hole boundary.
				common.PushFront(ea, hole, &nhole)
				common.PushFront(r, hreg, &nhreg)
				common.PushFront(a, harea, &nharea)
				add = true
			} else if hole[nhole-1] == ea {
				// The segment matches the end of the hole boundary.
				common.PushBack(eb, hole, &nhole)
				common.PushBack(r, hreg, &nhreg)
				common.PushBack(a, harea, &nharea)
				add = true
			}
			if add {
				// The edge segment was added, remove it.
				edges[i*4+0] = edges[(nedges-1)*4+0]
				edges[i*4+1] = edges[(nedges-1)*4+1]
				edges[i*4+2] = edges[(nedges-1)*4+2]
				edges[i*4+3] = edges[(nedges-1)*4+3]
				nedges--
				match = true
				i--
			}
		}

		if !match {
			break
		}

	}

	tris := make([]int32, nhole*3)
	tverts := make([]int32, nhole*4)
	thole := make([]int32, nhole)
	// Generate temp vertex array for triangulation.
	for i := int32(0); i < nhole; i++ {
		pi := hole[i]
		tverts[i*4+0] = int32(mesh.Verts[pi*3+0])
		tverts[i*4+1] = int32(mesh.Verts[pi*3+1])
		tverts[i*4+2] = int32(mesh.Verts[pi*3+2])
		tverts[i*4+3] = 0
		thole[i] = i
	}

	// Triangulate the hole.
	ntris := common.Triangulate(nhole, tverts, thole, tris)
	if ntris < 0 {
		ntris = -ntris
		log.Printf("removeVertex: triangulate() returned bad results.")
	}

	// Merge the hole triangles back to polygons.
	polys := make([]uint16, (ntris+1)*nvp)
	for i := range polys {
		polys[i] = 0xff
	}
	pregs := make([]uint16, ntris)
	pareas := make([]uint8, ntris)
	tmpPoly := polys[ntris*nvp:]
	// Build initial polygons.
	npolys := int32(0)
	for j := int32(0); j < ntris; j++ {
		t := common.GetVert3(tris, j)
		if t[0] != t[1] && t[0] != t[2] && t[1] != t[2] {
			polys[npolys*nvp+0] = uint16(hole[t[0]])
			polys[npolys*nvp+1] = uint16(hole[t[1]])
			polys[npolys*nvp+2] = uint16(hole[t[2]])

			// If this polygon covers multiple region types then
			// mark it as such
			if hreg[t[0]] != hreg[t[1]] || hreg[t[1]] != hreg[t[2]] {
				pregs[npolys] = RC_MULTIPLE_REGS
			} else {
				pregs[npolys] = uint16(hreg[t[0]])
			}

			pareas[npolys] = uint8(harea[t[0]])
			npolys++
		}
	}
	if npolys == 0 {
		return true
	}

	// Merge polygons.
	if nvp > 3 {
		for {
			// Find best polygons to merge.
			var bestMergeVal = int32(0)
			var bestPa = int32(0)
			var bestPb = int32(0)
			var bestEa = int32(0)
			var bestEb = int32(0)

			for j := int32(0); j < npolys-1; j++ {
				pj := polys[j*nvp:]
				for k := j + 1; k < npolys; k++ {
					pk := polys[k*nvp:]
					var ea, eb int32
					v := getPolyMergeValue(pj, pk, mesh.Verts, ea, eb, nvp)
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
				pa := polys[bestPa*nvp:]
				pb := polys[bestPb*nvp:]
				mergePolyVerts(pa, pb, bestEa, bestEb, tmpPoly, nvp)
				if pregs[bestPa] != pregs[bestPb] {
					pregs[bestPa] = RC_MULTIPLE_REGS
				}

				last := polys[(npolys-1)*nvp:]
				if !(reflect.DeepEqual(pb, last)) {
					copy(pb, last[:nvp])
				}

				pregs[bestPb] = pregs[npolys-1]
				pareas[bestPb] = pareas[npolys-1]
				npolys--
			} else {
				// Could not merge any polygons, stop.
				break
			}
		}
	}

	// Store polygons.
	for i := int32(0); i < npolys; i++ {
		if mesh.Npolys >= maxTris {
			break
		}
		p := mesh.Polys[mesh.Npolys*nvp*2:]
		for i := int32(0); i < nvp*2; i++ {
			p[i] = 0xff
		}
		for j := int32(0); j < nvp; j++ {
			p[j] = polys[i*nvp+j]
		}
		mesh.Regs[mesh.Npolys] = pregs[i]
		mesh.Areas[mesh.Npolys] = pareas[i]
		mesh.Npolys++
		if mesh.Npolys > maxTris {
			log.Printf("removeVertex: Too many polygons %d (max:%d).", mesh.Npolys, maxTris)
			return false
		}
	}

	return true
}

// / @par
// /
// / @note If the mesh data is to be used to construct a Detour navigation mesh, then the upper
// / limit must be restricted to <= #DT_VERTS_PER_POLYGON.
// /
// / @see rcAllocPolyMesh, RcContourSet, RcPolyMesh, RcConfig
func RcBuildPolyMesh(cset *RcContourSet, nvp int32, mesh *RcPolyMesh) bool {
	copy(mesh.Bmin, cset.Bmin)
	copy(mesh.Bmax, cset.Bmax)
	mesh.Cs = cset.Cs
	mesh.Ch = cset.Ch
	mesh.BorderSize = cset.BorderSize
	mesh.MaxEdgeError = cset.MaxError

	maxVertices := int32(0)
	maxTris := int32(0)
	maxVertsPerCont := int32(0)
	for i := int32(0); i < cset.Nconts; i++ {
		// Skip null contours.
		if cset.Conts[i].Nverts < 3 {
			continue
		}
		maxVertices += cset.Conts[i].Nverts
		maxTris += cset.Conts[i].Nverts - 2
		maxVertsPerCont = max(maxVertsPerCont, cset.Conts[i].Nverts)
	}

	if maxVertices >= 0xfffe {
		log.Printf("rcBuildPolyMesh: Too many vertices %d.", maxVertices)
		return false
	}

	vflags := make([]uint8, maxVertices)
	mesh.Verts = make([]uint16, maxVertices*3)

	mesh.Polys = make([]uint16, maxTris*nvp*2)
	for i := range mesh.Polys {
		mesh.Polys[i] = 0xff
	}
	mesh.Regs = make([]uint16, maxTris)
	mesh.Areas = make([]uint8, maxTris)

	mesh.Nverts = 0
	mesh.Npolys = 0
	mesh.Nvp = nvp
	mesh.Maxpolys = maxTris
	nextVert := make([]int32, maxVertices)
	firstVert := make([]int32, VERTEX_BUCKET_COUNT)
	for i := 0; i < VERTEX_BUCKET_COUNT; i++ {
		firstVert[i] = -1
	}

	indices := make([]int32, maxVertsPerCont)
	tris := make([]int32, maxVertsPerCont*3)
	polys := make([]uint16, (maxVertsPerCont+1)*nvp)
	tmpPoly := polys[maxVertsPerCont*nvp:]

	for i := int32(0); i < cset.Nconts; i++ {
		cont := cset.Conts[i]

		// Skip null contours.
		if cont.Nverts < 3 {
			continue
		}

		// Triangulate contour
		for j := int32(0); j < cont.Nverts; j++ {
			indices[j] = j
		}

		ntris := common.Triangulate(cont.Nverts, cont.Verts, indices, tris)
		if ntris <= 0 {
			// Bad triangulation, should not happen.
			/*			printf("\tconst float bmin[3] = {%ff,%ff,%ff};\n", cset.bmin[0], cset.bmin[1], cset.bmin[2]);
						printf("\tconst float cs = %ff;\n", cset.cs);
						printf("\tconst float ch = %ff;\n", cset.ch);
						printf("\tconst int verts[] = {\n");
						for (int k = 0; k < cont.nverts; ++k)
						{
							const int* v = &cont.verts[k*4];
							printf("\t\t%d,%d,%d,%d,\n", v[0], v[1], v[2], v[3]);
						}
						printf("\t};\n\tconst int nverts = sizeof(verts)/(sizeof(int)*4);\n");*/
			log.Printf("rcBuildPolyMesh: Bad triangulation Contour %d.", i)
			ntris = -ntris
		}

		// Add and merge vertices.
		for j := int32(0); j < cont.Nverts; j++ {
			v := common.GetVert4(cont.Verts, j)
			var tmp uint16
			mesh.Nverts, tmp = addVertex(uint16(v[0]), uint16(v[1]), uint16(v[2]), mesh.Verts, firstVert, nextVert)
			indices[j] = int32(tmp)
			if v[3]&RC_BORDER_VERTEX > 0 {
				// This vertex should be removed.
				vflags[indices[j]] = 1
			}
		}

		// Build initial polygons.
		npolys := int32(0)
		for i := range polys {
			polys[i] = 0xff
		}
		for j := int32(0); j < ntris; j++ {
			t := common.GetVert3(tris, j)
			if t[0] != t[1] && t[0] != t[2] && t[1] != t[2] {
				polys[npolys*nvp+0] = uint16(indices[t[0]])
				polys[npolys*nvp+1] = uint16(indices[t[1]])
				polys[npolys*nvp+2] = uint16(indices[t[2]])
				npolys++
			}
		}
		if npolys == 0 {
			continue
		}

		// Merge polygons.
		if nvp > 3 {
			for {
				// Find best polygons to merge.
				bestMergeVal := int32(0)
				bestPa := int32(0)
				bestPb := int32(0)
				bestEa := int32(0)
				bestEb := int32(0)

				for j := int32(0); j < npolys-1; j++ {
					pj := polys[j*nvp:]
					for k := j + 1; k < npolys; k++ {
						pk := polys[k*nvp:]
						var ea, eb int32
						v := getPolyMergeValue(pj, pk, mesh.Verts, ea, eb, nvp)
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
					pa := polys[bestPa*nvp:]
					pb := polys[bestPb*nvp:]
					mergePolyVerts(pa, pb, bestEa, bestEb, tmpPoly, nvp)
					lastPoly := polys[(npolys-1)*nvp:]
					if !reflect.DeepEqual(pb, lastPoly) {
						copy(pb, lastPoly[:nvp])
					}

					npolys--
				} else {
					// Could not merge any polygons, stop.
					break
				}
			}
		}

		// Store polygons.
		for j := int32(0); j < npolys; j++ {
			p := mesh.Polys[mesh.Npolys*nvp*2:]
			q := polys[j*nvp:]
			for k := int32(0); k < nvp; k++ {
				p[k] = q[k]
			}

			mesh.Regs[mesh.Npolys] = cont.Reg
			mesh.Areas[mesh.Npolys] = cont.Area
			mesh.Npolys++
			if mesh.Npolys > maxTris {
				log.Printf("rcBuildPolyMesh: Too many polygons %d (max:%d).", mesh.Npolys, maxTris)
				return false
			}
		}
	}

	// Remove edge vertices.
	for i := int32(0); i < mesh.Nverts; i++ {
		if vflags[i] != 0 {
			if !canRemoveVertex(mesh, uint16(i)) {
				continue
			}

			if !removeVertex(mesh, uint16(i), maxTris) {
				// Failed to remove vertex
				log.Printf("rcBuildPolyMesh: Failed to remove edge vertex %d.", i)
				return false
			}
			// Remove vertex
			// Note: mesh.nverts is already decremented inside removeVertex()!
			// Fixup vertex flags
			for j := i; j < mesh.Nverts; j++ {
				vflags[j] = vflags[j+1]
			}

			i--
		}
	}

	// Calculate adjacency.
	if !buildMeshAdjacency(mesh.Polys, mesh.Npolys, mesh.Nverts, nvp) {
		log.Printf("rcBuildPolyMesh: Adjacency failed.")
		return false
	}

	// Find portal edges
	if mesh.BorderSize > 0 {
		w := cset.Width
		h := cset.Height
		for i := int32(0); i < mesh.Npolys; i++ {
			p := mesh.Polys[i*2*nvp:]
			for j := int32(0); j < nvp; j++ {
				if p[j] == RC_MESH_NULL_IDX {
					break
				}
				// Skip connected edges.
				if p[nvp+j] != RC_MESH_NULL_IDX {
					continue
				}

				nj := j + 1
				if nj >= nvp || p[nj] == RC_MESH_NULL_IDX {
					nj = 0
				}
				va := common.GetVert3(mesh.Verts, p[j])
				vb := common.GetVert3(mesh.Verts, p[nj])

				if va[0] == 0 && vb[0] == 0 {
					p[nvp+j] = 0x8000 | 0
				} else if int32(va[2]) == h && int32(vb[2]) == h {
					p[nvp+j] = 0x8000 | 1
				} else if int32(va[0]) == w && int32(vb[0]) == w {
					p[nvp+j] = 0x8000 | 2
				} else if va[2] == 0 && vb[2] == 0 {
					p[nvp+j] = 0x8000 | 3
				}

			}
		}
	}

	// Just allocate the mesh flags array. The user is resposible to fill it.
	mesh.Flags = make([]uint16, mesh.Npolys)
	if mesh.Nverts > 0xffff {
		log.Printf("rcBuildPolyMesh: The resulting mesh has too many vertices %d (max %d). Data can be corrupted.", mesh.Nverts, 0xffff)
	}
	if mesh.Npolys > 0xffff {
		log.Printf("rcBuildPolyMesh: The resulting mesh has too many polygons %d (max %d). Data can be corrupted.", mesh.Npolys, 0xffff)
	}

	return true
}

// / @see rcAllocPolyMesh, RcPolyMesh
func rcMergePolyMeshes(meshes []*RcPolyMesh, nmeshes int32, mesh *RcPolyMesh) bool {

	if nmeshes == 0 || len(meshes) == 0 {
		return true
	}
	mesh.Nvp = meshes[0].Nvp
	mesh.Cs = meshes[0].Cs
	mesh.Ch = meshes[0].Ch
	copy(mesh.Bmin, meshes[0].Bmin)
	copy(mesh.Bmax, meshes[0].Bmax)

	maxVerts := int32(0)
	maxPolys := int32(0)
	maxVertsPerMesh := int32(0)
	for i := int32(0); i < nmeshes; i++ {
		common.Vmin(mesh.Bmin, meshes[i].Bmin)
		common.Vmax(mesh.Bmax, meshes[i].Bmax)
		maxVertsPerMesh = max(maxVertsPerMesh, meshes[i].Nverts)
		maxVerts += meshes[i].Nverts
		maxPolys += meshes[i].Npolys
	}

	mesh.Nverts = 0
	mesh.Verts = make([]uint16, maxVerts*3)

	mesh.Npolys = 0
	mesh.Polys = make([]uint16, maxPolys*2*mesh.Nvp)
	for i := range mesh.Polys {
		mesh.Polys[i] = 0xff
	}

	mesh.Regs = make([]uint16, maxPolys)

	mesh.Areas = make([]uint8, maxPolys)

	mesh.Flags = make([]uint16, maxPolys)

	nextVert := make([]int32, maxVerts)

	firstVert := make([]int32, VERTEX_BUCKET_COUNT)
	for i := 0; i < VERTEX_BUCKET_COUNT; i++ {
		firstVert[i] = -1
	}

	vremap := make([]uint16, maxVertsPerMesh)
	for i := int32(0); i < nmeshes; i++ {
		pmesh := meshes[i]

		ox := uint16(math.Floor(float64(pmesh.Bmin[0]-mesh.Bmin[0])/float64(mesh.Cs) + 0.5))
		oz := uint16(math.Floor(float64(pmesh.Bmin[2]-mesh.Bmin[2])/float64(mesh.Cs) + 0.5))

		isMinX := (ox == 0)
		isMinZ := (oz == 0)
		isMaxX := uint16(math.Floor(float64(mesh.Bmax[0]-pmesh.Bmax[0])/float64(mesh.Cs)+0.5)) == 0
		isMaxZ := uint16(math.Floor(float64(mesh.Bmax[2]-pmesh.Bmax[2])/float64(mesh.Cs)+0.5)) == 0
		isOnBorder := (isMinX || isMinZ || isMaxX || isMaxZ)

		for j := int32(0); j < pmesh.Nverts; j++ {
			v := common.GetVert3(pmesh.Verts, j)
			mesh.Nverts, vremap[j] = addVertex(v[0]+ox, v[1], v[2]+oz, mesh.Verts, firstVert, nextVert)
		}

		for j := int32(0); j < pmesh.Npolys; j++ {
			tgt := mesh.Polys[mesh.Npolys*2*mesh.Nvp:]
			src := pmesh.Polys[j*2*mesh.Nvp:]
			mesh.Regs[mesh.Npolys] = pmesh.Regs[j]
			mesh.Areas[mesh.Npolys] = pmesh.Areas[j]
			mesh.Flags[mesh.Npolys] = pmesh.Flags[j]
			mesh.Npolys++
			for k := int32(0); k < mesh.Nvp; k++ {
				if src[k] == RC_MESH_NULL_IDX {
					break
				}
				tgt[k] = vremap[src[k]]
			}

			if isOnBorder {
				for k := mesh.Nvp; k < mesh.Nvp*2; k++ {
					if src[k]&0x8000 > 0 && src[k] != 0xffff {
						dir := src[k] & 0xf
						switch dir {
						case 0: // Portal x-
							if isMinX {
								tgt[k] = src[k]
							}
							break
						case 1: // Portal z+
							if isMaxZ {
								tgt[k] = src[k]
							}
							break
						case 2: // Portal x+
							if isMaxX {
								tgt[k] = src[k]
							}

							break
						case 3: // Portal z-
							if isMinZ {
								tgt[k] = src[k]
							}

							break
						}
					}
				}
			}
		}
	}

	// Calculate adjacency.
	if !buildMeshAdjacency(mesh.Polys, mesh.Npolys, mesh.Nverts, mesh.Nvp) {
		log.Printf("rcMergePolyMeshes: Adjacency failed.")
		return false
	}

	if mesh.Nverts > 0xffff {
		log.Printf("rcMergePolyMeshes: The resulting mesh has too many vertices %d (max %d). Data can be corrupted.", mesh.Nverts, 0xffff)
	}
	if mesh.Npolys > 0xffff {
		log.Printf("rcMergePolyMeshes: The resulting mesh has too many polygons %d (max %d). Data can be corrupted.", mesh.Npolys, 0xffff)
	}

	return true
}

func rcCopyPolyMesh(src *RcPolyMesh, dst *RcPolyMesh) bool {
	// Destination must be empty.
	if len(dst.Verts) == 0 {
		panic("")
	}
	if len(dst.Polys) == 0 {
		panic("")
	}
	if len(dst.Regs) == 0 {
		panic("")
	}
	if len(dst.Areas) == 0 {
		panic("")
	}
	if len(dst.Flags) == 0 {
		panic("")
	}

	dst.Nverts = src.Nverts
	dst.Npolys = src.Npolys
	dst.Maxpolys = src.Npolys
	dst.Nvp = src.Nvp
	copy(dst.Bmin, src.Bmin)
	copy(dst.Bmax, src.Bmax)
	dst.Cs = src.Cs
	dst.Ch = src.Ch
	dst.BorderSize = src.BorderSize
	dst.MaxEdgeError = src.MaxEdgeError

	dst.Verts = make([]uint16, src.Nverts*3)
	copy(dst.Verts, src.Verts[:src.Nverts*3])

	dst.Polys = make([]uint16, src.Npolys*2*src.Nvp)
	copy(dst.Polys, src.Polys[:src.Npolys*2*src.Nvp])

	dst.Regs = make([]uint16, src.Npolys)
	copy(dst.Regs, src.Regs[:src.Npolys])

	dst.Areas = make([]uint8, src.Npolys)
	copy(dst.Areas, src.Areas[:src.Npolys])

	dst.Flags = make([]uint16, src.Npolys)

	copy(dst.Flags, src.Flags[:src.Npolys])

	return true
}
