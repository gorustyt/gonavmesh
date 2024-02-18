package recast

import (
	"github.com/gorustyt/gonavmesh/common"
	"sort"
)

const (
	/// Contour build flags.
	/// @see rcBuildContours
	/// Applied to the region id field of contour vertices in order to extract the region id.
	/// The region id field of a vertex may have several flags applied to it.  So the
	/// fields value can't be used directly.
	/// @see RcContour::verts, RcContour::rverts
	RC_CONTOUR_REG_MASK = 0xffff
	/// Area border flag.
	/// If a region ID has this bit set, then the associated element lies on
	/// the border of an area.
	/// (Used during the region and contour build process.)
	/// @see RcCompactSpan::reg, #RcContour::verts, #RcContour::rverts
	RC_AREA_BORDER = 0x20000
	/// Border vertex flag.
	/// If a region ID has this bit set, then the associated element lies on
	/// a tile border. If a contour vertex's region ID has this bit set, the
	/// vertex will later be removed in order to match the segments and vertices
	/// at tile boundaries.
	/// (Used during the build process.)
	/// @see RcCompactSpan::reg, #RcContour::verts, #RcContour::rverts
	RC_BORDER_VERTEX = 0x10000

	RC_CONTOUR_TESS_WALL_EDGES = 0x01 ///< Tessellate solid (impassable) edges during contour simplification.
	RC_CONTOUR_TESS_AREA_EDGES = 0x02 ///< Tessellate edges between areas during contour simplification.
)

func getCornerHeight(x, y, i, dir int32, chf *RcCompactHeightfield,
	isBorderVertex *bool) int32 {
	s := chf.Spans[i]
	ch := int32(s.Y)
	dirp := (dir + 1) & 0x3

	regs := []uint32{0, 0, 0, 0}

	// Combine region and area codes in order to prevent
	// border vertices which are in between two areas to be removed.
	regs[0] = uint32(chf.Spans[i].Reg) | (uint32(chf.Areas[i]) << 16)

	if rcGetCon(s, dir) != RC_NOT_CONNECTED {
		ax := x + common.GetDirOffsetX(dir)
		ay := y + common.GetDirOffsetY(dir)
		ai := int32(chf.Cells[ax+ay*chf.Width].Index) + rcGetCon(s, dir)
		as := chf.Spans[ai]
		ch = max(ch, int32(as.Y))
		regs[1] = uint32(chf.Spans[ai].Reg) | (uint32(chf.Areas[ai]) << 16)
		if rcGetCon(as, dirp) != RC_NOT_CONNECTED {
			ax2 := ax + common.GetDirOffsetX(dirp)
			ay2 := ay + common.GetDirOffsetY(dirp)
			ai2 := int32(chf.Cells[ax2+ay2*chf.Width].Index) + rcGetCon(as, dirp)
			as2 := chf.Spans[ai2]
			ch = max(ch, int32(as2.Y))
			regs[2] = uint32(chf.Spans[ai2].Reg) | (uint32(chf.Areas[ai2]) << 16)
		}
	}
	if rcGetCon(s, dirp) != RC_NOT_CONNECTED {
		ax := x + common.GetDirOffsetX(dirp)
		ay := y + common.GetDirOffsetY(dirp)
		ai := int32(chf.Cells[ax+ay*chf.Width].Index) + rcGetCon(s, dirp)
		as := chf.Spans[ai]
		ch = max(ch, int32(as.Y))
		regs[3] = uint32(chf.Spans[ai].Reg) | (uint32(chf.Areas[ai]) << 16)
		if rcGetCon(as, dir) != RC_NOT_CONNECTED {
			ax2 := ax + common.GetDirOffsetX(dir)
			ay2 := ay + common.GetDirOffsetY(dir)
			ai2 := int32(chf.Cells[ax2+ay2*chf.Width].Index) + rcGetCon(as, dir)
			as2 := chf.Spans[ai2]
			ch = max(ch, int32(as2.Y))
			regs[2] = uint32(chf.Spans[ai2].Reg) | (uint32(chf.Areas[ai2]) << 16)
		}
	}

	// Check if the vertex is special edge vertex, these vertices will be removed later.
	for j := 0; j < 4; j++ {
		a := j
		b := (j + 1) & 0x3
		c := (j + 2) & 0x3
		d := (j + 3) & 0x3

		// The vertex is a border vertex there are two same exterior cells in a row,
		// followed by two interior cells and none of the regions are out of bounds.
		twoSameExts := (regs[a]&regs[b]&RC_BORDER_REG) != 0 && regs[a] == regs[b]
		twoInts := ((regs[c] | regs[d]) & RC_BORDER_REG) == 0
		intsSameArea := (regs[c] >> 16) == (regs[d] >> 16)
		noZeros := regs[a] != 0 && regs[b] != 0 && regs[c] != 0 && regs[d] != 0
		if twoSameExts && twoInts && intsSameArea && noZeros {
			*isBorderVertex = true
			break
		}
	}

	return ch
}

func walkContour(x, y, i int32, chf *RcCompactHeightfield,
	flags []uint8, points Stack[int32]) {
	// Choose the first non-connected edge
	dir := int32(0)
	for (flags[i] & (1 << dir)) == 0 {
		dir++
	}

	startDir := dir
	starti := i

	area := chf.Areas[i]

	iter := 0
	for {
		iter++
		if iter >= 40000 {
			break
		}
		if flags[i]&(1<<dir) > 0 {
			// Choose the edge corner
			isBorderVertex := false
			isAreaBorder := false
			px := x
			py := getCornerHeight(x, y, i, dir, chf, &isBorderVertex)
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
			r := int32(0)
			s := chf.Spans[i]
			if rcGetCon(s, dir) != RC_NOT_CONNECTED {
				ax := x + common.GetDirOffsetX(dir)
				ay := y + common.GetDirOffsetY(dir)
				ai := int32(chf.Cells[ax+ay*chf.Width].Index) + rcGetCon(s, dir)
				r = int32(chf.Spans[ai].Reg)
				if area != chf.Areas[ai] {
					isAreaBorder = true
				}

			}
			if isBorderVertex {
				r |= RC_BORDER_VERTEX
			}

			if isAreaBorder {
				r |= RC_AREA_BORDER
			}

			points.Push(px)
			points.Push(py)
			points.Push(pz)
			points.Push(r)

			flags[i] &= ^(1 << dir) // Remove visited edges
			dir = (dir + 1) & 0x3   // Rotate CW
		} else {
			ni := int32(-1)
			nx := x + common.GetDirOffsetX(dir)
			ny := y + common.GetDirOffsetY(dir)
			s := chf.Spans[i]
			if rcGetCon(s, dir) != RC_NOT_CONNECTED {
				nc := chf.Cells[nx+ny*chf.Width]
				ni = int32(nc.Index) + rcGetCon(s, dir)
			}
			if ni == -1 {
				// Should not happen.
				return
			}
			x = nx
			y = ny
			i = ni
			dir = (dir + 3) & 0x3 // Rotate CCW
		}

		if starti == i && startDir == dir {
			break
		}
	}
}

func simplifyContour(points Stack[int32], simplified Stack[int32], maxError float32, maxEdgeLen, buildFlags int32) {
	// Add initial points.
	hasConnections := false
	for i := 0; i < points.Len(); i += 4 {
		if (points.Index(i+3) & RC_CONTOUR_REG_MASK) != 0 {
			hasConnections = true
			break
		}
	}

	if hasConnections {
		// The contour has some portals to other regions.
		// Add a new point to every location where the region changes.
		i := 0
		ni := points.Len() / 4
		for ; i < ni; i++ {
			ii := (i + 1) % ni
			differentRegs := (points.Index(i*4+3) & RC_CONTOUR_REG_MASK) != (points.Index(ii*4+3) & RC_CONTOUR_REG_MASK)
			areaBorders := (points.Index(i*4+3) & RC_AREA_BORDER) != (points.Index(ii*4+3) & RC_AREA_BORDER)
			if differentRegs || areaBorders {
				simplified.Push(points.Index(i*4 + 0))
				simplified.Push(points.Index(i*4 + 1))
				simplified.Push(points.Index(i*4 + 2))
				simplified.Push(int32(i))
			}
		}
	}

	if simplified.Len() == 0 {
		// If there is no connections at all,
		// create some initial points for the simplification process.
		// Find lower-left and upper-right vertices of the contour.
		llx := points.Index(0)
		lly := points.Index(1)
		llz := points.Index(2)
		lli := int32(0)
		urx := points.Index(0)
		ury := points.Index(1)
		urz := points.Index(2)
		uri := int32(0)
		for i := 0; i < points.Len(); i += 4 {
			x := points.Index(i + 0)
			y := points.Index(i + 1)
			z := points.Index(i + 2)
			if x < llx || (x == llx && z < llz) {
				llx = x
				lly = y
				llz = z
				lli = int32(i) / 4
			}
			if x > urx || (x == urx && z > urz) {
				urx = x
				ury = y
				urz = z
				uri = int32(i) / 4
			}
		}
		simplified.Push(llx)
		simplified.Push(lly)
		simplified.Push(llz)
		simplified.Push(lli)

		simplified.Push(urx)
		simplified.Push(ury)
		simplified.Push(urz)
		simplified.Push(uri)
	}

	// Add points until all raw points are within
	// error tolerance to the simplified shape.
	pn := points.Len() / 4
	for i := 0; i < simplified.Len()/4; {
		ii := (i + 1) % (simplified.Len() / 4)

		ax := simplified.Index(i*4 + 0)
		az := simplified.Index(i*4 + 2)
		ai := simplified.Index(i*4 + 3)

		bx := simplified.Index(ii*4 + 0)
		bz := simplified.Index(ii*4 + 2)
		bi := simplified.Index(ii*4 + 3)

		// Find maximum deviation from the segment.
		maxd := float32(0)
		maxi := int32(-1)
		var ci, cinc, endi int

		// Traverse the segment in lexilogical order so that the
		// max deviation is calculated similarly when traversing
		// opposite segments.
		if bx > ax || (bx == ax && bz > az) {
			cinc = 1
			ci = (int(ai) + cinc) % pn
			endi = int(bi)
		} else {
			cinc = pn - 1
			ci = (int(bi) + cinc) % pn
			endi = int(ai)
			ax, bx = bx, ax
			az, bz = bz, az
		}

		// Tessellate only outer edges or edges between areas.
		if (points.Index(ci*4+3)&RC_CONTOUR_REG_MASK) == 0 || (points.Index(ci*4+3)&RC_AREA_BORDER > 0) {
			for ci != endi {
				d := common.ContourdistancePtSeg(points.Index(ci*4+0), points.Index(ci*4+2), ax, az, bx, bz)
				if d > maxd {
					maxd = d
					maxi = int32(ci)
				}
				ci = (ci + cinc) % pn
			}
		}

		// If the max deviation is larger than accepted error,
		// add new point, else continue to next segment.
		if maxi != -1 && maxd > (maxError*maxError) {
			// Add space for the new point.
			simplified.Resize(simplified.Len() + 4)
			n := simplified.Len() / 4
			for j := n - 1; j > i; j-- {
				simplified.SetByIndex(j*4+0, simplified.Index((j-1)*4+0))
				simplified.SetByIndex(j*4+1, simplified.Index((j-1)*4+1))
				simplified.SetByIndex(j*4+2, simplified.Index((j-1)*4+2))
				simplified.SetByIndex(j*4+3, simplified.Index((j-1)*4+3))
			}
			// Add the point.
			simplified.SetByIndex((i+1)*4+0, points.Index(int(maxi)*4+0))
			simplified.SetByIndex((i+1)*4+1, points.Index(int(maxi)*4+1))
			simplified.SetByIndex((i+1)*4+2, points.Index(int(maxi)*4+2))
			simplified.SetByIndex((i+1)*4+3, maxi)
		} else {
			i++
		}
	}

	// Split too long edges.
	if maxEdgeLen > 0 && (buildFlags&(RC_CONTOUR_TESS_WALL_EDGES|RC_CONTOUR_TESS_AREA_EDGES)) != 0 {
		for i := 0; i < simplified.Len()/4; {
			ii := (i + 1) % (simplified.Len() / 4)

			ax := simplified.Index(i*4 + 0)
			az := simplified.Index(i*4 + 2)
			ai := simplified.Index(i*4 + 3)

			bx := simplified.Index(ii*4 + 0)
			bz := simplified.Index(ii*4 + 2)
			bi := simplified.Index(ii*4 + 3)

			// Find maximum deviation from the segment.
			maxi := -1
			ci := int(ai+1) % pn

			// Tessellate only outer edges or edges between areas.
			tess := false
			// Wall edges.
			if (buildFlags&RC_CONTOUR_TESS_WALL_EDGES > 0) && (points.Index(ci*4+3)&RC_CONTOUR_REG_MASK) == 0 {
				tess = true
			}

			// Edges between areas.
			if (buildFlags&RC_CONTOUR_TESS_AREA_EDGES > 0) && (points.Index(ci*4+3)&RC_AREA_BORDER > 0) {
				tess = true
			}

			if tess {
				dx := bx - ax
				dz := bz - az
				if dx*dx+dz*dz > maxEdgeLen*maxEdgeLen {
					// Round based on the segments in lexilogical order so that the
					// max tesselation is consistent regardless in which direction
					// segments are traversed.
					n := (bi - ai)
					if bi < ai {
						n = bi + int32(pn) - ai
					}
					if n > 1 {
						if bx > ax || (bx == ax && bz > az) {
							maxi = int(ai+n/2) % pn
						} else {
							maxi = int(ai+(n+1)/2) % pn
						}

					}
				}
			}

			// If the max deviation is larger than accepted error,
			// add new point, else continue to next segment.
			if maxi != -1 {
				// Add space for the new point.
				simplified.Resize(simplified.Len() + 4)
				n := simplified.Len() / 4
				for j := n - 1; j > i; j-- {
					simplified.SetByIndex(j*4+0, simplified.Index((j-1)*4+0))
					simplified.SetByIndex(j*4+1, simplified.Index((j-1)*4+1))
					simplified.SetByIndex(j*4+2, simplified.Index((j-1)*4+2))
					simplified.SetByIndex(j*4+3, simplified.Index((j-1)*4+3))
				}
				// Add the point.
				simplified.SetByIndex((i+1)*4+0, points.Index(maxi*4+0))
				simplified.SetByIndex((i+1)*4+1, points.Index(maxi*4+1))
				simplified.SetByIndex((i+1)*4+2, points.Index(maxi*4+2))
				simplified.SetByIndex((i+1)*4+3, int32(maxi))
			} else {
				i++
			}
		}
	}

	for i := 0; i < simplified.Len()/4; i++ {
		// The edge vertex flag is take from the current raw point,
		// and the neighbour region is take from the next raw point.
		ai := (simplified.Index(i*4+3) + 1) % int32(pn)
		bi := simplified.Index(i*4 + 3)
		v := (points.Index(int(ai*4+3)) & (RC_CONTOUR_REG_MASK | RC_AREA_BORDER)) | (points.Index(int(bi)*4+3) & RC_BORDER_VERTEX)
		simplified.SetByIndex(i*4+3, v)
	}

}

func removeDegenerateSegments(simplified Stack[int32]) {
	// Remove adjacent vertices which are equal on xz-plane,
	// or else the triangulator will get confused.
	npts := simplified.Len() / 4
	for i := 0; i < npts; i++ {
		ni := common.Next(i, npts)
		a := []int{int(simplified.Index(i * 4)), int(simplified.Index(i*4 + 1))}
		b := []int{int(simplified.Index(ni * 4)), int(simplified.Index(ni*4 + 1))}
		if common.RcVequal(a, b) {
			// Degenerate segment, remove.
			for j := i; j < simplified.Len()/4-1; j++ {
				simplified.SetByIndex(j*4+0, simplified.Index((j+1)*4+0))
				simplified.SetByIndex(j*4+1, simplified.Index((j+1)*4+1))
				simplified.SetByIndex(j*4+2, simplified.Index((j+1)*4+2))
				simplified.SetByIndex(j*4+3, simplified.Index((j+1)*4+3))
			}
			simplified.Resize(simplified.Len() - 4)
			npts--
		}
	}
}
func mergeContours(ca, cb *RcContour, ia, ib int32) bool {
	maxVerts := ca.Nverts + cb.Nverts + 2
	verts := make([]int32, maxVerts*4)

	nv := int32(0)

	// Copy contour A.
	for i := int32(0); i <= ca.Nverts; i++ {
		dst := common.GetVert4(verts, nv)
		src := common.GetVert4(ca.Verts, ((ia + i) % ca.Nverts))
		dst[0] = src[0]
		dst[1] = src[1]
		dst[2] = src[2]
		dst[3] = src[3]
		nv++
	}

	// Copy contour B
	for i := int32(0); i <= cb.Nverts; i++ {
		dst := common.GetVert4(verts, nv)
		src := common.GetVert4(ca.Verts, ((ib + i) % cb.Nverts))
		dst[0] = src[0]
		dst[1] = src[1]
		dst[2] = src[2]
		dst[3] = src[3]
		nv++
	}

	ca.Verts = verts
	ca.Nverts = nv

	cb.Verts = []int32{}
	cb.Nverts = 0

	return true
}

type RcContourHole struct {
	contour              *RcContour
	minx, minz, leftmost int32
}

type RcContourRegion struct {
	outline *RcContour
	holes   []*RcContourHole
	nholes  int32
}

type rcPotentialDiagonal struct {
	vert int32
	dist int32
}

// Finds the lowest leftmost vertex of a contour.
func findLeftMostVertex(contour *RcContour, minx, minz *int32, leftmost *int32) {
	*minx = contour.Verts[0]
	*minz = contour.Verts[2]
	*leftmost = 0
	for i := int32(1); i < contour.Nverts; i++ {
		x := contour.Verts[i*4+0]
		z := contour.Verts[i*4+2]
		if x < *minx || (x == *minx && z < *minz) {
			*minx = x
			*minz = z
			*leftmost = i
		}
	}
}

func compareHoles(va, vb *RcContourHole) int {
	a := va
	b := vb
	if a.minx == b.minx {
		if a.minz < b.minz {
			return -1
		}

		if a.minz > b.minz {
			return 1
		}

	} else {
		if a.minx < b.minx {
			return -1
		}

		if a.minx > b.minx {
			return 1
		}

	}
	return 0
}

func compareDiagDist(va, vb *rcPotentialDiagonal) int {
	a := va
	b := vb
	if a.dist < b.dist {
		return -1
	}

	if a.dist > b.dist {
		return 1
	}

	return 0
}

func mergeRegionHoles(region *RcContourRegion) {
	// Sort holes from left to right.
	for i := int32(0); i < region.nholes; i++ {
		findLeftMostVertex(region.holes[i].contour, &region.holes[i].minx, &region.holes[i].minz, &region.holes[i].leftmost)
	}
	sort.Slice(region.holes, func(i, j int) bool {
		return compareHoles(region.holes[i], region.holes[j]) > 0
	})
	maxVerts := region.outline.Nverts
	for i := int32(0); i < region.nholes; i++ {
		maxVerts += region.holes[i].contour.Nverts
	}

	diags := make([]*rcPotentialDiagonal, maxVerts)

	outline := region.outline

	// Merge holes into the outline one by one.
	for i := int32(0); i < region.nholes; i++ {
		hole := region.holes[i].contour

		index := int32(-1)
		bestVertex := region.holes[i].leftmost
		for iter := int32(0); iter < hole.Nverts; iter++ {
			// Find potential diagonals.
			// The 'best' vertex must be in the cone described by 3 consecutive vertices of the outline.
			// ..o j-1
			//   |
			//   |   * best
			//   |
			// j o-----o j+1
			//         :
			ndiags := 0
			corner := common.GetVert4(hole.Verts, bestVertex)
			for j := int32(0); j < outline.Nverts; j++ {
				if common.ContourInCone(j, outline.Nverts, outline.Verts, corner) {
					dx := outline.Verts[j*4+0] - corner[0]
					dz := outline.Verts[j*4+2] - corner[2]
					diags[ndiags].vert = j
					diags[ndiags].dist = dx*dx + dz*dz
					ndiags++
				}
			}
			// Sort potential diagonals by distance, we want to make the connection as short as possible.
			s := diags[:ndiags]
			sort.Slice(s, func(i, j int) bool {
				return compareDiagDist(s[i], s[j]) > 0
			})

			// Find a diagonal that is not intersecting the outline not the remaining holes.
			index = -1
			for j := 0; j < ndiags; j++ {
				pt := common.GetVert4(outline.Verts, diags[j].vert)
				intersect := common.IntersectSegContour(pt, corner, diags[i].vert, outline.Nverts, outline.Verts)
				for k := i; k < region.nholes && !intersect; k++ {
					if intersect || common.IntersectSegContour(pt, corner, -1, region.holes[k].contour.Nverts, region.holes[k].contour.Verts) {
						intersect = true
					}
				}

				if !intersect {
					index = diags[j].vert
					break
				}
			}
			// If found non-intersecting diagonal, stop looking.
			if index != -1 {
				break
			}

			// All the potential diagonals for the current vertex were intersecting, try next vertex.
			bestVertex = (bestVertex + 1) % hole.Nverts
		}

		if index == -1 {
			continue
		}
		if !mergeContours(region.outline, hole, index, bestVertex) {
			continue
		}
	}
}

// / @par
// /
// / The raw contours will match the region outlines exactly. The @p maxError and @p maxEdgeLen
// / parameters control how closely the simplified contours will match the raw contours.
// /
// / Simplified contours are generated such that the vertices for portals between areas match up.
// / (They are considered mandatory vertices.)
// /
// / Setting @p maxEdgeLength to zero will disabled the edge length feature.
// /
// / See the #RcConfig documentation for more information on the configuration parameters.
// /
// / @see rcAllocContourSet, RcCompactHeightfield, RcContourSet, RcConfig
func RcBuildContours(chf *RcCompactHeightfield,
	maxError float32, maxEdgeLen int32,
	cset *RcContourSet, buildFlagss ...int32) bool {
	buildFlags := int32(RC_CONTOUR_TESS_WALL_EDGES)
	if len(buildFlagss) > 0 {
		buildFlags = buildFlagss[0]
	}
	w := chf.Width
	h := chf.Height
	borderSize := chf.BorderSize

	copy(cset.Bmin, chf.Bmin[:])
	copy(cset.Bmax, chf.Bmax[:])
	if borderSize > 0 {
		// If the heightfield was build with bordersize, remove the offset.
		pad := float32(borderSize) * chf.Cs
		cset.Bmin[0] += pad
		cset.Bmin[2] += pad
		cset.Bmax[0] -= pad
		cset.Bmax[2] -= pad
	}
	cset.Cs = chf.Cs
	cset.Ch = chf.Ch
	cset.Width = chf.Width - chf.BorderSize*2
	cset.Height = chf.Height - chf.BorderSize*2
	cset.BorderSize = chf.BorderSize
	cset.MaxError = maxError

	maxContours := common.Max(chf.MaxRegions, 8)

	cset.Conts = make([]*RcContour, maxContours)
	cset.Nconts = 0
	flags := make([]uint8, chf.SpanCount)
	// Mark boundaries.
	for y := int32(0); y < h; y++ {
		for x := int32(0); x < w; x++ {
			c := chf.Cells[x+y*w]
			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni; i++ {
				res := uint8(0)
				s := chf.Spans[i]
				if chf.Spans[i].Reg == 0 || (chf.Spans[i].Reg&RC_BORDER_REG > 0) {
					flags[i] = 0
					continue
				}
				for dir := int32(0); dir < 4; dir++ {
					r := int32(0)
					if rcGetCon(s, dir) != RC_NOT_CONNECTED {
						ax := x + common.GetDirOffsetX(dir)
						ay := y + common.GetDirOffsetY(dir)
						ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(s, dir)
						r = int32(chf.Spans[ai].Reg)
					}
					if r == int32(chf.Spans[i].Reg) {
						res |= (1 << dir)
					}

				}
				flags[i] = res ^ 0xf // Inverse, mark non connected edges.
			}
		}
	}

	verts := NewStackArray(func() int32 { return 0 }, 256)
	simplified := NewStackArray(func() int32 { return 0 }, 64)

	for y := int32(0); y < h; y++ {
		for x := int32(0); x < w; x++ {
			c := chf.Cells[x+y*w]
			i := int32(c.Index)
			ni := int32(c.Index + c.Count)
			for ; i < ni; i++ {
				if flags[i] == 0 || flags[i] == 0xf {
					flags[i] = 0
					continue
				}
				reg := chf.Spans[i].Reg
				if reg == 0 || (reg&RC_BORDER_REG > 0) {
					continue
				}

				area := chf.Areas[i]

				verts.Clear()
				simplified.Clear()

				walkContour(x, y, i, chf, flags, verts)

				simplifyContour(verts, simplified, maxError, maxEdgeLen, buildFlags)
				removeDegenerateSegments(simplified)

				// Store region.contour remap info.
				// Create contour.
				if simplified.Len()/4 >= 3 {
					if cset.Nconts >= int32(maxContours) {
						// Allocate more contours.
						// This happens when a region has holes.
						maxContours *= 2
						newConts := make([]*RcContour, maxContours)
						for j := int32(0); j < cset.Nconts; j++ {
							newConts[j] = cset.Conts[j]
							// Reset source pointers to prevent data deletion.
							cset.Conts[j].Verts = []int32{}
							cset.Conts[j].Rverts = []int32{}
						}
						cset.Conts = newConts
					}

					cont := cset.Conts[cset.Nconts]
					cset.Nconts++

					cont.Nverts = int32(simplified.Len()) / 4
					cont.Verts = make([]int32, cont.Nverts*4)
					copy(cont.Verts, simplified.Slice(0, int(cont.Nverts)*4))
					if borderSize > 0 {
						// If the heightfield was build with bordersize, remove the offset.
						for j := int32(0); j < cont.Nverts; j++ {
							v := common.GetVert4(cont.Verts, j)
							v[0] -= borderSize
							v[2] -= borderSize
						}
					}

					cont.Nrverts = verts.Len() / 4
					cont.Rverts = make([]int32, cont.Nrverts*4)

					copy(cont.Rverts, verts.Slice(0, cont.Nrverts*4))
					if borderSize > 0 {
						// If the heightfield was build with bordersize, remove the offset.
						for j := 0; j < cont.Nrverts; j++ {
							v := common.GetVert4(cont.Rverts, j)
							v[0] -= borderSize
							v[2] -= borderSize
						}
					}

					cont.Reg = reg
					cont.Area = area
				}
			}
		}
	}

	// Merge holes if needed.
	if cset.Nconts > 0 {
		// Calculate winding of all polygons.
		winding := make([]int, cset.Nconts)
		nholes := 0
		for i := int32(0); i < cset.Nconts; i++ {
			cont := cset.Conts[i]
			// If the contour is wound backwards, it is a hole.
			winding[i] = 1
			if common.CalcAreaOfPolygon2D(cont.Verts, cont.Nverts) < 0 {
				winding[i] = -1
			}
			if winding[i] < 0 {
				nholes++
			}

		}

		if nholes > 0 {
			// Collect outline contour and holes contours per region.
			// We assume that there is one outline and multiple holes.
			nregions := chf.MaxRegions + 1
			regions := make([]*RcContourRegion, nregions)
			for i := range regions {
				regions[i] = &RcContourRegion{}
			}
			holes := make([]*RcContourHole, cset.Nconts)
			for i := range holes {
				holes[i] = &RcContourHole{}
			}
			for i := int32(0); i < cset.Nconts; i++ {
				cont := cset.Conts[i]
				// Positively would contours are outlines, negative holes.
				if winding[i] > 0 {
					if regions[cont.Reg].outline != nil {
						regions[cont.Reg].outline = cont
					}

				} else {
					regions[cont.Reg].nholes++
				}
			}
			index := int32(0)
			for i := uint16(0); i < nregions; i++ {
				if regions[i].nholes > 0 {
					regions[i].holes = holes[index:]
					index += regions[i].nholes
					regions[i].nholes = 0
				}
			}
			for i := int32(0); i < cset.Nconts; i++ {
				cont := cset.Conts[i]
				reg := regions[cont.Reg]
				if winding[i] < 0 {
					reg.holes[reg.nholes].contour = cont
					reg.nholes++
				}

			}

			// Finally merge each regions holes into the outline.
			for i := uint16(0); i < nregions; i++ {
				reg := regions[i]
				if reg.nholes == 0 {
					continue
				}

				if reg.outline == nil {
					mergeRegionHoles(reg)
				} else {
					// The region does not have an outline.
					// This can happen if the contour becaomes selfoverlapping because of
					// too aggressive simplification settings.
				}
			}
		}

	}

	return true
}
