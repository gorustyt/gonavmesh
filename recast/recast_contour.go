package recast

import "sort"

const (
	/// Contour build flags.
	/// @see rcBuildContours
	/// Applied to the region id field of contour vertices in order to extract the region id.
	/// The region id field of a vertex may have several flags applied to it.  So the
	/// fields value can't be used directly.
	/// @see rcContour::verts, rcContour::rverts
	RC_CONTOUR_REG_MASK = 0xffff
	/// Area border flag.
	/// If a region ID has this bit set, then the associated element lies on
	/// the border of an area.
	/// (Used during the region and contour build process.)
	/// @see rcCompactSpan::reg, #rcContour::verts, #rcContour::rverts
	RC_AREA_BORDER = 0x20000
	/// Border vertex flag.
	/// If a region ID has this bit set, then the associated element lies on
	/// a tile border. If a contour vertex's region ID has this bit set, the
	/// vertex will later be removed in order to match the segments and vertices
	/// at tile boundaries.
	/// (Used during the build process.)
	/// @see rcCompactSpan::reg, #rcContour::verts, #rcContour::rverts
	RC_BORDER_VERTEX = 0x10000

	RC_CONTOUR_TESS_WALL_EDGES = 0x01 ///< Tessellate solid (impassable) edges during contour simplification.
	RC_CONTOUR_TESS_AREA_EDGES = 0x02 ///< Tessellate edges between areas during contour simplification.
)

func contourInCone(i, n int, verts, pj []int) bool {
	pi := rcGetVert4(verts, i)
	pi1 := rcGetVert4(verts, next(i, n))
	pin1 := rcGetVert4(verts, prev(i, n))

	// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
	if leftOn(pin1, pi, pi1) {
		return left(pi, pj, pin1) && left(pj, pi, pi1)
	}

	// Assume (i-1,i,i+1) not collinear.
	// else P[i] is reflex.
	return !(leftOn(pi, pj, pi1) && leftOn(pj, pi, pin1))
}

func intersectSegContour(d0, d1 []int, i, n int, verts []int) bool {
	// For each edge (k,k+1) of P
	for k := 0; k < n; k++ {
		k1 := next(k, n)
		// Skip edges incident to i.
		if i == k || i == k1 {
			continue
		}

		p0 := rcGetVert4(verts, k)
		p1 := rcGetVert4(verts, k1)
		if vequal(d0, p0) || vequal(d1, p0) || vequal(d0, p1) || vequal(d1, p1) {
			continue
		}

		if intersect(d0, d1, p0, p1) {
			return true
		}

	}
	return false
}
func getCornerHeight(x, y, i, dir int, chf *rcCompactHeightfield,
	isBorderVertex *bool) int {
	s := chf.spans[i]
	ch := s.y
	dirp := (dir + 1) & 0x3

	regs := []int{0, 0, 0, 0}

	// Combine region and area codes in order to prevent
	// border vertices which are in between two areas to be removed.
	regs[0] = chf.spans[i].reg | (chf.areas[i] << 16)

	if rcGetCon(s, dir) != RC_NOT_CONNECTED {
		ax := x + rcGetDirOffsetX(dir)
		ay := y + rcGetDirOffsetY(dir)
		ai := chf.cells[ax+ay*chf.width].index + rcGetCon(s, dir)
		as := chf.spans[ai]
		ch = rcMax(ch, as.y)
		regs[1] = chf.spans[ai].reg | (chf.areas[ai] << 16)
		if rcGetCon(as, dirp) != RC_NOT_CONNECTED {
			ax2 := ax + rcGetDirOffsetX(dirp)
			ay2 := ay + rcGetDirOffsetY(dirp)
			ai2 := chf.cells[ax2+ay2*chf.width].index + rcGetCon(as, dirp)
			as2 := chf.spans[ai2]
			ch = rcMax(ch, as2.y)
			regs[2] = chf.spans[ai2].reg | (chf.areas[ai2] << 16)
		}
	}
	if rcGetCon(s, dirp) != RC_NOT_CONNECTED {
		ax := x + rcGetDirOffsetX(dirp)
		ay := y + rcGetDirOffsetY(dirp)
		ai := chf.cells[ax+ay*chf.width].index + rcGetCon(s, dirp)
		as := chf.spans[ai]
		ch = rcMax(ch, as.y)
		regs[3] = chf.spans[ai].reg | (chf.areas[ai] << 16)
		if rcGetCon(as, dir) != RC_NOT_CONNECTED {
			ax2 := ax + rcGetDirOffsetX(dir)
			ay2 := ay + rcGetDirOffsetY(dir)
			ai2 := chf.cells[ax2+ay2*chf.width].index + rcGetCon(as, dir)
			as2 := chf.spans[ai2]
			ch = rcMax(ch, as2.y)
			regs[2] = chf.spans[ai2].reg | (chf.areas[ai2] << 16)
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

func contourdistancePtSeg(x, z, px, pz, qx, qz int) float64 {
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
func walkContour(x, y, i int, chf *rcCompactHeightfield,
	flags []int, points Stack[int]) {
	// Choose the first non-connected edge
	dir := 0
	for (flags[i] & (1 << dir)) == 0 {
		dir++
	}

	startDir := dir
	starti := i

	area := chf.areas[i]

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
			r := 0
			s := chf.spans[i]
			if rcGetCon(s, dir) != RC_NOT_CONNECTED {
				ax := x + rcGetDirOffsetX(dir)
				ay := y + rcGetDirOffsetY(dir)
				ai := chf.cells[ax+ay*chf.width].index + rcGetCon(s, dir)
				r = chf.spans[ai].reg
				if area != chf.areas[ai] {
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
			ni := -1
			nx := x + rcGetDirOffsetX(dir)
			ny := y + rcGetDirOffsetY(dir)
			s := chf.spans[i]
			if rcGetCon(s, dir) != RC_NOT_CONNECTED {
				nc := chf.cells[nx+ny*chf.width]
				ni = nc.index + rcGetCon(s, dir)
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

func simplifyContour(points Stack[int], simplified Stack[int], maxError float64, maxEdgeLen, buildFlags int) {
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
				simplified.Push(i)
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
		lli := 0
		urx := points.Index(0)
		ury := points.Index(1)
		urz := points.Index(2)
		uri := 0
		for i := 0; i < points.Len(); i += 4 {
			x := points.Index(i + 0)
			y := points.Index(i + 1)
			z := points.Index(i + 2)
			if x < llx || (x == llx && z < llz) {
				llx = x
				lly = y
				llz = z
				lli = i / 4
			}
			if x > urx || (x == urx && z > urz) {
				urx = x
				ury = y
				urz = z
				uri = i / 4
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
		maxd := float64(0)
		maxi := -1
		var ci, cinc, endi int

		// Traverse the segment in lexilogical order so that the
		// max deviation is calculated similarly when traversing
		// opposite segments.
		if bx > ax || (bx == ax && bz > az) {
			cinc = 1
			ci = (ai + cinc) % pn
			endi = bi
		} else {
			cinc = pn - 1
			ci = (bi + cinc) % pn
			endi = ai
			ax, bx = bx, ax
			az, bz = bz, az
		}

		// Tessellate only outer edges or edges between areas.
		if (points.Index(ci*4+3)&RC_CONTOUR_REG_MASK) == 0 || (points.Index(ci*4+3)&RC_AREA_BORDER > 0) {
			for ci != endi {
				d := contourdistancePtSeg(points.Index(ci*4+0), points.Index(ci*4+2), ax, az, bx, bz)
				if d > maxd {
					maxd = d
					maxi = ci
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
			simplified.SetByIndex((i+1)*4+0, points.Index(maxi*4+0))
			simplified.SetByIndex((i+1)*4+1, points.Index(maxi*4+1))
			simplified.SetByIndex((i+1)*4+2, points.Index(maxi*4+2))
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
			ci := (ai + 1) % pn

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
						n = bi + pn - ai
					}
					if n > 1 {
						if bx > ax || (bx == ax && bz > az) {
							maxi = (ai + n/2) % pn
						} else {
							maxi = (ai + (n+1)/2) % pn
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
				simplified.SetByIndex((i+1)*4+3, maxi)
			} else {
				i++
			}
		}
	}

	for i := 0; i < simplified.Len()/4; i++ {
		// The edge vertex flag is take from the current raw point,
		// and the neighbour region is take from the next raw point.
		ai := (simplified.Index(i*4+3) + 1) % pn
		bi := simplified.Index(i*4 + 3)
		v := (points.Index(ai*4+3) & (RC_CONTOUR_REG_MASK | RC_AREA_BORDER)) | (points.Index(bi*4+3) & RC_BORDER_VERTEX)
		simplified.SetByIndex(i*4+3, v)
	}

}

func calcAreaOfPolygon2D(verts []int, nverts int) int {
	area := 0
	i := 0
	j := nverts - 1
	for i < nverts {
		vi := rcGetVert4(verts, i)
		vj := rcGetVert4(verts, j)
		area += vi[0]*vj[2] - vj[0]*vi[2]
		j = i
		i++
	}
	return (area + 1) / 2
}

func removeDegenerateSegments(simplified Stack[int]) {
	// Remove adjacent vertices which are equal on xz-plane,
	// or else the triangulator will get confused.
	npts := simplified.Len() / 4
	for i := 0; i < npts; i++ {
		ni := next(i, npts)
		a := []int{simplified.Index(i * 4), simplified.Index(i*4 + 1)}
		b := []int{simplified.Index(ni * 4), simplified.Index(ni*4 + 1)}
		if vequal(a, b) {
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
func mergeContours(ca, cb *rcContour, ia, ib int) bool {
	maxVerts := ca.nverts + cb.nverts + 2
	verts := make([]int, maxVerts*4)

	nv := 0

	// Copy contour A.
	for i := 0; i <= ca.nverts; i++ {
		dst := rcGetVert4(verts, nv)
		src := rcGetVert4(ca.verts, ((ia + i) % ca.nverts))
		dst[0] = src[0]
		dst[1] = src[1]
		dst[2] = src[2]
		dst[3] = src[3]
		nv++
	}

	// Copy contour B
	for i := 0; i <= cb.nverts; i++ {
		dst := rcGetVert4(verts, nv)
		src := rcGetVert4(ca.verts, ((ib + i) % cb.nverts))
		dst[0] = src[0]
		dst[1] = src[1]
		dst[2] = src[2]
		dst[3] = src[3]
		nv++
	}

	ca.verts = verts
	ca.nverts = nv

	cb.verts = []int{}
	cb.nverts = 0

	return true
}

type rcContourHole struct {
	contour              *rcContour
	minx, minz, leftmost int
}

type rcContourRegion struct {
	outline *rcContour
	holes   []*rcContourHole
	nholes  int
}

type rcPotentialDiagonal struct {
	vert int
	dist int
}

// Finds the lowest leftmost vertex of a contour.
func findLeftMostVertex(contour *rcContour, minx, minz *int, leftmost *int) {
	*minx = contour.verts[0]
	*minz = contour.verts[2]
	*leftmost = 0
	for i := 1; i < contour.nverts; i++ {
		x := contour.verts[i*4+0]
		z := contour.verts[i*4+2]
		if x < *minx || (x == *minx && z < *minz) {
			*minx = x
			*minz = z
			*leftmost = i
		}
	}
}

func compareHoles(va, vb *rcContourHole) int {
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

func mergeRegionHoles(region *rcContourRegion) {
	// Sort holes from left to right.
	for i := 0; i < region.nholes; i++ {
		findLeftMostVertex(region.holes[i].contour, &region.holes[i].minx, &region.holes[i].minz, &region.holes[i].leftmost)
	}
	sort.Slice(region.holes, func(i, j int) bool {
		return compareHoles(region.holes[i], region.holes[j]) > 0
	})
	maxVerts := region.outline.nverts
	for i := 0; i < region.nholes; i++ {
		maxVerts += region.holes[i].contour.nverts
	}

	diags := make([]*rcPotentialDiagonal, maxVerts)

	outline := region.outline

	// Merge holes into the outline one by one.
	for i := 0; i < region.nholes; i++ {
		hole := region.holes[i].contour

		index := -1
		bestVertex := region.holes[i].leftmost
		for iter := 0; iter < hole.nverts; iter++ {
			// Find potential diagonals.
			// The 'best' vertex must be in the cone described by 3 consecutive vertices of the outline.
			// ..o j-1
			//   |
			//   |   * best
			//   |
			// j o-----o j+1
			//         :
			ndiags := 0
			corner := rcGetVert4(hole.verts, bestVertex)
			for j := 0; j < outline.nverts; j++ {
				if contourInCone(j, outline.nverts, outline.verts, corner) {
					dx := outline.verts[j*4+0] - corner[0]
					dz := outline.verts[j*4+2] - corner[2]
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
				pt := rcGetVert4(outline.verts, diags[j].vert)
				intersect := intersectSegContour(pt, corner, diags[i].vert, outline.nverts, outline.verts)
				for k := i; k < region.nholes && !intersect; k++ {
					if intersect || intersectSegContour(pt, corner, -1, region.holes[k].contour.nverts, region.holes[k].contour.verts) {
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
			bestVertex = (bestVertex + 1) % hole.nverts
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
// / See the #rcConfig documentation for more information on the configuration parameters.
// /
// / @see rcAllocContourSet, rcCompactHeightfield, rcContourSet, rcConfig
func rcBuildContours(chf *rcCompactHeightfield,
	maxError float64, maxEdgeLen int,
	cset *rcContourSet, buildFlags int) bool {

	w := chf.width
	h := chf.height
	borderSize := chf.borderSize

	copy(cset.bmin, chf.bmin[:])
	copy(cset.bmax, chf.bmax[:])
	if borderSize > 0 {
		// If the heightfield was build with bordersize, remove the offset.
		pad := float64(borderSize) * chf.cs
		cset.bmin[0] += pad
		cset.bmin[2] += pad
		cset.bmax[0] -= pad
		cset.bmax[2] -= pad
	}
	cset.cs = chf.cs
	cset.ch = chf.ch
	cset.width = chf.width - chf.borderSize*2
	cset.height = chf.height - chf.borderSize*2
	cset.borderSize = chf.borderSize
	cset.maxError = maxError

	maxContours := rcMax(chf.maxRegions, 8)

	cset.conts = make([]*rcContour, maxContours)
	cset.nconts = 0
	flags := make([]int, chf.spanCount)
	// Mark boundaries.
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			c := chf.cells[x+y*w]
			i := c.index
			ni := (c.index + c.count)
			for ; i < ni; i++ {
				res := 0
				s := chf.spans[i]
				if chf.spans[i].reg == 0 || (chf.spans[i].reg&RC_BORDER_REG > 0) {
					flags[i] = 0
					continue
				}
				for dir := 0; dir < 4; dir++ {
					r := 0
					if rcGetCon(s, dir) != RC_NOT_CONNECTED {
						ax := x + rcGetDirOffsetX(dir)
						ay := y + rcGetDirOffsetY(dir)
						ai := chf.cells[ax+ay*w].index + rcGetCon(s, dir)
						r = chf.spans[ai].reg
					}
					if r == chf.spans[i].reg {
						res |= (1 << dir)
					}

				}
				flags[i] = res ^ 0xf // Inverse, mark non connected edges.
			}
		}
	}

	verts := NewStackArray(func() int { return 0 }, 256)
	simplified := NewStackArray(func() int { return 0 }, 64)

	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			c := chf.cells[x+y*w]
			i := c.index
			ni := (int)(c.index + c.count)
			for ; i < ni; i++ {
				if flags[i] == 0 || flags[i] == 0xf {
					flags[i] = 0
					continue
				}
				reg := chf.spans[i].reg
				if reg == 0 || (reg&RC_BORDER_REG > 0) {
					continue
				}

				area := chf.areas[i]

				verts.Clear()
				simplified.Clear()

				walkContour(x, y, i, chf, flags, verts)

				simplifyContour(verts, simplified, maxError, maxEdgeLen, buildFlags)
				removeDegenerateSegments(simplified)

				// Store region->contour remap info.
				// Create contour.
				if simplified.Len()/4 >= 3 {
					if cset.nconts >= maxContours {
						// Allocate more contours.
						// This happens when a region has holes.
						maxContours *= 2
						newConts := make([]*rcContour, maxContours)
						for j := 0; j < cset.nconts; j++ {
							newConts[j] = cset.conts[j]
							// Reset source pointers to prevent data deletion.
							cset.conts[j].verts = []int{}
							cset.conts[j].rverts = []int{}
						}
						cset.conts = newConts
					}

					cont := cset.conts[cset.nconts]
					cset.nconts++

					cont.nverts = simplified.Len() / 4
					cont.verts = make([]int, cont.nverts*4)
					copy(cont.verts, simplified.Slice(0, cont.nverts*4))
					if borderSize > 0 {
						// If the heightfield was build with bordersize, remove the offset.
						for j := 0; j < cont.nverts; j++ {
							v := rcGetVert4(cont.verts, j)
							v[0] -= borderSize
							v[2] -= borderSize
						}
					}

					cont.nrverts = verts.Len() / 4
					cont.rverts = make([]int, cont.nrverts*4)

					copy(cont.rverts, verts.Slice(0, cont.nrverts*4))
					if borderSize > 0 {
						// If the heightfield was build with bordersize, remove the offset.
						for j := 0; j < cont.nrverts; j++ {
							v := rcGetVert4(cont.rverts, j)
							v[0] -= borderSize
							v[2] -= borderSize
						}
					}

					cont.reg = reg
					cont.area = area
				}
			}
		}
	}

	// Merge holes if needed.
	if cset.nconts > 0 {
		// Calculate winding of all polygons.
		winding := make([]int, cset.nconts)
		nholes := 0
		for i := 0; i < cset.nconts; i++ {
			cont := cset.conts[i]
			// If the contour is wound backwards, it is a hole.
			winding[i] = 1
			if calcAreaOfPolygon2D(cont.verts, cont.nverts) < 0 {
				winding[i] = -1
			}
			if winding[i] < 0 {
				nholes++
			}

		}

		if nholes > 0 {
			// Collect outline contour and holes contours per region.
			// We assume that there is one outline and multiple holes.
			nregions := chf.maxRegions + 1
			regions := make([]*rcContourRegion, nregions)
			for i := range regions {
				regions[i] = &rcContourRegion{}
			}
			holes := make([]*rcContourHole, cset.nconts)
			for i := range holes {
				holes[i] = &rcContourHole{}
			}
			for i := 0; i < cset.nconts; i++ {
				cont := cset.conts[i]
				// Positively would contours are outlines, negative holes.
				if winding[i] > 0 {
					if regions[cont.reg].outline != nil {
						regions[cont.reg].outline = cont
					}

				} else {
					regions[cont.reg].nholes++
				}
			}
			index := 0
			for i := 0; i < nregions; i++ {
				if regions[i].nholes > 0 {
					regions[i].holes = holes[index:]
					index += regions[i].nholes
					regions[i].nholes = 0
				}
			}
			for i := 0; i < cset.nconts; i++ {
				cont := cset.conts[i]
				reg := regions[cont.reg]
				if winding[i] < 0 {
					reg.holes[reg.nholes].contour = cont
					reg.nholes++
				}

			}

			// Finally merge each regions holes into the outline.
			for i := 0; i < nregions; i++ {
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
