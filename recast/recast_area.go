package recast

import (
	"gonavamesh/common"
	"sort"
)

const (
	RC_NOT_CONNECTED = 0x3f
)

type rcCompactCell struct {
	index int ///< Index to the first span in the column.
	count int ///< Number of spans in the column.
}

func newRcCompactCell() *rcCompactCell {
	return &rcCompactCell{
		index: 24,
		count: 8,
	}
}

// / Represents a span of unobstructed space within a compact heightfield.
type rcCompactSpan struct {
	y   int ///< The lower extent of the span. (Measured from the heightfield's base.)
	reg int ///< The id of the region the span belongs to. (Or zero if not in a region.)
	con int ///< Packed neighbor connection data.
	h   int ///< The height of the span.  (Measured from #y.)
}

func newRcCompactSpan() *rcCompactSpan {
	return &rcCompactSpan{
		h:   8,
		con: 24,
	}
}

func rcGetCon(span *rcCompactSpan, direction int) int {
	shift := direction * 6
	return (span.con >> shift) & 0x3f
}

// / A compact, static heightfield representing unobstructed space.
// / @ingroup recast
type RcCompactHeightfield struct {
	width          int              ///< The width of the heightfield. (Along the x-axis in cell units.)
	height         int              ///< The height of the heightfield. (Along the z-axis in cell units.)
	spanCount      int              ///< The number of spans in the heightfield.
	walkableHeight int              ///< The walkable height used during the build of the field.  (See: RcConfig::walkableHeight)
	walkableClimb  int              ///< The walkable climb used during the build of the field. (See: RcConfig::walkableClimb)
	borderSize     int              ///< The AABB border size used during the build of the field. (See: RcConfig::borderSize)
	maxDistance    int              ///< The maximum distance value of any span within the field.
	maxRegions     int              ///< The maximum region id of any span within the field.
	bmin           [3]float64       ///< The minimum bounds in world space. [(x, y, z)]
	bmax           [3]float64       ///< The maximum bounds in world space. [(x, y, z)]
	cs             float64          ///< The size of each cell. (On the xz-plane.)
	ch             float64          ///< The height of each cell. (The minimum increment along the y-axis.)
	cells          []*rcCompactCell ///< Array of cells. [Size: #width*#height]
	spans          []*rcCompactSpan ///< Array of spans. [Size: #spanCount]
	dist           []int            ///< Array containing border distance data. [Size: #spanCount]
	areas          []int            ///< Array containing area id data. [Size: #spanCount]
}

// TODO (graham): This is duplicated in the ConvexVolumeTool in RecastDemo
// / Checks if a point is contained within a polygon
// /
// / @param[in]	numVerts	Number of vertices in the polygon
// / @param[in]	verts		The polygon vertices
// / @param[in]	point		The point to check
// / @returns true if the point lies within the polygon, false otherwise.
func pointInPoly(numVerts int, verts []float64, point []float64) bool {
	inPoly := false
	i := 0
	j := numVerts - 1
	for i < numVerts {
		vi := verts[i*3 : i*3+3]
		vj := verts[j*3 : j*3+3]

		if (vi[2] > point[2]) == (vj[2] > point[2]) {
			j = i
			i++
			continue
		}

		if point[0] >= (vj[0]-vi[0])*(point[2]-vi[2])/(vj[2]-vi[2])+vi[0] {
			j = i
			i++
			continue
		}
		inPoly = !inPoly
		j = i
		i++
	}
	return inPoly
}

func RcErodeWalkableArea(erosionRadius int, compactHeightfield *RcCompactHeightfield) bool {

	xSize := compactHeightfield.width
	zSize := compactHeightfield.height
	zStride := xSize // For readability

	distanceToBoundary := make([]int, compactHeightfield.spanCount)
	for i := range distanceToBoundary {
		distanceToBoundary[i] = 0xff
	}
	// Mark boundary cells.
	for z := 0; z < zSize; z++ {
		for x := 0; x < xSize; x++ {
			cell := compactHeightfield.cells[x+z*zStride]
			spanIndex := cell.index
			maxSpanIndex := cell.index + cell.count
			for ; spanIndex < maxSpanIndex; spanIndex++ {
				if compactHeightfield.areas[spanIndex] == RC_NULL_AREA {
					distanceToBoundary[spanIndex] = 0
					continue
				}
				span := compactHeightfield.spans[spanIndex]

				// Check that there is a non-null adjacent span in each of the 4 cardinal directions.
				neighborCount := 0
				for direction := 0; direction < 4; direction++ {
					neighborConnection := rcGetCon(span, direction)
					if neighborConnection == RC_NOT_CONNECTED {
						break
					}

					neighborX := x + common.GetDirOffsetX(direction)
					neighborZ := z + common.GetDirOffsetY(direction)
					neighborSpanIndex := compactHeightfield.cells[neighborX+neighborZ*zStride].index + neighborConnection

					if compactHeightfield.areas[neighborSpanIndex] == RC_NULL_AREA {
						break
					}
					neighborCount++
				}

				// At least one missing neighbour, so this is a boundary cell.
				if neighborCount != 4 {
					distanceToBoundary[spanIndex] = 0
				}
			}
		}
	}
	// Pass 1
	for z := 0; z < zSize; z++ {
		for x := 0; x < xSize; x++ {
			cell := compactHeightfield.cells[x+z*zStride]
			maxSpanIndex := cell.index + cell.count
			for spanIndex := cell.index; spanIndex < maxSpanIndex; spanIndex++ {
				span := compactHeightfield.spans[spanIndex]

				if rcGetCon(span, 0) != RC_NOT_CONNECTED {
					// (-1,0)
					aX := x + common.GetDirOffsetX(0)
					aY := z + common.GetDirOffsetY(0)
					aIndex := compactHeightfield.cells[aX+aY*xSize].index + rcGetCon(span, 0)
					aSpan := compactHeightfield.spans[aIndex]
					newDistance := common.Min(distanceToBoundary[aIndex]+2, 255)
					if newDistance < distanceToBoundary[spanIndex] {
						distanceToBoundary[spanIndex] = newDistance
					}

					// (-1,-1)
					if rcGetCon(aSpan, 3) != RC_NOT_CONNECTED {
						bX := aX + common.GetDirOffsetX(3)
						bY := aY + common.GetDirOffsetY(3)
						bIndex := compactHeightfield.cells[bX+bY*xSize].index + rcGetCon(aSpan, 3)
						newDistance = common.Min(distanceToBoundary[bIndex]+3, 255)
						if newDistance < distanceToBoundary[spanIndex] {
							distanceToBoundary[spanIndex] = newDistance
						}
					}
				}
				if rcGetCon(span, 3) != RC_NOT_CONNECTED {
					// (0,-1)
					aX := x + common.GetDirOffsetX(3)
					aY := z + common.GetDirOffsetY(3)
					aIndex := compactHeightfield.cells[aX+aY*xSize].index + rcGetCon(span, 3)
					aSpan := compactHeightfield.spans[aIndex]
					newDistance := common.Min(distanceToBoundary[aIndex]+2, 255)
					if newDistance < distanceToBoundary[spanIndex] {
						distanceToBoundary[spanIndex] = newDistance
					}

					// (1,-1)
					if rcGetCon(aSpan, 2) != RC_NOT_CONNECTED {
						bX := aX + common.GetDirOffsetX(2)
						bY := aY + common.GetDirOffsetY(2)
						bIndex := compactHeightfield.cells[bX+bY*xSize].index + rcGetCon(aSpan, 2)
						newDistance := common.Min(distanceToBoundary[bIndex]+3, 255)
						if newDistance < distanceToBoundary[spanIndex] {
							distanceToBoundary[spanIndex] = newDistance
						}
					}
				}
			}
		}
	}

	// Pass 2
	for z := zSize - 1; z >= 0; z-- {
		for x := xSize - 1; x >= 0; x-- {
			cell := compactHeightfield.cells[x+z*zStride]
			maxSpanIndex := cell.index + cell.count
			for spanIndex := cell.index; spanIndex < maxSpanIndex; spanIndex++ {
				span := compactHeightfield.spans[spanIndex]

				if rcGetCon(span, 2) != RC_NOT_CONNECTED {
					// (1,0)
					aX := x + common.GetDirOffsetX(2)
					aY := z + common.GetDirOffsetY(2)
					aIndex := compactHeightfield.cells[aX+aY*xSize].index + rcGetCon(span, 2)
					aSpan := compactHeightfield.spans[aIndex]
					newDistance := common.Min(distanceToBoundary[aIndex]+2, 255)
					if newDistance < distanceToBoundary[spanIndex] {
						distanceToBoundary[spanIndex] = newDistance
					}

					// (1,1)
					if rcGetCon(aSpan, 1) != RC_NOT_CONNECTED {
						bX := aX + common.GetDirOffsetX(1)
						bY := aY + common.GetDirOffsetY(1)
						bIndex := compactHeightfield.cells[bX+bY*xSize].index + rcGetCon(aSpan, 1)
						newDistance = common.Min(distanceToBoundary[bIndex]+3, 255)
						if newDistance < distanceToBoundary[spanIndex] {
							distanceToBoundary[spanIndex] = newDistance
						}
					}
				}
				if rcGetCon(span, 1) != RC_NOT_CONNECTED {
					// (0,1)
					aX := x + common.GetDirOffsetX(1)
					aY := z + common.GetDirOffsetY(1)
					aIndex := compactHeightfield.cells[aX+aY*xSize].index + rcGetCon(span, 1)
					aSpan := compactHeightfield.spans[aIndex]
					newDistance := common.Min(distanceToBoundary[aIndex]+2, 255)
					if newDistance < distanceToBoundary[spanIndex] {
						distanceToBoundary[spanIndex] = newDistance
					}

					// (-1,1)
					if rcGetCon(aSpan, 0) != RC_NOT_CONNECTED {
						bX := aX + common.GetDirOffsetX(0)
						bY := aY + common.GetDirOffsetY(0)
						bIndex := compactHeightfield.cells[bX+bY*xSize].index + rcGetCon(aSpan, 0)
						newDistance := common.Min(distanceToBoundary[bIndex]+3, 255)
						if newDistance < distanceToBoundary[spanIndex] {
							distanceToBoundary[spanIndex] = newDistance
						}
					}
				}
			}
		}
	}

	minBoundaryDistance := erosionRadius * 2
	for spanIndex := 0; spanIndex < compactHeightfield.spanCount; spanIndex++ {
		if distanceToBoundary[spanIndex] < minBoundaryDistance {
			compactHeightfield.areas[spanIndex] = RC_NULL_AREA
		}
	}

	return true
}

func rcMedianFilterWalkableArea(compactHeightfield *RcCompactHeightfield) bool {

	xSize := compactHeightfield.width
	zSize := compactHeightfield.height
	zStride := xSize // For readability

	areas := make([]int, compactHeightfield.spanCount)

	for i := range areas {
		areas[i] = 0xff
	}

	for z := 0; z < zSize; z++ {
		for x := 0; x < xSize; x++ {
			cell := compactHeightfield.cells[x+z*zStride]
			maxSpanIndex := cell.index + cell.count
			for spanIndex := cell.index; spanIndex < maxSpanIndex; spanIndex++ {
				span := compactHeightfield.spans[spanIndex]
				if compactHeightfield.areas[spanIndex] == RC_NULL_AREA {
					areas[spanIndex] = compactHeightfield.areas[spanIndex]
					continue
				}

				var neighborAreas [9]int
				for neighborIndex := 0; neighborIndex < 9; neighborIndex++ {
					neighborAreas[neighborIndex] = compactHeightfield.areas[spanIndex]
				}

				for dir := 0; dir < 4; dir++ {
					if rcGetCon(span, dir) == RC_NOT_CONNECTED {
						continue
					}

					aX := x + common.GetDirOffsetX(dir)
					aZ := z + common.GetDirOffsetY(dir)
					aIndex := compactHeightfield.cells[aX+aZ*zStride].index + rcGetCon(span, dir)
					if compactHeightfield.areas[aIndex] != RC_NULL_AREA {
						neighborAreas[dir*2+0] = compactHeightfield.areas[aIndex]
					}

					aSpan := compactHeightfield.spans[aIndex]
					dir2 := (dir + 1) & 0x3
					neighborConnection2 := rcGetCon(aSpan, dir2)
					if neighborConnection2 != RC_NOT_CONNECTED {
						bX := aX + common.GetDirOffsetX(dir2)
						bZ := aZ + common.GetDirOffsetY(dir2)
						bIndex := compactHeightfield.cells[bX+bZ*zStride].index + neighborConnection2
						if compactHeightfield.areas[bIndex] != RC_NULL_AREA {
							neighborAreas[dir*2+1] = compactHeightfield.areas[bIndex]
						}
					}
				}
				//recast 这里用的插入排序
				var tmp = neighborAreas[:]
				sort.Slice(tmp, func(i, j int) bool {
					return tmp[i] < tmp[j]
				})
				for i := range tmp {
					neighborAreas[i] = tmp[i]
				}
				areas[spanIndex] = neighborAreas[4]
			}
		}
	}
	compactHeightfield.areas = areas
	return true
}

func rcMarkBoxArea(boxMinBounds, boxMaxBounds []float64, areaId int, compactHeightfield *RcCompactHeightfield) {

	xSize := compactHeightfield.width
	zSize := compactHeightfield.height
	zStride := xSize // For readability

	// Find the footprint of the box area in grid cell coordinates.
	minX := int((boxMinBounds[0] - compactHeightfield.bmin[0]) / compactHeightfield.cs)
	minY := int((boxMinBounds[1] - compactHeightfield.bmin[1]) / compactHeightfield.ch)
	minZ := int((boxMinBounds[2] - compactHeightfield.bmin[2]) / compactHeightfield.cs)
	maxX := int((boxMaxBounds[0] - compactHeightfield.bmin[0]) / compactHeightfield.cs)
	maxY := int((boxMaxBounds[1] - compactHeightfield.bmin[1]) / compactHeightfield.ch)
	maxZ := int((boxMaxBounds[2] - compactHeightfield.bmin[2]) / compactHeightfield.cs)

	// Early-out if the box is outside the bounds of the grid.
	if maxX < 0 {
		return
	}
	if minX >= xSize {
		return
	}
	if maxZ < 0 {
		return
	}
	if minZ >= zSize {
		return
	}

	// Clamp relevant bound coordinates to the grid.
	if minX < 0 {
		minX = 0
	}
	if maxX >= xSize {
		maxX = xSize - 1
	}
	if minZ < 0 {
		minZ = 0
	}
	if maxZ >= zSize {
		maxZ = zSize - 1
	}

	// Mark relevant cells.
	for z := minZ; z <= maxZ; z++ {
		for x := minX; x <= maxX; x++ {
			cell := compactHeightfield.cells[x+z*zStride]
			maxSpanIndex := cell.index + cell.count
			for spanIndex := cell.index; spanIndex < maxSpanIndex; spanIndex++ {
				span := compactHeightfield.spans[spanIndex]

				// Skip if the span is outside the box extents.
				if span.y < minY || span.y > maxY {
					continue
				}

				// Skip if the span has been removed.
				if compactHeightfield.areas[spanIndex] == RC_NULL_AREA {
					continue
				}

				// Mark the span.
				compactHeightfield.areas[spanIndex] = areaId
			}
		}
	}
}

func RcMarkConvexPolyArea(verts []float64, numVerts int, minY, maxY float64, areaId int, compactHeightfield *RcCompactHeightfield) {

	xSize := compactHeightfield.width
	zSize := compactHeightfield.height
	zStride := xSize // For readability

	// Compute the bounding box of the polygon
	bmin := make([]float64, 3)
	bmax := make([]float64, 3)
	copy(bmin, verts)
	copy(bmax, verts)
	for i := 1; i < numVerts; i++ {
		common.Vmin(bmin, verts[i*3:i*3+3])
		common.Vmax(bmax, verts[i*3:i*3+3])
	}
	bmin[1] = minY
	bmax[1] = maxY

	// Compute the grid footprint of the polygon
	minx := int((bmin[0] - compactHeightfield.bmin[0]) / compactHeightfield.cs)
	miny := int((bmin[1] - compactHeightfield.bmin[1]) / compactHeightfield.ch)
	minz := int((bmin[2] - compactHeightfield.bmin[2]) / compactHeightfield.cs)
	maxx := int((bmax[0] - compactHeightfield.bmin[0]) / compactHeightfield.cs)
	maxy := int((bmax[1] - compactHeightfield.bmin[1]) / compactHeightfield.ch)
	maxz := int((bmax[2] - compactHeightfield.bmin[2]) / compactHeightfield.cs)

	// Early-out if the polygon lies entirely outside the grid.
	if maxx < 0 {
		return
	}
	if minx >= xSize {
		return
	}
	if maxz < 0 {
		return
	}
	if minz >= zSize {
		return
	}

	// Clamp the polygon footprint to the grid
	if minx < 0 {
		minx = 0
	}
	if maxx >= xSize {
		maxx = xSize - 1
	}
	if minz < 0 {
		minz = 0
	}
	if maxz >= zSize {
		maxz = zSize - 1
	}

	// TODO: Optimize.
	for z := minz; z <= maxz; z++ {
		for x := minx; x <= maxx; x++ {
			cell := compactHeightfield.cells[x+z*zStride]
			maxSpanIndex := (int)(cell.index + cell.count)
			for spanIndex := cell.index; spanIndex < maxSpanIndex; spanIndex++ {
				span := compactHeightfield.spans[spanIndex]

				// Skip if span is removed.
				if compactHeightfield.areas[spanIndex] == RC_NULL_AREA {
					continue
				}

				// Skip if y extents don't overlap.
				if span.y < miny || span.y > maxy {
					continue
				}

				point := []float64{
					compactHeightfield.bmin[0] + (float64(x)+0.5)*compactHeightfield.cs,
					0,
					compactHeightfield.bmin[2] + (float64(z)+0.5)*compactHeightfield.cs,
				}

				if pointInPoly(numVerts, verts, point) {
					compactHeightfield.areas[spanIndex] = areaId
				}
			}
		}
	}
}

const (
	EPSILON = 1e-6
)

// / Normalizes the vector if the length is greater than zero.
// / If the magnitude is zero, the vector is unchanged.
// / @param[in,out]	v	The vector to normalize. [(x, y, z)]
func rcVsafeNormalize(v []float64) {
	sqMag := common.Sqr(v[0]) + common.Sqr(v[1]) + common.Sqr(v[2])
	if sqMag > EPSILON {
		inverseMag := 1.0 / common.Sqrt(sqMag)
		v[0] *= inverseMag
		v[1] *= inverseMag
		v[2] *= inverseMag
	}
}

func RcOffsetPoly(verts []float64, numVerts int, offset float64, outVerts []float64, maxOutVerts int) int {
	// Defines the limit at which a miter becomes a bevel.
	// Similar in behavior to https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/stroke-miterlimit
	const MITER_LIMIT float64 = 1.20
	numOutVerts := 0

	for vertIndex := 0; vertIndex < numVerts; vertIndex++ {
		// Grab three vertices of the polygon.
		vertIndexA := (vertIndex + numVerts - 1) % numVerts
		vertIndexB := vertIndex
		vertIndexC := (vertIndex + 1) % numVerts
		vertA := verts[vertIndexA*3 : vertIndexA*3+3]
		vertB := verts[vertIndexB*3 : vertIndexB*3+3]
		vertC := verts[vertIndexC*3 : vertIndexC*3+3]

		// From A to B on the x/z plane

		prevSegmentDir := common.Vsub(vertB, vertA)
		prevSegmentDir[1] = 0 // Squash onto x/z plane
		rcVsafeNormalize(prevSegmentDir)

		// From B to C on the x/z plane

		currSegmentDir := common.Vsub(vertC, vertB)
		currSegmentDir[1] = 0 // Squash onto x/z plane
		rcVsafeNormalize(currSegmentDir)

		// The y component of the cross product of the two normalized segment directions.
		// The X and Z components of the cross product are both zero because the two
		// segment direction vectors fall within the x/z plane.
		cross := currSegmentDir[0]*prevSegmentDir[2] - prevSegmentDir[0]*currSegmentDir[2]

		// CCW perpendicular vector to AB.  The segment normal.
		prevSegmentNormX := -prevSegmentDir[2]
		prevSegmentNormZ := prevSegmentDir[0]

		// CCW perpendicular vector to BC.  The segment normal.
		currSegmentNormX := -currSegmentDir[2]
		currSegmentNormZ := currSegmentDir[0]

		// Average the two segment normals to get the proportional miter offset for B.
		// This isn't normalized because it's defining the distance and direction the corner will need to be
		// adjusted proportionally to the edge offsets to properly miter the adjoining edges.
		cornerMiterX := (prevSegmentNormX + currSegmentNormX) * 0.5
		cornerMiterZ := (prevSegmentNormZ + currSegmentNormZ) * 0.5
		cornerMiterSqMag := common.Sqr(cornerMiterX) + common.Sqr(cornerMiterZ)

		// If the magnitude of the segment normal average is less than about .69444,
		// the corner is an acute enough angle that the result should be beveled.
		bevel := cornerMiterSqMag*MITER_LIMIT*MITER_LIMIT < 1.0

		// Scale the corner miter so it's proportional to how much the corner should be offset compared to the edges.
		if cornerMiterSqMag > EPSILON {
			scale := 1.0 / cornerMiterSqMag
			cornerMiterX *= scale
			cornerMiterZ *= scale
		}

		if bevel && cross < 0.0 { // If the corner is convex and an acute enough angle, generate a bevel.

			if numOutVerts+2 > maxOutVerts {
				return 0
			}

			// Generate two bevel vertices at a distances from B proportional to the angle between the two segments.
			// Move each bevel vertex out proportional to the given offset.
			d := 1.0 - (prevSegmentDir[0]*currSegmentDir[0]+prevSegmentDir[2]*currSegmentDir[2])*0.5

			outVerts[numOutVerts*3+0] = vertB[0] + (-prevSegmentNormX+prevSegmentDir[0]*d)*offset
			outVerts[numOutVerts*3+1] = vertB[1]
			outVerts[numOutVerts*3+2] = vertB[2] + (-prevSegmentNormZ+prevSegmentDir[2]*d)*offset
			numOutVerts++

			outVerts[numOutVerts*3+0] = vertB[0] + (-currSegmentNormX-currSegmentDir[0]*d)*offset
			outVerts[numOutVerts*3+1] = vertB[1]
			outVerts[numOutVerts*3+2] = vertB[2] + (-currSegmentNormZ-currSegmentDir[2]*d)*offset
			numOutVerts++
		} else {
			if numOutVerts+1 > maxOutVerts {
				return 0
			}

			// Move B along the miter direction by the specified offset.
			outVerts[numOutVerts*3+0] = vertB[0] - cornerMiterX*offset
			outVerts[numOutVerts*3+1] = vertB[1]
			outVerts[numOutVerts*3+2] = vertB[2] - cornerMiterZ*offset
			numOutVerts++
		}
	}

	return numOutVerts
}

func rcMarkCylinderArea(position []float64, radius float64, height float64, areaId int, compactHeightfield *RcCompactHeightfield) {

	xSize := compactHeightfield.width
	zSize := compactHeightfield.height
	zStride := xSize // For readability

	// Compute the bounding box of the cylinder
	cylinderBBMin := []float64{
		position[0] - radius,
		position[1],
		position[2] - radius}
	cylinderBBMax := []float64{
		position[0] + radius,
		position[1] + height,
		position[2] + radius}

	// Compute the grid footprint of the cylinder
	minx := int((cylinderBBMin[0] - compactHeightfield.bmin[0]) / compactHeightfield.cs)
	miny := int((cylinderBBMin[1] - compactHeightfield.bmin[1]) / compactHeightfield.ch)
	minz := int((cylinderBBMin[2] - compactHeightfield.bmin[2]) / compactHeightfield.cs)
	maxx := int((cylinderBBMax[0] - compactHeightfield.bmin[0]) / compactHeightfield.cs)
	maxy := int((cylinderBBMax[1] - compactHeightfield.bmin[1]) / compactHeightfield.ch)
	maxz := int((cylinderBBMax[2] - compactHeightfield.bmin[2]) / compactHeightfield.cs)

	// Early-out if the cylinder is completely outside the grid bounds.
	if maxx < 0 {
		return
	}
	if minx >= xSize {
		return
	}
	if maxz < 0 {
		return
	}
	if minz >= zSize {
		return
	}

	// Clamp the cylinder bounds to the grid.
	if minx < 0 {
		minx = 0
	}
	if maxx >= xSize {
		maxx = xSize - 1
	}
	if minz < 0 {
		minz = 0
	}
	if maxz >= zSize {
		maxz = zSize - 1
	}

	radiusSq := radius * radius

	for z := minz; z <= maxz; z++ {
		for x := minx; x <= maxx; x++ {
			cell := compactHeightfield.cells[x+z*zStride]
			maxSpanIndex := cell.index + cell.count

			cellX := compactHeightfield.bmin[0] + (float64(x)+0.5)*compactHeightfield.cs
			cellZ := compactHeightfield.bmin[2] + (float64(z)+0.5)*compactHeightfield.cs
			deltaX := cellX - position[0]
			deltaZ := cellZ - position[2]

			// Skip this column if it's too far from the center point of the cylinder.
			if common.Sqr(deltaX)+common.Sqr(deltaZ) >= radiusSq {
				continue
			}

			// Mark all overlapping spans
			for spanIndex := cell.index; spanIndex < maxSpanIndex; spanIndex++ {
				span := compactHeightfield.spans[spanIndex]

				// Skip if span is removed.
				if compactHeightfield.areas[spanIndex] == RC_NULL_AREA {
					continue
				}

				// Mark if y extents overlap.
				if span.y >= miny && span.y <= maxy {
					compactHeightfield.areas[spanIndex] = areaId
				}
			}
		}
	}
}
