package recast

import (
	"github.com/gorustyt/gonavmesh/common"
	"sort"
)

const (
	RC_NOT_CONNECTED = 0x3f
)

type RcCompactCell struct {
	Index uint32 ///< Index to the first span in the column.
	Count uint32 ///< Number of spans in the column.
}

func newRcCompactCell() *RcCompactCell {
	return &RcCompactCell{
		Index: 24,
		Count: 8,
	}
}

// / Represents a span of unobstructed space within a compact heightfield.
type RcCompactSpan struct {
	Y   uint16 ///< The lower extent of the span. (Measured from the heightfield's base.)
	Reg uint16 ///< The id of the region the span belongs to. (Or zero if not in a region.)
	Con int32  ///< Packed neighbor connection data.
	H   int32  ///< The height of the span.  (Measured from #y.)
}

func newRcCompactSpan() *RcCompactSpan {
	return &RcCompactSpan{
		H:   8,
		Con: 24,
	}
}

func rcGetCon(span *RcCompactSpan, direction int32) int32 {
	shift := direction * 6
	return (span.Con >> shift) & 0x3f
}

// / A compact, static heightfield representing unobstructed space.
// / @ingroup recast
type RcCompactHeightfield struct {
	Width          int32            ///< The width of the heightfield. (Along the x-axis in cell units.)
	Height         int32            ///< The height of the heightfield. (Along the z-axis in cell units.)
	SpanCount      int32            ///< The number of spans in the heightfield.
	WalkableHeight int32            ///< The walkable height used during the build of the field.  (See: RcConfig::walkableHeight)
	WalkableClimb  int32            ///< The walkable climb used during the build of the field. (See: RcConfig::walkableClimb)
	BorderSize     int32            ///< The AABB border size used during the build of the field. (See: RcConfig::borderSize)
	MaxDistance    uint16           ///< The maximum distance value of any span within the field.
	MaxRegions     uint16           ///< The maximum region id of any span within the field.
	Bmin           [3]float32       ///< The minimum bounds in world space. [(x, y, z)]
	Bmax           [3]float32       ///< The maximum bounds in world space. [(x, y, z)]
	Cs             float32          ///< The size of each cell. (On the xz-plane.)
	Ch             float32          ///< The height of each cell. (The minimum increment along the y-axis.)
	Cells          []*RcCompactCell ///< Array of cells. [Size: #width*#height]
	Spans          []*RcCompactSpan ///< Array of spans. [Size: #spanCount]
	Dist           []uint16         ///< Array containing border distance data. [Size: #spanCount]
	Areas          []uint8          ///< Array containing area id data. [Size: #spanCount]
}

func RcErodeWalkableArea(erosionRadius int, compactHeightfield *RcCompactHeightfield) bool {

	xSize := compactHeightfield.Width
	zSize := compactHeightfield.Height
	zStride := xSize // For readability

	distanceToBoundary := make([]int, compactHeightfield.SpanCount)
	for i := range distanceToBoundary {
		distanceToBoundary[i] = 0xff
	}
	// Mark boundary cells.
	for z := int32(0); z < zSize; z++ {
		for x := int32(0); x < xSize; x++ {
			cell := compactHeightfield.Cells[x+z*zStride]
			spanIndex := cell.Index
			maxSpanIndex := cell.Index + cell.Count
			for ; spanIndex < maxSpanIndex; spanIndex++ {
				if compactHeightfield.Areas[spanIndex] == RC_NULL_AREA {
					distanceToBoundary[spanIndex] = 0
					continue
				}
				span := compactHeightfield.Spans[spanIndex]

				// Check that there is a non-null adjacent span in each of the 4 cardinal directions.
				neighborCount := 0
				for direction := int32(0); direction < 4; direction++ {
					neighborConnection := rcGetCon(span, direction)
					if neighborConnection == RC_NOT_CONNECTED {
						break
					}

					neighborX := x + common.GetDirOffsetX(direction)
					neighborZ := z + common.GetDirOffsetY(direction)
					neighborSpanIndex := int32(compactHeightfield.Cells[neighborX+neighborZ*zStride].Index) + neighborConnection

					if compactHeightfield.Areas[neighborSpanIndex] == RC_NULL_AREA {
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
	for z := int32(0); z < zSize; z++ {
		for x := int32(0); x < xSize; x++ {
			cell := compactHeightfield.Cells[x+z*zStride]
			maxSpanIndex := cell.Index + cell.Count
			for spanIndex := cell.Index; spanIndex < maxSpanIndex; spanIndex++ {
				span := compactHeightfield.Spans[spanIndex]

				if rcGetCon(span, 0) != RC_NOT_CONNECTED {
					// (-1,0)
					aX := x + common.GetDirOffsetX(0)
					aY := z + common.GetDirOffsetY(0)
					aIndex := int32(compactHeightfield.Cells[aX+aY*xSize].Index) + rcGetCon(span, 0)
					aSpan := compactHeightfield.Spans[aIndex]
					newDistance := common.Min(distanceToBoundary[aIndex]+2, 255)
					if newDistance < distanceToBoundary[spanIndex] {
						distanceToBoundary[spanIndex] = newDistance
					}

					// (-1,-1)
					if rcGetCon(aSpan, 3) != RC_NOT_CONNECTED {
						bX := aX + common.GetDirOffsetX(3)
						bY := aY + common.GetDirOffsetY(3)
						bIndex := int32(compactHeightfield.Cells[bX+bY*xSize].Index) + rcGetCon(aSpan, 3)
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
					aIndex := int32(compactHeightfield.Cells[aX+aY*xSize].Index) + rcGetCon(span, 3)
					aSpan := compactHeightfield.Spans[aIndex]
					newDistance := common.Min(distanceToBoundary[aIndex]+2, 255)
					if newDistance < distanceToBoundary[spanIndex] {
						distanceToBoundary[spanIndex] = newDistance
					}

					// (1,-1)
					if rcGetCon(aSpan, 2) != RC_NOT_CONNECTED {
						bX := aX + common.GetDirOffsetX(2)
						bY := aY + common.GetDirOffsetY(2)
						bIndex := int32(compactHeightfield.Cells[bX+bY*xSize].Index) + rcGetCon(aSpan, 2)
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
			cell := compactHeightfield.Cells[x+z*zStride]
			maxSpanIndex := cell.Index + cell.Count
			for spanIndex := cell.Index; spanIndex < maxSpanIndex; spanIndex++ {
				span := compactHeightfield.Spans[spanIndex]

				if rcGetCon(span, 2) != RC_NOT_CONNECTED {
					// (1,0)
					aX := x + common.GetDirOffsetX(2)
					aY := z + common.GetDirOffsetY(2)
					aIndex := int32(compactHeightfield.Cells[aX+aY*xSize].Index) + rcGetCon(span, 2)
					aSpan := compactHeightfield.Spans[aIndex]
					newDistance := common.Min(distanceToBoundary[aIndex]+2, 255)
					if newDistance < distanceToBoundary[spanIndex] {
						distanceToBoundary[spanIndex] = newDistance
					}

					// (1,1)
					if rcGetCon(aSpan, 1) != RC_NOT_CONNECTED {
						bX := aX + common.GetDirOffsetX(1)
						bY := aY + common.GetDirOffsetY(1)
						bIndex := int32(compactHeightfield.Cells[bX+bY*xSize].Index) + rcGetCon(aSpan, 1)
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
					aIndex := int32(compactHeightfield.Cells[aX+aY*xSize].Index) + rcGetCon(span, 1)
					aSpan := compactHeightfield.Spans[aIndex]
					newDistance := common.Min(distanceToBoundary[aIndex]+2, 255)
					if newDistance < distanceToBoundary[spanIndex] {
						distanceToBoundary[spanIndex] = newDistance
					}

					// (-1,1)
					if rcGetCon(aSpan, 0) != RC_NOT_CONNECTED {
						bX := aX + common.GetDirOffsetX(0)
						bY := aY + common.GetDirOffsetY(0)
						bIndex := int32(compactHeightfield.Cells[bX+bY*xSize].Index) + rcGetCon(aSpan, 0)
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
	for spanIndex := int32(0); spanIndex < compactHeightfield.SpanCount; spanIndex++ {
		if distanceToBoundary[spanIndex] < minBoundaryDistance {
			compactHeightfield.Areas[spanIndex] = RC_NULL_AREA
		}
	}

	return true
}

func rcMedianFilterWalkableArea(compactHeightfield *RcCompactHeightfield) bool {

	xSize := compactHeightfield.Width
	zSize := compactHeightfield.Height
	zStride := xSize // For readability

	areas := make([]uint8, compactHeightfield.SpanCount)

	for i := range areas {
		areas[i] = 0xff
	}

	for z := int32(0); z < zSize; z++ {
		for x := int32(0); x < xSize; x++ {
			cell := compactHeightfield.Cells[x+z*zStride]
			maxSpanIndex := cell.Index + cell.Count
			for spanIndex := cell.Index; spanIndex < maxSpanIndex; spanIndex++ {
				span := compactHeightfield.Spans[spanIndex]
				if compactHeightfield.Areas[spanIndex] == RC_NULL_AREA {
					areas[spanIndex] = compactHeightfield.Areas[spanIndex]
					continue
				}

				var neighborAreas [9]uint8
				for neighborIndex := 0; neighborIndex < 9; neighborIndex++ {
					neighborAreas[neighborIndex] = compactHeightfield.Areas[spanIndex]
				}

				for dir := int32(0); dir < 4; dir++ {
					if rcGetCon(span, dir) == RC_NOT_CONNECTED {
						continue
					}

					aX := x + common.GetDirOffsetX(dir)
					aZ := z + common.GetDirOffsetY(dir)
					aIndex := int32(compactHeightfield.Cells[aX+aZ*zStride].Index) + rcGetCon(span, dir)
					if compactHeightfield.Areas[aIndex] != RC_NULL_AREA {
						neighborAreas[dir*2+0] = compactHeightfield.Areas[aIndex]
					}

					aSpan := compactHeightfield.Spans[aIndex]
					dir2 := (dir + 1) & 0x3
					neighborConnection2 := rcGetCon(aSpan, dir2)
					if neighborConnection2 != RC_NOT_CONNECTED {
						bX := aX + common.GetDirOffsetX(dir2)
						bZ := aZ + common.GetDirOffsetY(dir2)
						bIndex := int32(compactHeightfield.Cells[bX+bZ*zStride].Index) + neighborConnection2
						if compactHeightfield.Areas[bIndex] != RC_NULL_AREA {
							neighborAreas[dir*2+1] = compactHeightfield.Areas[bIndex]
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
	compactHeightfield.Areas = areas
	return true
}

func rcMarkBoxArea(boxMinBounds, boxMaxBounds []float32, areaId uint8, compactHeightfield *RcCompactHeightfield) {

	xSize := compactHeightfield.Width
	zSize := compactHeightfield.Height
	zStride := xSize // For readability

	// Find the footprint of the box area in grid cell coordinates.
	minX := int32((boxMinBounds[0] - compactHeightfield.Bmin[0]) / compactHeightfield.Cs)
	minY := int32((boxMinBounds[1] - compactHeightfield.Bmin[1]) / compactHeightfield.Ch)
	minZ := int32((boxMinBounds[2] - compactHeightfield.Bmin[2]) / compactHeightfield.Cs)
	maxX := int32((boxMaxBounds[0] - compactHeightfield.Bmin[0]) / compactHeightfield.Cs)
	maxY := int32((boxMaxBounds[1] - compactHeightfield.Bmin[1]) / compactHeightfield.Ch)
	maxZ := int32((boxMaxBounds[2] - compactHeightfield.Bmin[2]) / compactHeightfield.Cs)

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
			cell := compactHeightfield.Cells[x+z*zStride]
			maxSpanIndex := cell.Index + cell.Count
			for spanIndex := cell.Index; spanIndex < maxSpanIndex; spanIndex++ {
				span := compactHeightfield.Spans[spanIndex]

				// Skip if the span is outside the box extents.
				if int32(span.Y) < minY || int32(span.Y) > maxY {
					continue
				}

				// Skip if the span has been removed.
				if compactHeightfield.Areas[spanIndex] == RC_NULL_AREA {
					continue
				}

				// Mark the span.
				compactHeightfield.Areas[spanIndex] = areaId
			}
		}
	}
}

func RcMarkConvexPolyArea(verts []float32, numVerts int32, minY, maxY float32, areaId uint8, compactHeightfield *RcCompactHeightfield) {

	xSize := compactHeightfield.Width
	zSize := compactHeightfield.Height
	zStride := xSize // For readability

	// Compute the bounding box of the polygon
	bmin := make([]float32, 3)
	bmax := make([]float32, 3)
	copy(bmin, verts)
	copy(bmax, verts)
	for i := int32(1); i < numVerts; i++ {
		common.Vmin(bmin, verts[i*3:i*3+3])
		common.Vmax(bmax, verts[i*3:i*3+3])
	}
	bmin[1] = minY
	bmax[1] = maxY

	// Compute the grid footprint of the polygon
	minx := int32((bmin[0] - compactHeightfield.Bmin[0]) / compactHeightfield.Cs)
	miny := int32((bmin[1] - compactHeightfield.Bmin[1]) / compactHeightfield.Ch)
	minz := int32((bmin[2] - compactHeightfield.Bmin[2]) / compactHeightfield.Cs)
	maxx := int32((bmax[0] - compactHeightfield.Bmin[0]) / compactHeightfield.Cs)
	maxy := int32((bmax[1] - compactHeightfield.Bmin[1]) / compactHeightfield.Ch)
	maxz := int32((bmax[2] - compactHeightfield.Bmin[2]) / compactHeightfield.Cs)

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
			cell := compactHeightfield.Cells[x+z*zStride]
			maxSpanIndex := int32(cell.Index + cell.Count)
			for spanIndex := cell.Index; int(spanIndex) < int(maxSpanIndex); spanIndex++ {
				span := compactHeightfield.Spans[spanIndex]

				// Skip if span is removed.
				if compactHeightfield.Areas[spanIndex] == RC_NULL_AREA {
					continue
				}

				// Skip if y extents don't overlap.
				if int32(span.Y) < miny || int32(span.Y) > maxy {
					continue
				}

				point := []float32{
					compactHeightfield.Bmin[0] + (float32(x)+0.5)*compactHeightfield.Cs,
					0,
					compactHeightfield.Bmin[2] + (float32(z)+0.5)*compactHeightfield.Cs,
				}

				if common.PointInPoly(numVerts, verts, point) {
					compactHeightfield.Areas[spanIndex] = areaId
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
func rcVsafeNormalize(v []float32) {
	sqMag := common.Sqr(v[0]) + common.Sqr(v[1]) + common.Sqr(v[2])
	if sqMag > EPSILON {
		inverseMag := float32(1.0 / common.Sqrt(float64(sqMag)))
		v[0] *= inverseMag
		v[1] *= inverseMag
		v[2] *= inverseMag
	}
}

func RcOffsetPoly(verts []float32, numVerts int, offset float32, outVerts []float32, maxOutVerts int) int {
	// Defines the limit at which a miter becomes a bevel.
	// Similar in behavior to https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/stroke-miterlimit
	const MITER_LIMIT float32 = 1.20
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
		prevSegmentDir := make([]float32, 3)
		common.Vsub(prevSegmentDir, vertB, vertA)
		prevSegmentDir[1] = 0 // Squash onto x/z plane
		rcVsafeNormalize(prevSegmentDir)

		// From B to C on the x/z plane
		currSegmentDir := make([]float32, 3)
		common.Vsub(currSegmentDir, vertC, vertB)
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

func rcMarkCylinderArea(position []float32, radius float32, height float32, areaId uint8, compactHeightfield *RcCompactHeightfield) {

	xSize := compactHeightfield.Width
	zSize := compactHeightfield.Height
	zStride := xSize // For readability

	// Compute the bounding box of the cylinder
	cylinderBBMin := []float32{
		position[0] - radius,
		position[1],
		position[2] - radius}
	cylinderBBMax := []float32{
		position[0] + radius,
		position[1] + height,
		position[2] + radius}

	// Compute the grid footprint of the cylinder
	minx := int32((cylinderBBMin[0] - compactHeightfield.Bmin[0]) / compactHeightfield.Cs)
	miny := int32((cylinderBBMin[1] - compactHeightfield.Bmin[1]) / compactHeightfield.Ch)
	minz := int32((cylinderBBMin[2] - compactHeightfield.Bmin[2]) / compactHeightfield.Cs)
	maxx := int32((cylinderBBMax[0] - compactHeightfield.Bmin[0]) / compactHeightfield.Cs)
	maxy := int32((cylinderBBMax[1] - compactHeightfield.Bmin[1]) / compactHeightfield.Ch)
	maxz := int32((cylinderBBMax[2] - compactHeightfield.Bmin[2]) / compactHeightfield.Cs)

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
			cell := compactHeightfield.Cells[x+z*zStride]
			maxSpanIndex := cell.Index + cell.Count

			cellX := compactHeightfield.Bmin[0] + (float32(x)+0.5)*compactHeightfield.Cs
			cellZ := compactHeightfield.Bmin[2] + (float32(z)+0.5)*compactHeightfield.Cs
			deltaX := cellX - position[0]
			deltaZ := cellZ - position[2]

			// Skip this column if it's too far from the center point of the cylinder.
			if common.Sqr(deltaX)+common.Sqr(deltaZ) >= radiusSq {
				continue
			}

			// Mark all overlapping spans
			for spanIndex := cell.Index; spanIndex < maxSpanIndex; spanIndex++ {
				span := compactHeightfield.Spans[spanIndex]

				// Skip if span is removed.
				if compactHeightfield.Areas[spanIndex] == RC_NULL_AREA {
					continue
				}

				// Mark if y extents overlap.
				if int32(span.Y) >= miny && int32(span.Y) <= maxy {
					compactHeightfield.Areas[spanIndex] = areaId
				}
			}
		}
	}
}
