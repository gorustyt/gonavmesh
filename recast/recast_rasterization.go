package recast

import (
	"github.com/gorustyt/gonavmesh/common"
	"log"
	"math"
)

// / Check whether two bounding boxes overlap
// /
// / @param[in]	aMin	Min axis extents of bounding box A
// / @param[in]	aMax	Max axis extents of bounding box A
// / @param[in]	bMin	Min axis extents of bounding box B
// / @param[in]	bMax	Max axis extents of bounding box B
// / @returns true if the two bounding boxes overlap.  False otherwise.
func overlapBounds(aMin, aMax, bMin, bMax []float32) bool {
	return aMin[0] <= bMax[0] && aMax[0] >= bMin[0] &&
		aMin[1] <= bMax[1] && aMax[1] >= bMin[1] &&
		aMin[2] <= bMax[2] && aMax[2] >= bMin[2]
}

// /	Rasterize a single triangle to the heightfield.
// /
// /	This code is extremely hot, so much care should be given to maintaining maximum perf here.
// /
// / @param[in] 	v0					Triangle vertex 0
// / @param[in] 	v1					Triangle vertex 1
// / @param[in] 	v2					Triangle vertex 2
// / @param[in] 	areaID				The area ID to assign to the rasterized spans
// / @param[in] 	heightfield			Heightfield to rasterize into
// / @param[in] 	heightfieldBBMin	The min extents of the heightfield bounding box
// / @param[in] 	heightfieldBBMax	The max extents of the heightfield bounding box
// / @param[in] 	cellSize			The x and z axis size of a voxel in the heightfield
// / @param[in] 	inverseCellSize		1 / cellSize
// / @param[in] 	inverseCellHeight	1 / cellHeight
// / @param[in] 	flagMergeThreshold	The threshold in which area flags will be merged
// / @returns true if the operation completes successfully.  false if there was an error adding spans to the heightfield.
func rasterizeTri(v0, v1, v2 []float32,
	areaID uint8, heightfield *RcHeightfield,
	heightfieldBBMin, heightfieldBBMax []float32,
	cellSize, inverseCellSize, inverseCellHeight float32,
	flagMergeThreshold int32) bool {
	// Calculate the bounding box of the triangle.
	triBBMin := make([]float32, 3)
	copy(triBBMin, v0)
	common.Vmin(triBBMin, v1)
	common.Vmin(triBBMin, v2)

	triBBMax := make([]float32, 3)
	copy(triBBMax, v0)
	common.Vmax(triBBMax, v1)
	common.Vmax(triBBMax, v2)

	// If the triangle does not touch the bounding box of the heightfield, skip the triangle.
	if !overlapBounds(triBBMin, triBBMax, heightfieldBBMin, heightfieldBBMax) {
		return true
	}

	w := heightfield.Width
	h := heightfield.Height
	by := heightfieldBBMax[1] - heightfieldBBMin[1]

	// Calculate the footprint of the triangle on the grid's z-axis
	z0 := int32((triBBMin[2] - heightfieldBBMin[2]) * inverseCellSize)
	z1 := int32((triBBMax[2] - heightfieldBBMin[2]) * inverseCellSize)

	// use -1 rather than 0 to cut the polygon properly at the start of the tile
	z0 = common.Clamp(z0, -1, h-1)
	z1 = common.Clamp(z1, 0, h-1)

	// Clip the triangle into all grid cells it touches.
	buf := make([]float32, 7*3*4)
	in := buf
	inRow := buf[7*3:]
	p1 := inRow[7*3:]
	p2 := p1[7*3:]

	copy(in[0:], v0)
	copy(in[1*3:], v1)
	copy(in[2*3:], v2)
	var nvRow int32
	nvIn := int32(3)

	for z := z0; z <= z1; z++ {
		// Clip polygon to row. Store the remaining polygon as well
		cellZ := heightfieldBBMin[2] + float32(z)*cellSize
		dividePoly(in, nvIn, inRow, &nvRow, p1, &nvIn, cellZ+cellSize, RC_AXIS_Z)
		in, p1 = p1, in

		if nvRow < 3 {
			continue
		}
		if z < 0 {
			continue
		}

		// find X-axis bounds of the row
		minX := inRow[0]
		maxX := inRow[0]
		for vert := int32(1); vert < nvRow; vert++ {
			if minX > inRow[vert*3] {
				minX = inRow[vert*3]
			}
			if maxX < inRow[vert*3] {
				maxX = inRow[vert*3]
			}
		}
		x0 := int32((minX - heightfieldBBMin[0]) * inverseCellSize)
		x1 := int32((maxX - heightfieldBBMin[0]) * inverseCellSize)
		if x1 < 0 || x0 >= w {
			continue
		}
		x0 = common.Clamp(x0, -1, w-1)
		x1 = common.Clamp(x1, 0, w-1)

		var nv int32
		nv2 := nvRow

		for x := x0; x <= x1; x++ {
			// Clip polygon to column. store the remaining polygon as well
			cx := heightfieldBBMin[0] + float32(x)*cellSize
			dividePoly(inRow, nv2, p1, &nv, p2, &nv2, cx+cellSize, RC_AXIS_X)
			inRow, p2 = p2, inRow

			if nv < 3 {
				continue
			}
			if x < 0 {
				continue
			}

			// Calculate min and max of the span.
			spanMin := p1[1]
			spanMax := p1[1]
			for vert := int32(1); vert < nv; vert++ {
				spanMin = min(spanMin, p1[vert*3+1])
				spanMax = max(spanMax, p1[vert*3+1])
			}
			spanMin -= heightfieldBBMin[1]
			spanMax -= heightfieldBBMin[1]

			// Skip the span if it's completely outside the heightfield bounding box
			if spanMax < 0.0 {
				continue
			}
			if spanMin > by {
				continue
			}

			// Clamp the span to the heightfield bounding box.
			if spanMin < 0.0 {
				spanMin = 0
			}
			if spanMax > by {
				spanMax = by
			}

			// Snap the span to the heightfield height grid.
			spanMinCellIndex := common.Clamp(uint16(math.Floor(float64(spanMin*inverseCellHeight))), 0, RC_SPAN_MAX_HEIGHT)
			spanMaxCellIndex := common.Clamp(uint16(math.Ceil(float64(spanMax*inverseCellHeight))), spanMinCellIndex+1, RC_SPAN_MAX_HEIGHT)

			if !addSpan(heightfield, x, z, spanMinCellIndex, spanMaxCellIndex, areaID, flagMergeThreshold) {
				return false
			}
		}
	}

	return true
}
func RcRasterizeTriangles(verts []float32, nv int32, tris []int32, triAreaIDs []uint8, numTris int32,
	heightfield *RcHeightfield, flagMergeThreshold int32) bool {

	// Rasterize the triangles.
	inverseCellSize := 1.0 / heightfield.Cs
	inverseCellHeight := 1.0 / heightfield.Ch
	for triIndex := int32(0); triIndex < numTris; triIndex++ {
		v0 := common.GetVert3(verts, tris[triIndex*3+0])
		v1 := common.GetVert3(verts, tris[triIndex*3+1])
		v2 := common.GetVert3(verts, tris[triIndex*3+2])
		if !rasterizeTri(v0, v1, v2, triAreaIDs[triIndex], heightfield, heightfield.Bmin[:], heightfield.Bmax[:], heightfield.Cs, inverseCellSize, inverseCellHeight, flagMergeThreshold) {
			return false
		}
	}

	return true
}

func RcRasterizeTriangles1(verts []float32, triAreaIDs []uint8, numTris int,
	heightfield *RcHeightfield, flagMergeThreshold int32) bool {
	// Rasterize the triangles.
	inverseCellSize := 1.0 / heightfield.Cs
	inverseCellHeight := 1.0 / heightfield.Ch
	for triIndex := 0; triIndex < numTris; triIndex++ {
		v0 := common.GetVert3(verts, triIndex*3+0)
		v1 := common.GetVert3(verts, triIndex*3+1)
		v2 := common.GetVert3(verts, triIndex*3+2)
		if !rasterizeTri(v0, v1, v2, triAreaIDs[triIndex], heightfield, heightfield.Bmin[:], heightfield.Bmax[:], heightfield.Cs, inverseCellSize, inverseCellHeight, flagMergeThreshold) {
			log.Printf("rcRasterizeTriangles: Out of memory.")
			return false
		}
	}

	return true
}

type rcAxis int

const (
	RC_AXIS_X = 0
	RC_AXIS_Y = 1
	RC_AXIS_Z = 2
)

// / Divides a convex polygon of max 12 vertices into two convex polygons
// / across a separating axis.
// /
// / @param[in]	inVerts			The input polygon vertices
// / @param[in]	inVertsCount	The number of input polygon vertices
// / @param[out]	outVerts1		Resulting polygon 1's vertices
// / @param[out]	outVerts1Count	The number of resulting polygon 1 vertices
// / @param[out]	outVerts2		Resulting polygon 2's vertices
// / @param[out]	outVerts2Count	The number of resulting polygon 2 vertices
// / @param[in]	axisOffset		THe offset along the specified axis
// / @param[in]	axis			The separating axis
func dividePoly(inVerts []float32, inVertsCount int32,
	outVerts1 []float32, outVerts1Count *int32,
	outVerts2 []float32, outVerts2Count *int32,
	axisOffset float32, axis rcAxis) {
	common.AssertTrue(inVertsCount <= 12)

	// How far positive or negative away from the separating axis is each vertex.
	inVertAxisDelta := make([]float32, 12)
	for inVert := int32(0); inVert < inVertsCount; inVert++ {
		inVertAxisDelta[inVert] = axisOffset - inVerts[inVert*3+int32(axis)]
	}

	poly1Vert := int32(0)
	poly2Vert := int32(0)
	inVertA := int32(0)
	inVertB := inVertsCount - 1
	for inVertA < inVertsCount {
		// If the two vertices are on the same side of the separating axis
		sameSide := (inVertAxisDelta[inVertA] >= 0) == (inVertAxisDelta[inVertB] >= 0)

		if !sameSide {
			s := inVertAxisDelta[inVertB] / (inVertAxisDelta[inVertB] - inVertAxisDelta[inVertA])
			outVerts1[poly1Vert*3+0] = inVerts[inVertB*3+0] + (inVerts[inVertA*3+0]-inVerts[inVertB*3+0])*s
			outVerts1[poly1Vert*3+1] = inVerts[inVertB*3+1] + (inVerts[inVertA*3+1]-inVerts[inVertB*3+1])*s
			outVerts1[poly1Vert*3+2] = inVerts[inVertB*3+2] + (inVerts[inVertA*3+2]-inVerts[inVertB*3+2])*s

			copy(common.GetVert3(outVerts2, poly2Vert), common.GetVert3(outVerts1, poly1Vert))
			poly1Vert++
			poly2Vert++

			// add the inVertA point to the right polygon. Do NOT add points that are on the dividing line
			// since these were already added above
			if inVertAxisDelta[inVertA] > 0 {

				copy(common.GetVert3(outVerts1, poly1Vert), common.GetVert3(inVerts, inVertA))
				poly1Vert++
			} else if inVertAxisDelta[inVertA] < 0 {
				copy(common.GetVert3(outVerts1, poly2Vert), common.GetVert3(inVerts, inVertA))
				poly2Vert++
			}
		} else {
			// add the inVertA point to the right polygon. Addition is done even for points on the dividing line
			if inVertAxisDelta[inVertA] >= 0 {
				copy(common.GetVert3(outVerts1, poly1Vert), common.GetVert3(inVerts, inVertA))
				poly1Vert++
				if inVertAxisDelta[inVertA] != 0 {
					inVertB = inVertA
					inVertA++
					continue
				}
			}

			copy(common.GetVert3(outVerts2, poly2Vert), common.GetVert3(inVerts, inVertA))
			poly2Vert++
		}
		inVertB = inVertA
		inVertA++
	}

	*outVerts1Count = poly1Vert
	*outVerts2Count = poly2Vert
}

// / Adds a span to the heightfield.  If the new span overlaps existing spans,
// / it will merge the new span with the existing ones.
// /
// / @param[in]	heightfield					Heightfield to add spans to
// / @param[in]	x					The new span's column cell x index
// / @param[in]	z					The new span's column cell z index
// / @param[in]	min					The new span's minimum cell index
// / @param[in]	max					The new span's maximum cell index
// / @param[in]	areaID				The new span's area type ID
// / @param[in]	flagMergeThreshold	How close two spans maximum extents need to be to merge area type IDs
func addSpan(heightfield *RcHeightfield,
	x, z int32,
	minValue, maxValue uint16, areaID uint8, flagMergeThreshold int32) bool {
	// Create the new span.
	newSpan := allocSpan(heightfield)
	if newSpan == nil {
		return false
	}
	newSpan.Smin = uint32(minValue)
	newSpan.Smax = uint32(maxValue)
	newSpan.Area = uint32(areaID)
	newSpan.Next = nil

	columnIndex := x + z*heightfield.Width
	var previousSpan *RcSpan
	currentSpan := heightfield.Spans[columnIndex]

	// Insert the new span, possibly merging it with existing spans.
	for currentSpan != nil {
		if currentSpan.Smin > newSpan.Smax {
			// Current span is completely after the new span, break.
			break
		}

		if currentSpan.Smax < newSpan.Smin {
			// Current span is completely before the new span.  Keep going.
			previousSpan = currentSpan
			currentSpan = currentSpan.Next
		} else {
			// The new span overlaps with an existing span.  Merge them.
			if currentSpan.Smin < newSpan.Smin {
				newSpan.Smin = currentSpan.Smin
			}
			if currentSpan.Smax > newSpan.Smax {
				newSpan.Smax = currentSpan.Smax
			}

			// Merge flags.
			if common.Abs(float64(newSpan.Smax-currentSpan.Smax)) <= float64(flagMergeThreshold) {
				// Higher area ID numbers indicate higher resolution priority.
				newSpan.Area = max(newSpan.Area, currentSpan.Area)
			}

			// Remove the current span since it's now merged with newSpan.
			// Keep going because there might be other overlapping spans that also need to be merged.
			next := currentSpan.Next
			freeSpan(heightfield, currentSpan)
			if previousSpan != nil {
				previousSpan.Next = next
			} else {
				heightfield.Spans[columnIndex] = next
			}
			currentSpan = next
		}
	}

	// Insert new span after prev
	if previousSpan != nil {
		newSpan.Next = previousSpan.Next
		previousSpan.Next = newSpan
	} else {
		// This span should go before the others in the list
		newSpan.Next = heightfield.Spans[columnIndex]
		heightfield.Spans[columnIndex] = newSpan
	}

	return true
}

// / Releases the memory used by the span back to the heightfield, so it can be re-used for new spans.
// / @param[in]	heightfield		The heightfield.
// / @param[in]	span	A pointer to the span to free
func freeSpan(heightfield *RcHeightfield, span *RcSpan) {
	if span == nil {
		return
	}
	// Add the span to the front of the free list.
	span.Next = heightfield.Freelist
	heightfield.Freelist = span
}

// / Allocates a new span in the heightfield.
// / Use a memory pool and free list to minimize actual allocations.
// /
// / @param[in]	heightfield		The heightfield
// / @returns A pointer to the allocated or re-used span memory.
func allocSpan(heightfield *RcHeightfield) *RcSpan {
	// If necessary, allocate new page and update the freelist.
	if heightfield.Freelist == nil || heightfield.Freelist.Next == nil {
		// Create new page.
		// Allocate memory for the new pool.
		spanPool := &RcSpanPool{}

		// Add the pool into the list of pools.
		spanPool.next = heightfield.Pools
		heightfield.Pools = spanPool

		// Add new spans to the free list.
		freeList := heightfield.Freelist
		head := &spanPool.items[0]
		it := RC_SPANS_PER_POOL - 1
		it--
		spanPool.items[it].Next = freeList
		freeList = &spanPool.items[it]
		for &spanPool.items[it] != head {
			it--
			spanPool.items[it].Next = freeList
			freeList = &spanPool.items[it]
		}

		heightfield.Freelist = &spanPool.items[it]
	}

	// Pop item from the front of the free list.
	newSpan := heightfield.Freelist
	heightfield.Freelist = heightfield.Freelist.Next
	return newSpan
}

func rcAddSpan(heightfield *RcHeightfield,
	x, z int32, spanMin, spanMax uint16, areaID uint8, flagMergeThreshold int32) bool {

	if !addSpan(heightfield, x, z, spanMin, spanMax, areaID, flagMergeThreshold) {
		return false
	}

	return true
}

func rcRasterizeTriangle(
	v0, v1, v2 []float32,
	areaID uint8, heightfield *RcHeightfield, flagMergeThreshold int32) bool {

	// Rasterize the single triangle.
	inverseCellSize := 1.0 / heightfield.Cs
	inverseCellHeight := 1.0 / heightfield.Ch
	if !rasterizeTri(v0, v1, v2, areaID, heightfield, heightfield.Bmin[:], heightfield.Bmax[:], heightfield.Cs, inverseCellSize, inverseCellHeight, flagMergeThreshold) {
		return false
	}
	return true
}
