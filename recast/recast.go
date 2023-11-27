package recast

import (
	"gonavamesh/common"
	"log"
	"math"
)

// / Specifies a configuration to use when performing Recast builds.
// / @ingroup recast
type RcConfig struct {
	/// The width of the field along the x-axis. [Limit: >= 0] [Units: vx]
	Width int

	/// The height of the field along the z-axis. [Limit: >= 0] [Units: vx]
	Height int

	/// The width/height size of tile's on the xz-plane. [Limit: >= 0] [Units: vx]
	TileSize int

	/// The size of the non-navigable border around the heightfield. [Limit: >=0] [Units: vx]
	BorderSize int

	/// The xz-plane cell size to use for fields. [Limit: > 0] [Units: wu]
	Cs float64

	/// The y-axis cell size to use for fields. [Limit: > 0] [Units: wu]
	Ch float64

	/// The minimum bounds of the field's AABB. [(x, y, z)] [Units: wu]
	Bmin [3]float64

	/// The maximum bounds of the field's AABB. [(x, y, z)] [Units: wu]
	Bmax [3]float64

	/// The maximum slope that is considered walkable. [Limits: 0 <= value < 90] [Units: Degrees]
	WalkableSlopeAngle float64

	/// Minimum floor to 'ceiling' height that will still allow the floor area to
	/// be considered walkable. [Limit: >= 3] [Units: vx]
	WalkableHeight int

	/// Maximum ledge height that is considered to still be traversable. [Limit: >=0] [Units: vx]
	WalkableClimb int

	/// The distance to erode/shrink the walkable area of the heightfield away from
	/// obstructions.  [Limit: >=0] [Units: vx]
	WalkableRadius int

	/// The maximum allowed length for contour edges along the border of the mesh. [Limit: >=0] [Units: vx]
	MaxEdgeLen int

	/// The maximum distance a simplified contour's border edges should deviate
	/// the original raw contour. [Limit: >=0] [Units: vx]
	MaxSimplificationError float64

	/// The minimum number of cells allowed to form isolated island areas. [Limit: >=0] [Units: vx]
	MinRegionArea int

	/// Any regions with a span count smaller than this value will, if possible,
	/// be merged with larger regions. [Limit: >=0] [Units: vx]
	MergeRegionArea int

	/// The maximum number of vertices allowed for polygons generated during the
	/// contour to polygon conversion process. [Limit: >= 3]
	MaxVertsPerPoly int

	/// Sets the sampling distance to use when generating the detail mesh.
	/// (For height detail only.) [Limits: 0 or >= 0.9] [Units: wu]
	DetailSampleDist float64

	/// The maximum distance the detail mesh surface should deviate from heightfield
	/// data. (For height detail only.) [Limit: >=0] [Units: wu]
	DetailSampleMaxError float64
}

// / The default area id used to indicate a walkable polygon.
// / This is also the maximum allowed area id, and the only non-null area id
// / recognized by some steps in the build process.
const RC_WALKABLE_AREA = 63

func rcCalcBounds(verts []float64, numVerts int, minBounds []float64, maxBounds []float64) {
	// Calculate bounding box.
	copy(minBounds, verts)
	copy(maxBounds, verts)
	for i := 1; i < numVerts; i++ {
		v := rcGetVert(verts, i)
		common.Vmin(minBounds, v)
		common.Vmax(maxBounds, v)
	}
}

func RcCalcGridSize(minBounds, maxBounds []float64, cellSize float64, sizeX, sizeZ *int) {
	*sizeX = int((maxBounds[0]-minBounds[0])/cellSize + 0.5)
	*sizeZ = int((maxBounds[2]-minBounds[2])/cellSize + 0.5)
}

func calcTriNormal(v0, v1, v2 []float64, faceNormal []float64) {
	e0 := make([]float64, 3)
	e1 := make([]float64, 3)
	common.Vsub(e0, v1, v0)
	common.Vsub(e1, v2, v0)
	common.Vcross(faceNormal, e0, e1)
	common.Vnormalize(faceNormal)
}

func RcMarkWalkableTriangles(walkableSlopeAngle float64,
	verts []float64, numVerts int,
	tris []int, numTris int,
	triAreaIDs []int) {

	walkableThr := math.Cos(walkableSlopeAngle / 180.0 * math.Pi)

	norm := make([]float64, 3)

	for i := 0; i < numTris; i++ {
		tri := rcGetVert(tris, i)

		calcTriNormal(rcGetVert(verts, tri[0]), rcGetVert(verts, tri[1]), rcGetVert(verts, tri[2]), norm)
		// Check if the face is walkable.
		if norm[1] > walkableThr {
			triAreaIDs[i] = RC_WALKABLE_AREA
		}
	}
}

func rcClearUnwalkableTriangles(walkableSlopeAngle float64,
	verts []float64, numVerts int,
	tris []int, numTris int,
	triAreaIDs []int) {
	// The minimum Y value for a face normal of a triangle with a walkable slope.
	walkableLimitY := math.Cos(walkableSlopeAngle / 180.0 * math.Pi)

	faceNormal := make([]float64, 3)
	for i := 0; i < numTris; i++ {
		tri := rcGetVert(tris, i)
		calcTriNormal(rcGetVert(verts, tri[0]), rcGetVert(verts, tri[1]), rcGetVert(verts, tri[2]), faceNormal)
		// Check if the face is walkable.
		if faceNormal[1] <= walkableLimitY {
			triAreaIDs[i] = RC_NULL_AREA
		}
	}
}

func rcGetHeightFieldSpanCount(heightfield *RcHeightfield) int {
	numCols := heightfield.Width * heightfield.Height
	spanCount := 0
	for columnIndex := 0; columnIndex < numCols; columnIndex++ {
		for span := heightfield.Spans[columnIndex]; span != nil; span = span.Next {
			if span.Area != RC_NULL_AREA {
				spanCount++
			}
		}
	}
	return spanCount
}

func RcBuildCompactHeightfield(walkableHeight, walkableClimb int,
	heightfield *RcHeightfield, compactHeightfield *RcCompactHeightfield) bool {
	xSize := heightfield.Width
	zSize := heightfield.Height
	spanCount := rcGetHeightFieldSpanCount(heightfield)

	// Fill in header.
	compactHeightfield.Width = xSize
	compactHeightfield.Height = zSize
	compactHeightfield.SpanCount = spanCount
	compactHeightfield.WalkableHeight = walkableHeight
	compactHeightfield.WalkableClimb = walkableClimb
	compactHeightfield.MaxRegions = 0
	copy(compactHeightfield.Bmin[:], heightfield.Bmin[:])
	copy(compactHeightfield.Bmax[:], heightfield.Bmax[:])
	compactHeightfield.Bmax[1] += float64(walkableHeight) * heightfield.Ch
	compactHeightfield.Cs = heightfield.Cs
	compactHeightfield.Ch = heightfield.Ch
	compactHeightfield.Cells = make([]*RcCompactCell, xSize*zSize)
	for i := range compactHeightfield.Cells {
		compactHeightfield.Cells[i] = &RcCompactCell{}
	}
	compactHeightfield.Spans = make([]*RcCompactSpan, spanCount)
	for i := range compactHeightfield.Spans {
		compactHeightfield.Spans[i] = &RcCompactSpan{}
	}
	compactHeightfield.Areas = make([]int, spanCount)
	for i := range compactHeightfield.Areas {
		compactHeightfield.Areas[i] = RC_NULL_AREA
	}
	MAX_HEIGHT := 0xffff

	// Fill in cells and spans.
	currentCellIndex := 0
	numColumns := xSize * zSize
	for columnIndex := 0; columnIndex < numColumns; columnIndex++ {
		span := heightfield.Spans[columnIndex]

		// If there are no spans at this cell, just leave the data to index=0, count=0.
		if span == nil {
			continue
		}

		cell := compactHeightfield.Cells[columnIndex]
		cell.Index = currentCellIndex
		cell.Count = 0

		for ; span != nil; span = span.Next {
			if span.Area != RC_NULL_AREA {
				bot := span.Smax
				top := MAX_HEIGHT
				if span.Next != nil {
					top = span.Next.Smin
				}
				compactHeightfield.Spans[currentCellIndex].Y = common.Clamp(bot, 0, 0xffff)
				compactHeightfield.Spans[currentCellIndex].H = common.Clamp(top-bot, 0, 0xff)
				compactHeightfield.Areas[currentCellIndex] = span.Area
				currentCellIndex++
				cell.Count++
			}
		}
	}

	// Find neighbour connections.
	MAX_LAYERS := RC_NOT_CONNECTED - 1
	maxLayerIndex := 0
	zStride := xSize // for readability
	for z := 0; z < zSize; z++ {
		for x := 0; x < xSize; x++ {
			cell := compactHeightfield.Cells[x+z*zStride]
			i := cell.Index
			ni := (cell.Index + cell.Count)
			for ; i < ni; i++ {
				span := compactHeightfield.Spans[i]

				for dir := 0; dir < 4; dir++ {
					rcSetCon(span, dir, RC_NOT_CONNECTED)
					neighborX := x + common.GetDirOffsetX(dir)
					neighborZ := z + common.GetDirOffsetY(dir)
					// First check that the neighbour cell is in bounds.
					if neighborX < 0 || neighborZ < 0 || neighborX >= xSize || neighborZ >= zSize {
						continue
					}

					// Iterate over all neighbour spans and check if any of the is
					// accessible from current cell.
					neighborCell := compactHeightfield.Cells[neighborX+neighborZ*zStride]
					k := neighborCell.Index
					nk := (neighborCell.Index + neighborCell.Count)
					for ; k < nk; k++ {
						neighborSpan := compactHeightfield.Spans[k]
						bot := common.Max(span.Y, neighborSpan.Y)
						top := common.Min(span.Y+span.H, neighborSpan.Y+neighborSpan.H)

						// Check that the gap between the spans is walkable,
						// and that the climb height between the gaps is not too high.
						if (top-bot) >= walkableHeight && common.Abs(neighborSpan.Y-span.Y) <= walkableClimb {
							// Mark direction as walkable.
							layerIndex := k - neighborCell.Index
							if layerIndex < 0 || layerIndex > MAX_LAYERS {
								maxLayerIndex = common.Max(maxLayerIndex, layerIndex)
								continue
							}
							rcSetCon(span, dir, layerIndex)
							break
						}
					}
				}
			}
		}
	}

	if maxLayerIndex > MAX_LAYERS {
		log.Printf("rcBuildCompactHeightfield: Heightfield has too many layers %d (max: %d)", maxLayerIndex, MAX_LAYERS)
	}

	return true
}

func RcCreateHeightfield(sizeX, sizeZ int,
	minBounds, maxBounds []float64,
	cellSize, cellHeight float64) (heightfield *RcHeightfield) {
	heightfield = &RcHeightfield{}
	heightfield.Width = sizeX
	heightfield.Height = sizeZ
	copy(heightfield.Bmin[:], minBounds)
	copy(heightfield.Bmax[:], maxBounds)
	heightfield.Cs = cellSize
	heightfield.Ch = cellHeight
	heightfield.Spans = make([]*RcSpan, heightfield.Width*heightfield.Height)
	return heightfield
}
