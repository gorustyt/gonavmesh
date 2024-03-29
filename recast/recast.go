package recast

import (
	"github.com/gorustyt/gonavmesh/common"
	"log"
	"math"
)

// / Specifies a configuration to use when performing Recast builds.
// / @ingroup recast
type RcConfig struct {
	/// The width of the field along the x-axis. [Limit: >= 0] [Units: vx]
	Width int32

	/// The height of the field along the z-axis. [Limit: >= 0] [Units: vx]
	Height int32

	/// The width/height size of tile's on the xz-plane. [Limit: >= 0] [Units: vx]
	TileSize int32

	/// The size of the non-navigable border around the heightfield. [Limit: >=0] [Units: vx]
	BorderSize int32

	/// The xz-plane cell size to use for fields. [Limit: > 0] [Units: wu]
	Cs float32

	/// The y-axis cell size to use for fields. [Limit: > 0] [Units: wu]
	Ch float32

	/// The minimum bounds of the field's AABB. [(x, y, z)] [Units: wu]
	Bmin [3]float32

	/// The maximum bounds of the field's AABB. [(x, y, z)] [Units: wu]
	Bmax [3]float32

	/// The maximum slope that is considered walkable. [Limits: 0 <= value < 90] [Units: Degrees]
	WalkableSlopeAngle float32

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
	MaxSimplificationError float32

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
	DetailSampleDist float32

	/// The maximum distance the detail mesh surface should deviate from heightfield
	/// data. (For height detail only.) [Limit: >=0] [Units: wu]
	DetailSampleMaxError float32
}

// / The default area id used to indicate a walkable polygon.
// / This is also the maximum allowed area id, and the only non-null area id
// / recognized by some steps in the build process.
const RC_WALKABLE_AREA = 63

func RcCalcGridSize(minBounds, maxBounds []float32, cellSize float32, sizeX, sizeZ *int32) {
	*sizeX = int32((maxBounds[0]-minBounds[0])/cellSize + 0.5)
	*sizeZ = int32((maxBounds[2]-minBounds[2])/cellSize + 0.5)
}

func calcTriNormal(v0, v1, v2 []float32, faceNormal []float32) {
	e0 := make([]float32, 3)
	e1 := make([]float32, 3)
	common.Vsub(e0, v1, v0)
	common.Vsub(e1, v2, v0)
	common.Vcross(faceNormal, e0, e1)
	common.Vnormalize(faceNormal)
}

func RcMarkWalkableTriangles(walkableSlopeAngle float32,
	verts []float32, numVerts int,
	tris []int, numTris int,
	triAreaIDs []int) {

	walkableThr := math.Cos(float64(walkableSlopeAngle) / 180.0 * math.Pi)

	norm := make([]float32, 3)

	for i := 0; i < numTris; i++ {
		tri := common.GetVert3(tris, i)

		calcTriNormal(common.GetVert3(verts, tri[0]), common.GetVert3(verts, tri[1]), common.GetVert3(verts, tri[2]), norm)
		// Check if the face is walkable.
		if float64(norm[1]) > walkableThr {
			triAreaIDs[i] = RC_WALKABLE_AREA
		}
	}
}

func rcClearUnwalkableTriangles(walkableSlopeAngle float32,
	verts []float32, numVerts int,
	tris []int, numTris int,
	triAreaIDs []int) {
	// The minimum Y value for a face normal of a triangle with a walkable slope.
	walkableLimitY := math.Cos(float64(walkableSlopeAngle) / 180.0 * math.Pi)

	faceNormal := make([]float32, 3)
	for i := 0; i < numTris; i++ {
		tri := common.GetVert3(tris, i)
		calcTriNormal(common.GetVert3(verts, tri[0]), common.GetVert3(verts, tri[1]), common.GetVert3(verts, tri[2]), faceNormal)
		// Check if the face is walkable.
		if float64(faceNormal[1]) <= walkableLimitY {
			triAreaIDs[i] = RC_NULL_AREA
		}
	}
}

func rcGetHeightFieldSpanCount(heightfield *RcHeightfield) int32 {
	numCols := heightfield.Width * heightfield.Height
	spanCount := int32(0)
	for columnIndex := int32(0); columnIndex < numCols; columnIndex++ {
		for span := heightfield.Spans[columnIndex]; span != nil; span = span.Next {
			if span.Area != RC_NULL_AREA {
				spanCount++
			}
		}
	}
	return spanCount
}

func RcBuildCompactHeightfield(walkableHeight, walkableClimb int32,
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
	compactHeightfield.Bmax[1] += float32(walkableHeight) * heightfield.Ch
	compactHeightfield.Cs = heightfield.Cs
	compactHeightfield.Ch = heightfield.Ch
	compactHeightfield.Cells = make([]*RcCompactCell, xSize*zSize)
	for i := range compactHeightfield.Cells {
		compactHeightfield.Cells[i] = newRcCompactCell()
	}
	compactHeightfield.Spans = make([]*RcCompactSpan, spanCount)
	for i := range compactHeightfield.Spans {
		compactHeightfield.Spans[i] = newRcCompactSpan()
	}
	compactHeightfield.Areas = make([]uint8, spanCount)
	for i := range compactHeightfield.Areas {
		compactHeightfield.Areas[i] = RC_NULL_AREA
	}
	MAX_HEIGHT := 0xffff

	// Fill in cells and spans.
	currentCellIndex := uint32(0)
	numColumns := xSize * zSize
	for columnIndex := int32(0); columnIndex < numColumns; columnIndex++ {
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
				top := int32(MAX_HEIGHT)
				if span.Next != nil {
					top = int32(span.Next.Smin)
				}
				compactHeightfield.Spans[currentCellIndex].Y = common.Clamp(uint16(bot), 0, 0xffff)
				compactHeightfield.Spans[currentCellIndex].H = common.Clamp(top-int32(bot), 0, 0xff)
				compactHeightfield.Areas[currentCellIndex] = uint8(span.Area)
				currentCellIndex++
				cell.Count++
			}
		}
	}

	// Find neighbour connections.
	MAX_LAYERS := RC_NOT_CONNECTED - 1
	maxLayerIndex := 0
	zStride := xSize // for readability
	for z := int32(0); z < zSize; z++ {
		for x := int32(0); x < xSize; x++ {
			cell := compactHeightfield.Cells[x+z*zStride]
			i := cell.Index
			ni := (cell.Index + cell.Count)
			for ; i < ni; i++ {
				span := compactHeightfield.Spans[i]

				for dir := int32(0); dir < 4; dir++ {
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
						bot := max(int32(span.Y), int32(neighborSpan.Y))
						top := min(int32(span.Y)+span.H, int32(neighborSpan.Y)+neighborSpan.H)

						// Check that the gap between the spans is walkable,
						// and that the climb height between the gaps is not too high.
						if (top-bot) >= walkableHeight && common.Abs(int32(neighborSpan.Y)-int32(span.Y)) <= walkableClimb {
							// Mark direction as walkable.
							layerIndex := k - neighborCell.Index
							if layerIndex < 0 || int(layerIndex) > MAX_LAYERS {
								maxLayerIndex = max(maxLayerIndex, int(layerIndex))
								continue
							}
							rcSetCon(span, dir, int32(layerIndex))
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

func RcCreateHeightfield(sizeX, sizeZ int32,
	minBounds, maxBounds []float32,
	cellSize, cellHeight float32) (heightfield *RcHeightfield) {
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
