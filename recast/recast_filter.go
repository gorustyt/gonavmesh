package recast

import "github.com/gorustyt/gonavmesh/common"

const (
	/// The number of spans allocated per span spool.
	/// @see RcSpanPool
	RC_SPANS_PER_POOL = 2048
	/// Defines the number of bits allocated to RcSpan::smin and RcSpan::smax.
	RC_SPAN_HEIGHT_BITS = 13
	/// Defines the maximum value for RcSpan::smin and RcSpan::smax.
	RC_SPAN_MAX_HEIGHT = (1 << RC_SPAN_HEIGHT_BITS) - 1
	/// Represents the null area.
	/// When a data element is given this value it is considered to no longer be
	/// assigned to a usable area.  (E.g. It is un-walkable.)
	RC_NULL_AREA = 0
)
const MAX_HEIGHTFIELD_HEIGHT = 0xffff // TODO (graham): Move this to a more visible constant and update usages.

type RcSpan struct {
	Smin uint32  ///< The lower limit of the span. [Limit: < #smax]
	Smax uint32  ///< The upper limit of the span. [Limit: <= #RC_SPAN_MAX_HEIGHT]
	Area uint32  ///< The area id assigned to the span.
	Next *RcSpan ///< The next span higher up in column.
}

func NewRcSpan() *RcSpan {
	return &RcSpan{
		Smin: RC_SPAN_HEIGHT_BITS,
		Smax: RC_SPAN_HEIGHT_BITS,
		Area: 6,
	}
}

// / A memory pool used for quick allocation of spans within a heightfield.
// / @see RcHeightfield
type RcSpanPool struct {
	next  *RcSpanPool               ///< The next span pool.
	items [RC_SPANS_PER_POOL]RcSpan ///< Array of spans in the pool.
}

// / A dynamic heightfield representing obstructed space.
// / @ingroup recast
type RcHeightfield struct {
	Width    int32       ///< The width of the heightfield. (Along the x-axis in cell units.)
	Height   int32       ///< The height of the heightfield. (Along the z-axis in cell units.)
	Bmin     [3]float32  ///< The minimum bounds in world space. [(x, y, z)]
	Bmax     [3]float32  ///< The maximum bounds in world space. [(x, y, z)]
	Cs       float32     ///< The size of each cell. (On the xz-plane.)
	Ch       float32     ///< The height of each cell. (The minimum increment along the y-axis.)
	Spans    []*RcSpan   ///< Heightfield of spans (width*height).
	Pools    *RcSpanPool ///< Linked list of span pools.
	Freelist *RcSpan     ///< The next free span.
	// Explicitly-disabled copy constructor and copy assignment operator.
}

func RcFilterLowHangingWalkableObstacles(walkableClimb int32, heightfield *RcHeightfield) {
	xSize := heightfield.Width
	zSize := heightfield.Height

	for z := int32(0); z < zSize; z++ {
		for x := int32(0); x < xSize; x++ {
			var previousSpan *RcSpan
			previousWasWalkable := false
			previousArea := uint32(RC_NULL_AREA)

			for span := heightfield.Spans[x+z*xSize]; span != nil; {
				walkable := span.Area != RC_NULL_AREA
				// If current span is not walkable, but there is walkable
				// span just below it, mark the span above it walkable too.
				if !walkable && previousWasWalkable {
					if common.Abs(float64(span.Smax-previousSpan.Smax)) <= float64(walkableClimb) {
						span.Area = previousArea
					}
				}
				// Copy walkable flag so that it cannot propagate
				// past multiple non-walkable objects.
				previousWasWalkable = walkable
				previousArea = span.Area
				previousSpan = span
				span = span.Next
			}
		}
	}
}

func RcFilterLedgeSpans(walkableHeight int32, walkableClimb int32, heightfield *RcHeightfield) {
	xSize := heightfield.Width;
	zSize := heightfield.Height;

	// Mark spans that are adjacent to a ledge as unwalkable..
	for z := int32(0); z < zSize; z++ {
		for x := int32(0); x < xSize; x++ {
			for span := heightfield.Spans[x+z*xSize]; span != nil; span = span.Next {
				// Skip non-walkable spans.
				if (span.Area == RC_NULL_AREA) {
					continue;
				}

				floor := int32(span.Smax);
				ceiling := int32(MAX_HEIGHTFIELD_HEIGHT);
				if span.Next != nil {
					ceiling = int32(span.Next.Smin)
				}

				// The difference between this walkable area and the lowest neighbor walkable area.
				// This is the difference between the current span and all neighbor spans that have
				// enough space for an agent to move between, but not accounting at all for surface slope.
				lowestNeighborFloorDifference := int32(MAX_HEIGHTFIELD_HEIGHT)

				// Min and max height of accessible neighbours.
				lowestTraversableNeighborFloor := int32(span.Smax);
				highestTraversableNeighborFloor := int32(span.Smax);

				for direction := int32(0); direction < 4; direction++ {
					neighborX := x + common.GetDirOffsetX(direction);
					neighborZ := z + common.GetDirOffsetY(direction);

					// Skip neighbours which are out of bounds.
					if (neighborX < 0 || neighborZ < 0 || neighborX >= xSize || neighborZ >= zSize) {
						lowestNeighborFloorDifference = -walkableClimb - 1;
						break;
					}

					neighborSpan := heightfield.Spans[neighborX+neighborZ*xSize];

					// The most we can step down to the neighbor is the walkableClimb distance.
					// Start with the area under the neighbor span
					neighborCeiling := int32(MAX_HEIGHTFIELD_HEIGHT)
					if neighborSpan != nil {
						neighborCeiling = int32(neighborSpan.Smin)
					}
					// Skip neighbour if the gap between the spans is too small.
					if (min(ceiling, neighborCeiling)-floor >= walkableHeight) {
						lowestNeighborFloorDifference = (-walkableClimb - 1);
						break;
					}

					// For each span in the neighboring column...
					for ; neighborSpan != nil; neighborSpan = neighborSpan.Next {
						neighborFloor := int32(neighborSpan.Smax)
						neighborCeiling = MAX_HEIGHTFIELD_HEIGHT;
						if neighborSpan.Next != nil {
							neighborCeiling = int32(neighborSpan.Next.Smin)
						}
						// Only consider neighboring areas that have enough overlap to be potentially traversable.
						if (min(ceiling, neighborCeiling)-max(floor, neighborFloor) < walkableHeight) {
							// No space to traverse between them.
							continue;
						}

						neighborFloorDifference := neighborFloor - floor;
						lowestNeighborFloorDifference = min(lowestNeighborFloorDifference, neighborFloorDifference);

						// Find min/max accessible neighbor height.
						// Only consider neighbors that are at most walkableClimb away.
						if (common.Abs(neighborFloorDifference) <= walkableClimb) {
							// There is space to move to the neighbor cell and the slope isn't too much.
							lowestTraversableNeighborFloor = min(lowestTraversableNeighborFloor, neighborFloor);
							highestTraversableNeighborFloor = max(highestTraversableNeighborFloor, neighborFloor);
						} else if (neighborFloorDifference < -walkableClimb) {
							// We already know this will be considered a ledge span so we can early-out
							break;
						}
					}
				}

				// The current span is close to a ledge if the magnitude of the drop to any neighbour span is greater than the walkableClimb distance.
				// That is, there is a gap that is large enough to let an agent move between them, but the drop (surface slope) is too large to allow it.
				// (If this is the case, then biggestNeighborStepDown will be negative, so compare against the negative walkableClimb as a means of checking
				// the magnitude of the delta)
				if (lowestNeighborFloorDifference < -walkableClimb) {
					span.Area = RC_NULL_AREA;
					// If the difference between all neighbor floors is too large, this is a steep slope, so mark the span as an unwalkable ledge.
				} else if (highestTraversableNeighborFloor-lowestTraversableNeighborFloor > walkableClimb) {

					span.Area = RC_NULL_AREA;
				}
			}
		}
	}
}

func RcFilterWalkableLowHeightSpans(walkableHeight int32, heightfield *RcHeightfield) {
	xSize := heightfield.Width
	zSize := heightfield.Height

	// Remove walkable flag from spans which do not have enough
	// space above them for the agent to stand there.
	for z := int32(0); z < zSize; z++ {
		for x := int32(0); x < xSize; x++ {
			for span := heightfield.Spans[x+z*xSize]; span != nil; span = span.Next {
				floor := int32(span.Smax)
				ceiling := int32(MAX_HEIGHTFIELD_HEIGHT)
				if span.Next != nil {
					ceiling = int32(span.Next.Smin)
				}
				if ceiling-floor < walkableHeight {
					span.Area = RC_NULL_AREA
				}
			}
		}
	}
}
