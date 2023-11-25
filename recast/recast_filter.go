package recast

import "gonavamesh/common"

const (
	/// The number of spans allocated per span spool.
	/// @see rcSpanPool
	RC_SPANS_PER_POOL = 2048
	/// Defines the number of bits allocated to rcSpan::smin and rcSpan::smax.
	RC_SPAN_HEIGHT_BITS = 13
	/// Defines the maximum value for rcSpan::smin and rcSpan::smax.
	RC_SPAN_MAX_HEIGHT = (1 << RC_SPAN_HEIGHT_BITS) - 1
	/// Represents the null area.
	/// When a data element is given this value it is considered to no longer be
	/// assigned to a usable area.  (E.g. It is un-walkable.)
	RC_NULL_AREA = 0
)

type rcSpan struct {
	smin int     ///< The lower limit of the span. [Limit: < #smax]
	smax int     ///< The upper limit of the span. [Limit: <= #RC_SPAN_MAX_HEIGHT]
	area int     ///< The area id assigned to the span.
	next *rcSpan ///< The next span higher up in column.
}

func NewRcSpan() *rcSpan {
	return &rcSpan{
		smin: RC_SPAN_HEIGHT_BITS,
		smax: RC_SPAN_HEIGHT_BITS,
		area: 6,
	}
}

// / A memory pool used for quick allocation of spans within a heightfield.
// / @see RcHeightfield
type rcSpanPool struct {
	next  *rcSpanPool               ///< The next span pool.
	items [RC_SPANS_PER_POOL]rcSpan ///< Array of spans in the pool.
}

// / A dynamic heightfield representing obstructed space.
// / @ingroup recast
type RcHeightfield struct {
	width    int         ///< The width of the heightfield. (Along the x-axis in cell units.)
	height   int         ///< The height of the heightfield. (Along the z-axis in cell units.)
	bmin     [3]float64  ///< The minimum bounds in world space. [(x, y, z)]
	bmax     [3]float64  ///< The maximum bounds in world space. [(x, y, z)]
	cs       float64     ///< The size of each cell. (On the xz-plane.)
	ch       float64     ///< The height of each cell. (The minimum increment along the y-axis.)
	spans    []*rcSpan   ///< Heightfield of spans (width*height).
	pools    *rcSpanPool ///< Linked list of span pools.
	freelist *rcSpan     ///< The next free span.
	// Explicitly-disabled copy constructor and copy assignment operator.
}

func RcFilterLowHangingWalkableObstacles(walkableClimb int, heightfield *RcHeightfield) {
	xSize := heightfield.width
	zSize := heightfield.height

	for z := 0; z < zSize; z++ {
		for x := 0; x < xSize; x++ {
			var previousSpan *rcSpan
			previousWasWalkable := false
			previousArea := RC_NULL_AREA

			for span := heightfield.spans[x+z*xSize]; span != nil; {
				walkable := span.area != RC_NULL_AREA
				// If current span is not walkable, but there is walkable
				// span just below it, mark the span above it walkable too.
				if !walkable && previousWasWalkable {
					if common.Abs(span.smax-previousSpan.smax) <= walkableClimb {
						span.area = previousArea
					}
				}
				// Copy walkable flag so that it cannot propagate
				// past multiple non-walkable objects.
				previousWasWalkable = walkable
				previousArea = span.area
				previousSpan = span
				span = span.next
			}
		}
	}
}

func RcFilterLedgeSpans(walkableHeight int, walkableClimb int, heightfield *RcHeightfield) {
	xSize := heightfield.width
	zSize := heightfield.height
	MAX_HEIGHT := 0xffff // TODO (graham): Move this to a more visible constant and update usages.

	// Mark border spans.
	for z := 0; z < zSize; z++ {
		for x := 0; x < xSize; x++ {
			for span := heightfield.spans[x+z*xSize]; span != nil; span = span.next {
				// Skip non walkable spans.
				if span.area == RC_NULL_AREA {
					continue
				}

				bot := span.smax
				top := MAX_HEIGHT
				if span.next != nil {
					top = span.next.smin
				}
				// Find neighbours minimum height.
				minNeighborHeight := MAX_HEIGHT

				// Min and max height of accessible neighbours.
				accessibleNeighborMinHeight := span.smax
				accessibleNeighborMaxHeight := span.smax

				for direction := 0; direction < 4; direction++ {
					dx := x + common.GetDirOffsetX(direction)
					dy := z + common.GetDirOffsetY(direction)
					// Skip neighbours which are out of bounds.
					if dx < 0 || dy < 0 || dx >= xSize || dy >= zSize {
						minNeighborHeight = common.Min(minNeighborHeight, -walkableClimb-bot)
						continue
					}

					// From minus infinity to the first span.
					neighborSpan := heightfield.spans[dx+dy*xSize]
					neighborBot := -walkableClimb
					neighborTop := MAX_HEIGHT
					if neighborSpan != nil {
						neighborTop = neighborSpan.smin
					}
					// Skip neighbour if the gap between the spans is too small.
					if common.Min(top, neighborTop)-common.Max(bot, neighborBot) > walkableHeight {
						minNeighborHeight = common.Min(minNeighborHeight, neighborBot-bot)
					}

					// Rest of the spans.
					for neighborSpan := heightfield.spans[dx+dy*xSize]; neighborSpan != nil; neighborSpan = neighborSpan.next {
						neighborBot = neighborSpan.smax
						neighborTop := MAX_HEIGHT
						if neighborSpan.next != nil {
							neighborTop = neighborSpan.next.smin
						}

						// Skip neighbour if the gap between the spans is too small.
						if common.Min(top, neighborTop)-common.Max(bot, neighborBot) > walkableHeight {
							minNeighborHeight = common.Min(minNeighborHeight, neighborBot-bot)

							// Find min/max accessible neighbour height.
							if common.Abs(neighborBot-bot) <= walkableClimb {
								if neighborBot < accessibleNeighborMinHeight {
									accessibleNeighborMinHeight = neighborBot
								}
								if neighborBot > accessibleNeighborMaxHeight {
									accessibleNeighborMaxHeight = neighborBot
								}
							}

						}
					}
				}

				// The current span is close to a ledge if the drop to any
				// neighbour span is less than the walkableClimb.
				if minNeighborHeight < -walkableClimb {
					span.area = RC_NULL_AREA
					// If the difference between all neighbours is too large,
					// we are at steep slope, mark the span as ledge.
				} else if (accessibleNeighborMaxHeight - accessibleNeighborMinHeight) > walkableClimb {
					span.area = RC_NULL_AREA
				}
			}
		}
	}
}

func RcFilterWalkableLowHeightSpans(walkableHeight int, heightfield *RcHeightfield) {
	var (
		xSize      = heightfield.width
		zSize      = heightfield.height
		MAX_HEIGHT = 0xffff
	)
	// Remove walkable flag from spans which do not have enough
	// space above them for the agent to stand there.
	for z := 0; z < zSize; z++ {
		for x := 0; x < xSize; x++ {
			for span := heightfield.spans[x+z*xSize]; span != nil; span = span.next {
				bot := span.smax
				top := MAX_HEIGHT
				if span.next != nil {
					top = span.next.smin
				}
				if (top - bot) < walkableHeight {
					span.area = RC_NULL_AREA
				}
			}
		}
	}
}
