package recast

import "gonavamesh/common"

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

type RcSpan struct {
	Smin int     ///< The lower limit of the span. [Limit: < #smax]
	Smax int     ///< The upper limit of the span. [Limit: <= #RC_SPAN_MAX_HEIGHT]
	Area int     ///< The area id assigned to the span.
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
	Width    int         ///< The width of the heightfield. (Along the x-axis in cell units.)
	Height   int         ///< The height of the heightfield. (Along the z-axis in cell units.)
	Bmin     [3]float64  ///< The minimum bounds in world space. [(x, y, z)]
	Bmax     [3]float64  ///< The maximum bounds in world space. [(x, y, z)]
	Cs       float64     ///< The size of each cell. (On the xz-plane.)
	Ch       float64     ///< The height of each cell. (The minimum increment along the y-axis.)
	Spans    []*RcSpan   ///< Heightfield of spans (width*height).
	Pools    *RcSpanPool ///< Linked list of span pools.
	Freelist *RcSpan     ///< The next free span.
	// Explicitly-disabled copy constructor and copy assignment operator.
}

func RcFilterLowHangingWalkableObstacles(walkableClimb int, heightfield *RcHeightfield) {
	xSize := heightfield.Width
	zSize := heightfield.Height

	for z := 0; z < zSize; z++ {
		for x := 0; x < xSize; x++ {
			var previousSpan *RcSpan
			previousWasWalkable := false
			previousArea := RC_NULL_AREA

			for span := heightfield.Spans[x+z*xSize]; span != nil; {
				walkable := span.Area != RC_NULL_AREA
				// If current span is not walkable, but there is walkable
				// span just below it, mark the span above it walkable too.
				if !walkable && previousWasWalkable {
					if common.Abs(span.Smax-previousSpan.Smax) <= walkableClimb {
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

func RcFilterLedgeSpans(walkableHeight int, walkableClimb int, heightfield *RcHeightfield) {
	xSize := heightfield.Width
	zSize := heightfield.Height
	MAX_HEIGHT := 0xffff // TODO (graham): Move this to a more visible constant and update usages.

	// Mark border spans.
	for z := 0; z < zSize; z++ {
		for x := 0; x < xSize; x++ {
			for span := heightfield.Spans[x+z*xSize]; span != nil; span = span.Next {
				// Skip non walkable spans.
				if span.Area == RC_NULL_AREA {
					continue
				}

				bot := span.Smax
				top := MAX_HEIGHT
				if span.Next != nil {
					top = span.Next.Smin
				}
				// Find neighbours minimum height.
				minNeighborHeight := MAX_HEIGHT

				// Min and max height of accessible neighbours.
				accessibleNeighborMinHeight := span.Smax
				accessibleNeighborMaxHeight := span.Smax

				for direction := 0; direction < 4; direction++ {
					dx := x + common.GetDirOffsetX(direction)
					dy := z + common.GetDirOffsetY(direction)
					// Skip neighbours which are out of bounds.
					if dx < 0 || dy < 0 || dx >= xSize || dy >= zSize {
						minNeighborHeight = common.Min(minNeighborHeight, -walkableClimb-bot)
						continue
					}

					// From minus infinity to the first span.
					neighborSpan := heightfield.Spans[dx+dy*xSize]
					neighborBot := -walkableClimb
					neighborTop := MAX_HEIGHT
					if neighborSpan != nil {
						neighborTop = neighborSpan.Smin
					}
					// Skip neighbour if the gap between the spans is too small.
					if common.Min(top, neighborTop)-common.Max(bot, neighborBot) > walkableHeight {
						minNeighborHeight = common.Min(minNeighborHeight, neighborBot-bot)
					}

					// Rest of the spans.
					for neighborSpan := heightfield.Spans[dx+dy*xSize]; neighborSpan != nil; neighborSpan = neighborSpan.Next {
						neighborBot = neighborSpan.Smax
						neighborTop := MAX_HEIGHT
						if neighborSpan.Next != nil {
							neighborTop = neighborSpan.Next.Smin
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
					span.Area = RC_NULL_AREA
					// If the difference between all neighbours is too large,
					// we are at steep slope, mark the span as ledge.
				} else if (accessibleNeighborMaxHeight - accessibleNeighborMinHeight) > walkableClimb {
					span.Area = RC_NULL_AREA
				}
			}
		}
	}
}

func RcFilterWalkableLowHeightSpans(walkableHeight int, heightfield *RcHeightfield) {
	var (
		xSize      = heightfield.Width
		zSize      = heightfield.Height
		MAX_HEIGHT = 0xffff
	)
	// Remove walkable flag from spans which do not have enough
	// space above them for the agent to stand there.
	for z := 0; z < zSize; z++ {
		for x := 0; x < xSize; x++ {
			for span := heightfield.Spans[x+z*xSize]; span != nil; span = span.Next {
				bot := span.Smax
				top := MAX_HEIGHT
				if span.Next != nil {
					top = span.Next.Smin
				}
				if (top - bot) < walkableHeight {
					span.Area = RC_NULL_AREA
				}
			}
		}
	}
}
