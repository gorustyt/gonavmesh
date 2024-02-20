package mesh

import (
	"github.com/gorustyt/gonavmesh/common"
	"github.com/gorustyt/gonavmesh/debug_utils"
	"github.com/gorustyt/gonavmesh/recast"
	"math"
)

// Quick and dirty convex hull.

// Returns true if 'c' is left of line 'a'-'b'.
func left(a, b, c []float64) bool {
	u1 := b[0] - a[0]
	v1 := b[2] - a[2]
	u2 := c[0] - a[0]
	v2 := c[2] - a[2]
	return u1*v2-v1*u2 < 0
}

// Returns true if 'a' is more lower-left than 'b'.
func cmppt(a, b []float64) bool {
	if a[0] < b[0] {
		return true
	}
	if a[0] > b[0] {
		return false
	}
	if a[2] < b[2] {
		return true
	}
	if a[2] > b[2] {
		return false
	}
	return false
}

// Calculates convex hull on xz-plane of points on 'pts',
// stores the indices of the resulting hull in 'out' and
// returns number of points on hull.
func convexhull(pts []float64, npts int, out []int) int {
	// Find lower-leftmost point.
	hull := 0
	for i := 1; i < npts; i++ {
		if cmppt(pts[i*3:], pts[hull*3:]) {
			hull = i
		}
	}

	// Gift wrap hull.
	endpt := 0
	i := 0
	common.DoWhile(func() (stop bool) {

		out[i] = hull
		i++
		endpt = 0
		for j := 1; j < npts; j++ {
			if hull == endpt || left(pts[hull*3:], pts[endpt*3:], pts[j*3:]) {
				endpt = j
			}
		}
		hull = endpt
		return false
	}, func() bool {
		return endpt != out[0]
	})

	return i
}

func pointInPoly(nvert int, verts []float64, p []float64) int {
	var i, j, c int
	i = 0
	j = nvert - 1
	for i < nvert {
		vi := verts[i*3:]
		vj := verts[j*3:]
		if ((vi[2] > p[2]) != (vj[2] > p[2])) && (p[0] < (vj[0]-vi[0])*(p[2]-vi[2])/(vj[2]-vi[2])+vi[0]) {
			c = 1
		}

		j = i
		i++
	}
	return c
}

const (
	MAX_PTS = 12
)

type ConvexVolumeTool struct {
	m_sample     *Sample
	m_areaType   int
	m_polyOffset float64
	m_boxHeight  float64
	m_boxDescent float64

	m_pts   []float64
	m_npts  int
	m_hull  []int
	m_nhull int

	ctx *Content
}

func newConvexVolumeTool(ctx *Content) *ConvexVolumeTool {
	c := &ConvexVolumeTool{
		ctx:    ctx,
		m_pts:  make([]float64, MAX_PTS*3),
		m_hull: make([]int, MAX_PTS),
	}
	ctx.GetConfig().ToolsConfig.OnClearShapeClick = func() {
		c.reset()
	}
	return c
}

func (c *ConvexVolumeTool) Type() int           { return TOOL_CONVEX_VOLUME }
func (c *ConvexVolumeTool) init(sample *Sample) { c.m_sample = sample }
func (c *ConvexVolumeTool) reset() {
	c.m_npts = 0
	c.m_nhull = 0
}

func (c *ConvexVolumeTool) handleClick(s []float64, p []float64, shift bool) {
	if c.m_sample == nil {
		return
	}
	geom := c.m_sample.getInputGeom()
	if geom == nil {
		return
	}

	if shift {
		// Delete
		nearestIndex := -1
		vols := geom.getConvexVolumes()
		for i := 0; i < geom.getConvexVolumeCount(); i++ {
			if pointInPoly(vols[i].nverts, vols[i].verts[:], p) > 0 &&
				p[1] >= vols[i].hmin && p[1] <= vols[i].hmax {
				nearestIndex = i
			}
		}
		// If end point close enough, delete it.
		if nearestIndex != -1 {
			geom.deleteConvexVolume(nearestIndex)
		}
	} else {
		// Create

		// If clicked on that last pt, create the shape.
		if c.m_npts != 0 && common.VdistSqr(p, c.m_pts[(c.m_npts-1)*3:]) < common.Sqr(0.2) {
			if c.m_nhull > 2 {
				// Create shape.
				verts := make([]float64, MAX_PTS*3)
				for i := 0; i < c.m_nhull; i++ {
					copy(common.GetVert3(verts, i), common.GetVert3(c.m_pts, c.m_hull[i]))
				}

				minh := math.MaxFloat64
				maxh := 0.0
				for i := 0; i < c.m_nhull; i++ {
					minh = min(minh, verts[i*3+1])
				}

				minh -= c.m_boxDescent
				maxh = minh + c.m_boxHeight

				if c.m_polyOffset > 0.01 {
					offset := make([]float64, MAX_PTS*2*3)
					noffset := recast.RcOffsetPoly(verts, c.m_nhull, c.m_polyOffset, offset, MAX_PTS*2)
					if noffset > 0 {
						geom.addConvexVolume(offset, noffset, minh, maxh, c.m_areaType)
					}
				} else {
					geom.addConvexVolume(verts, c.m_nhull, minh, maxh, c.m_areaType)
				}
			}

			c.m_npts = 0
			c.m_nhull = 0
		} else {
			// Add new point
			if c.m_npts < MAX_PTS {
				copy(common.GetVert3(c.m_pts, c.m_npts), p)
				c.m_npts++
				// Update hull.
				if c.m_npts > 1 {
					c.m_nhull = convexhull(c.m_pts, c.m_npts, c.m_hull)
				} else {
					c.m_nhull = 0
				}

			}
		}
	}
}
func (c *ConvexVolumeTool) handleRender() {
	dd := c.m_sample.getDebugDraw()

	// Find height extent of the shape.
	minh := math.MaxFloat64
	maxh := 0.0
	for i := 0; i < c.m_npts; i++ {
		minh = min(minh, c.m_pts[i*3+1])
	}

	minh -= c.m_boxDescent
	maxh = minh + c.m_boxHeight

	dd.Begin(debug_utils.DU_DRAW_POINTS, 4.0)
	for i := 0; i < c.m_npts; i++ {
		col := debug_utils.DuRGBA(255, 255, 255, 255)
		if i == c.m_npts-1 {
			col = debug_utils.DuRGBA(240, 32, 16, 255)
		}

		dd.Vertex1(c.m_pts[i*3+0], c.m_pts[i*3+1]+0.1, c.m_pts[i*3+2], col)
	}
	dd.End()

	dd.Begin(debug_utils.DU_DRAW_LINES, 2.0)
	i := 0
	j := c.m_nhull - 1
	for i < c.m_nhull {
		vi := c.m_pts[c.m_hull[j]*3:]
		vj := c.m_pts[c.m_hull[i]*3:]
		dd.Vertex1(vj[0], minh, vj[2], debug_utils.DuRGBA(255, 255, 255, 64))
		dd.Vertex1(vi[0], minh, vi[2], debug_utils.DuRGBA(255, 255, 255, 64))
		dd.Vertex1(vj[0], maxh, vj[2], debug_utils.DuRGBA(255, 255, 255, 64))
		dd.Vertex1(vi[0], maxh, vi[2], debug_utils.DuRGBA(255, 255, 255, 64))
		dd.Vertex1(vj[0], minh, vj[2], debug_utils.DuRGBA(255, 255, 255, 64))
		dd.Vertex1(vj[0], maxh, vj[2], debug_utils.DuRGBA(255, 255, 255, 64))
		j = i
		i++
	}
	dd.End()
}
