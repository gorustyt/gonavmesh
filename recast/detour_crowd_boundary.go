package recast

import (
	"gonavamesh/common"
	"math"
)

type Segment struct {
	s [6]float64 ///< Segment start/end
	d float64    ///< Distance for pruning.

}

type dtLocalBoundary struct {
	MAX_LOCAL_SEGS  int
	MAX_LOCAL_POLYS int
	m_center        [3]float64
	m_segs          []*Segment
	m_nsegs         int

	m_polys  []DtPolyRef
	m_npolys int
}

func (d *dtLocalBoundary) GetCenter() [3]float64       { return d.m_center }
func (d *dtLocalBoundary) GetSegmentCount() int        { return d.m_nsegs }
func (d *dtLocalBoundary) GetSegment(i int) [6]float64 { return d.m_segs[i].s }
func newDtLocalBoundary() *dtLocalBoundary {
	d := &dtLocalBoundary{
		MAX_LOCAL_SEGS:  8,
		MAX_LOCAL_POLYS: 16,
	}
	d.m_polys = make([]DtPolyRef, d.MAX_LOCAL_POLYS)
	d.m_segs = make([]*Segment, d.MAX_LOCAL_SEGS)
	return d
}

func (d *dtLocalBoundary) reset() {
	common.Vset(d.m_center[:], math.MaxFloat64, math.MaxFloat64, math.MaxFloat64)
	d.m_npolys = 0
	d.m_nsegs = 0
}
func (d *dtLocalBoundary) addSegment(dist float64, s []float64) {
	// Insert neighbour based on the distance.
	var seg *Segment
	if d.m_nsegs == 0 {
		// First, trivial accept.
		seg = d.m_segs[0]
	} else if dist >= d.m_segs[d.m_nsegs-1].d {
		// Further than the last segment, skip.
		if d.m_nsegs >= d.MAX_LOCAL_SEGS {
			return
		}

		// Last, trivial accept.
		seg = d.m_segs[d.m_nsegs]
	} else {
		// Insert inbetween.
		var i int
		for i = 0; i < d.m_nsegs; i++ {
			if dist <= d.m_segs[i].d {
				break
			}
		}

		tgt := i + 1
		n := common.Min(d.m_nsegs-i, d.MAX_LOCAL_SEGS-tgt)
		common.AssertTrue(tgt+n <= d.MAX_LOCAL_SEGS)
		if n > 0 {
			copy(d.m_segs[tgt:], d.m_segs[i:i+n])
		}

		seg = d.m_segs[i]
	}
	seg.d = dist
	copy(seg.s[:], s[:6])
	if d.m_nsegs < d.MAX_LOCAL_SEGS {
		d.m_nsegs++
	}

}

func (d *dtLocalBoundary) update(ref DtPolyRef, pos []float64, collisionQueryRange float64,
	navquery NavMeshQuery, filter *DtQueryFilter) {
	MAX_SEGS_PER_POLY := DT_VERTS_PER_POLYGON * 3

	if ref == 0 {
		common.Vset(d.m_center[:], math.MaxFloat64, math.MaxFloat64, math.MaxFloat64)
		d.m_nsegs = 0
		d.m_npolys = 0
		return
	}

	copy(d.m_center[:], pos)

	// First query non-overlapping polygons.
	d.m_npolys, _ = navquery.FindLocalNeighbourhood(ref, pos, collisionQueryRange, filter, d.m_polys, nil, d.MAX_LOCAL_POLYS)

	// Secondly, store all polygon edges.
	segs := make([]float64, MAX_SEGS_PER_POLY*6)
	nsegs := 0
	for j := 0; j < d.m_npolys; j++ {
		nsegs, _ = navquery.GetPolyWallSegments(d.m_polys[j], filter, segs, nil, MAX_SEGS_PER_POLY)
		for k := 0; k < nsegs; k++ {
			s := segs[k*6:]
			// Skip too distant segments.

			_, distSqr := dtDistancePtSegSqr2D(pos, s, s[3:])
			if distSqr > common.Sqr(collisionQueryRange) {
				continue
			}

			d.addSegment(distSqr, s)
		}
	}
}

func (d *dtLocalBoundary) isValid(navquery NavMeshQuery, filter *DtQueryFilter) bool {
	if d.m_npolys == 0 {
		return false
	}
	// Check that all polygons still pass query filter.
	for i := 0; i < d.m_npolys; i++ {
		if !navquery.IsValidPolyRef(d.m_polys[i], filter) {
			return false
		}

	}

	return true
}
