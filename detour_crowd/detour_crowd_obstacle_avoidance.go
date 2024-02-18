package detour_crowd

import (
	"github.com/gorustyt/gonavmesh/common"
	"github.com/gorustyt/gonavmesh/detour"
	"math"
)

type dtObstacleCircle struct {
	p      [3]float32 ///< Position of the obstacle
	vel    [3]float32 ///< Velocity of the obstacle
	dvel   [3]float32 ///< Velocity of the obstacle
	rad    float32    ///< Radius of the obstacle
	dp, np [3]float32 ///< Use for side selection during sampling.
}

type dtObstacleSegment struct {
	p, q  [3]float32 ///< End points of the obstacle segment
	touch bool
}

type DtObstacleAvoidanceDebugData struct {
	m_nsamples   int
	m_maxSamples int
	m_vel        []float32
	m_ssize      []float32
	m_pen        []float32
	m_vpen       []float32
	m_vcpen      []float32
	m_spen       []float32
	m_tpen       []float32
}

func (d *DtObstacleAvoidanceDebugData) GetSampleCount() int               { return d.m_nsamples }
func (d *DtObstacleAvoidanceDebugData) GetSampleVelocity(i int) []float32 { return d.m_vel[i*3:] }
func (d *DtObstacleAvoidanceDebugData) GetSampleSize(i int) float32       { return d.m_ssize[i] }
func (d *DtObstacleAvoidanceDebugData) GetSamplePenalty(i int) float32    { return d.m_pen[i] }
func (d *DtObstacleAvoidanceDebugData) GetSampleDesiredVelocityPenalty(i int) float32 {
	return d.m_vpen[i]
}
func (d *DtObstacleAvoidanceDebugData) GetSampleCurrentVelocityPenalty(i int) float32 {
	return d.m_vcpen[i]
}
func (d *DtObstacleAvoidanceDebugData) GetSamplePreferredSidePenalty(i int) float32 {
	return d.m_spen[i]
}
func (d *DtObstacleAvoidanceDebugData) GetSampleCollisionTimePenalty(i int) float32 {
	return d.m_tpen[i]
}

const DT_MAX_PATTERN_DIVS = 32 ///< Max numver of adaptive divs.
const DT_MAX_PATTERN_RINGS = 4 ///< Max number of adaptive rings.

type DtObstacleAvoidanceParams struct {
	VelBias       float32
	WeightDesVel  float32
	WeightCurVel  float32
	WeightSide    float32
	WeightToi     float32
	HorizTime     float32
	GridSize      int ///< grid
	AdaptiveDivs  int ///< adaptive
	AdaptiveRings int ///< adaptive
	AdaptiveDepth int ///< adaptive
}

type dtObstacleAvoidanceQuery struct {
	m_params       DtObstacleAvoidanceParams
	m_invHorizTime float32
	m_vmax         float32
	m_invVmax      float32

	m_maxCircles int
	m_circles    []*dtObstacleCircle
	m_ncircles   int

	m_maxSegments int
	m_segments    []*dtObstacleSegment
	m_nsegments   int
}

const DT_PI = 3.14159265

func sweepCircleCircle(c0 []float32, r0 float32, v, c1 []float32, r1 float32, tmin, tmax *float32) int {
	EPS := float32(0.0001)
	s := make([]float32, 3)
	common.Vsub(s, c1, c0)
	r := r0 + r1
	c := common.Vdot2D(s, s) - r*r
	a := common.Vdot2D(v, v)
	if a < EPS {
		return 0
	} // not moving

	// Overlap, calc time to exit.
	b := common.Vdot2D(v, s)
	d := b*b - a*c
	if d < 0.0 {
		return 0
	} // no intersection.
	a = 1.0 / a
	rd := math.Sqrt(float64(d))
	*tmin = float32(float64(b)-rd) * a
	*tmax = float32(float64(b)+rd) * a
	return 1
}

func isectRaySeg(ap, u []float32, bp, bq []float32, t *float32) int {
	v := make([]float32, 3)
	w := make([]float32, 3)
	common.Vsub(v, bq, bp)
	common.Vsub(w, ap, bp)
	d := common.Vperp2D(u, v)
	if math.Abs(float64(d)) < 1e-6 {
		return 0
	}
	d = 1.0 / d
	*t = common.Vperp2D(v, w) * d
	if *t < 0 || *t > 1 {
		return 0
	}
	s := common.Vperp2D(u, w) * d
	if s < 0 || s > 1 {
		return 0
	}
	return 1
}
func NewDtObstacleAvoidanceDebugData(maxSamples int) *DtObstacleAvoidanceDebugData {
	d := &DtObstacleAvoidanceDebugData{}
	d.m_maxSamples = maxSamples
	d.m_vel = make([]float32, 3*d.m_maxSamples)
	d.m_pen = make([]float32, d.m_maxSamples)
	d.m_ssize = make([]float32, d.m_maxSamples)
	d.m_vpen = make([]float32, d.m_maxSamples)
	d.m_vcpen = make([]float32, d.m_maxSamples)
	d.m_spen = make([]float32, d.m_maxSamples)
	d.m_tpen = make([]float32, d.m_maxSamples)
	return d
}
func (d *DtObstacleAvoidanceDebugData) reset() {
	d.m_nsamples = 0
}

func (d *DtObstacleAvoidanceDebugData) addSample(vel []float32, ssize float32, pen float32,
	vpen, vcpen, spen, tpen float32) {
	if d.m_nsamples >= d.m_maxSamples {
		return
	}
	copy(d.m_vel[d.m_nsamples*3:], vel)
	d.m_ssize[d.m_nsamples] = ssize
	d.m_pen[d.m_nsamples] = pen
	d.m_vpen[d.m_nsamples] = vpen
	d.m_vcpen[d.m_nsamples] = vcpen
	d.m_spen[d.m_nsamples] = spen
	d.m_tpen[d.m_nsamples] = tpen
	d.m_nsamples++
}

func (d *DtObstacleAvoidanceDebugData) NormalizeArray(arr []float32, n int) {
	// Normalize penaly range.
	minPen := float32(math.MaxFloat32)
	maxPen := float32(math.SmallestNonzeroFloat32)
	for i := 0; i < n; i++ {
		minPen = min(minPen, arr[i])
		maxPen = min(maxPen, arr[i])
	}
	penRange := maxPen - minPen
	s := float32(1)
	if penRange > 0.001 {
		s = 1.0 / penRange
	}
	for i := 0; i < n; i++ {
		arr[i] = common.Clamp((arr[i]-minPen)*s, 0.0, 1.0)
	}

}

func (d *DtObstacleAvoidanceDebugData) NormalizeSamples() {
	d.NormalizeArray(d.m_pen, d.m_nsamples)
	d.NormalizeArray(d.m_vpen, d.m_nsamples)
	d.NormalizeArray(d.m_vcpen, d.m_nsamples)
	d.NormalizeArray(d.m_spen, d.m_nsamples)
	d.NormalizeArray(d.m_tpen, d.m_nsamples)
}

func (d *dtObstacleAvoidanceQuery) init(maxCircles, maxSegments int) bool {
	d.m_maxCircles = maxCircles
	d.m_ncircles = 0
	d.m_circles = make([]*dtObstacleCircle, d.m_maxCircles)
	d.m_maxSegments = maxSegments
	d.m_nsegments = 0
	d.m_segments = make([]*dtObstacleSegment, d.m_maxSegments)
	return true
}

func (d *dtObstacleAvoidanceQuery) reset() {
	d.m_ncircles = 0
	d.m_nsegments = 0
}

func (d *dtObstacleAvoidanceQuery) addCircle(pos []float32, rad float32,
	vel, dvel []float32) {
	if d.m_ncircles >= d.m_maxCircles {
		return
	}

	cir := d.m_circles[d.m_ncircles]
	d.m_ncircles++
	copy(cir.p[:], pos)
	cir.rad = rad
	copy(cir.vel[:], vel)
	copy(cir.dvel[:], dvel)
}

func (d *dtObstacleAvoidanceQuery) addSegment(p, q []float32) {
	if d.m_nsegments >= d.m_maxSegments {
		return
	}

	seg := d.m_segments[d.m_nsegments]
	d.m_nsegments++
	copy(seg.p[:], p)
	copy(seg.q[:], q)
}

func (d *dtObstacleAvoidanceQuery) prepare(pos, dvel []float32) {
	// Prepare obstacles
	for i := 0; i < d.m_ncircles; i++ {
		cir := d.m_circles[i]

		// Side
		pa := pos
		pb := cir.p
		dv := make([]float32, 3)
		orig := []float32{0, 0, 0}
		common.Vsub(cir.dp[:], pb[:], pa)
		common.Vnormalize(cir.dp[:])
		common.Vsub(dv, cir.dvel[:], dvel)

		a := common.TriArea2D(orig, cir.dp[:], dv)
		if a < 0.01 {
			cir.np[0] = -cir.dp[2]
			cir.np[2] = cir.dp[0]
		} else {
			cir.np[0] = cir.dp[2]
			cir.np[2] = -cir.dp[0]
		}
	}

	for i := 0; i < d.m_nsegments; i++ {
		seg := d.m_segments[i]

		// Precalc if the agent is really close to the segment.
		r := 0.01
		_, res := detour.DtDistancePtSegSqr2D(pos, seg.p[:], seg.q[:])
		seg.touch = float64(res) < common.Sqr(r)
	}
}

/* Calculate the collision penalty for a given velocity vector
 *
 * @param vcand sampled velocity
 * @param dvel desired velocity
 * @param minPenalty threshold penalty for early out
 */
func (d *dtObstacleAvoidanceQuery) processSample(vcand []float32, cs float32,
	pos []float32, rad float32,
	vel, dvel []float32,
	minPenalty float32,
	debug *DtObstacleAvoidanceDebugData) float32 {
	// penalty for straying away from the desired and current velocities
	vpen := d.m_params.WeightDesVel * (common.Vdist2D(vcand, dvel) * d.m_invVmax)
	vcpen := d.m_params.WeightCurVel * (common.Vdist2D(vcand, vel) * d.m_invVmax)

	// find the threshold hit time to bail out based on the early out penalty
	// (see how the penalty is calculated below to understand)
	minPen := minPenalty - vpen - vcpen
	tThresold := (d.m_params.WeightToi/minPen - 0.1) * d.m_params.HorizTime
	if tThresold-d.m_params.HorizTime > math.MaxFloat32 {
		return minPenalty
	} // already too much
	// Find min time of impact and exit amongst all obstacles.
	tmin := d.m_params.HorizTime
	side := float32(0)
	nside := 0

	for i := 0; i < d.m_ncircles; i++ {
		cir := d.m_circles[i]

		// RVO
		vab := make([]float32, 3)
		common.Vscale(vab, vcand, 2)
		common.Vsub(vab, vab, vel)
		common.Vsub(vab, vab, cir.vel[:])

		// Side
		side += common.Clamp(common.Min(common.Vdot2D(cir.dp[:], vab)*0.5+0.5, common.Vdot2D(cir.np[:], vab)*2), 0.0, 1.0)
		nside++

		var htmin = float32(0)
		var htmax = float32(0)
		if sweepCircleCircle(pos, rad, vab, cir.p[:], cir.rad, &htmin, &htmax) == 0 {
			continue
		}

		// Handle overlapping obstacles.
		if htmin < 0.0 && htmax > 0.0 {
			// Avoid more when overlapped.
			htmin = -htmin * 0.5
		}

		if htmin >= 0.0 {
			// The closest obstacle is somewhere ahead of us, keep track of nearest obstacle.
			if htmin < tmin {
				tmin = htmin
				if tmin < tThresold {
					return minPenalty
				}

			}
		}
	}

	for i := 0; i < d.m_nsegments; i++ {
		seg := d.m_segments[i]
		htmin := float32(0)

		if seg.touch {
			// Special case when the agent is very close to the segment.
			sdir := make([]float32, 3)
			snorm := make([]float32, 3)
			common.Vsub(sdir, seg.q[:], seg.p[:])
			snorm[0] = -sdir[2]
			snorm[2] = sdir[0]
			// If the velocity is pointing towards the segment, no collision.
			if common.Vdot2D(snorm, vcand) < 0.0 {
				continue
			}
			// Else immediate collision.
			htmin = 0.0
		} else {
			if isectRaySeg(pos, vcand, seg.p[:], seg.q[:], &htmin) == 0 {
				continue
			}

		}

		// Avoid less when facing walls.
		htmin *= 2.0

		// The closest obstacle is somewhere ahead of us, keep track of nearest obstacle.
		if htmin < tmin {
			tmin = htmin
			if tmin < tThresold {
				return minPenalty
			}

		}
	}

	// Normalize side bias, to prevent it dominating too much.
	if nside != 0 {
		side /= float32(nside)
	}

	spen := d.m_params.WeightSide * side
	tpen := d.m_params.WeightToi * (1.0 / (0.1 + tmin*d.m_invHorizTime))

	penalty := vpen + vcpen + spen + tpen

	// Store different penalties for debug viewing
	if debug != nil {
		debug.addSample(vcand, cs, penalty, vpen, vcpen, spen, tpen)
	}

	return penalty
}

func (d *dtObstacleAvoidanceQuery) sampleVelocityGrid(pos []float32, rad float32, vmax float32, vel, dvel, nvel []float32,
	params *DtObstacleAvoidanceParams,
	debug *DtObstacleAvoidanceDebugData) int32 {
	d.prepare(pos, dvel)
	d.m_params = *params
	d.m_invHorizTime = 1.0 / d.m_params.HorizTime
	d.m_vmax = vmax
	d.m_invVmax = math.MaxFloat32
	if vmax > 0 {
		d.m_invVmax = 1.0 / vmax
	}

	common.Vset(nvel, 0, 0, 0)

	if debug != nil {
		debug.reset()
	}

	cvx := dvel[0] * d.m_params.VelBias
	cvz := dvel[2] * d.m_params.VelBias
	cs := vmax * 2 * (1 - d.m_params.VelBias) / float32((d.m_params.GridSize - 1))
	half := float32(d.m_params.GridSize-1) * cs * 0.5

	minPenalty := float32(math.MaxFloat32)
	ns := int32(0)

	for y := 0; y < d.m_params.GridSize; y++ {
		for x := 0; x < d.m_params.GridSize; x++ {
			vcand := make([]float32, 3)
			vcand[0] = cvx + float32(x)*cs - half
			vcand[1] = 0
			vcand[2] = cvz + float32(y)*cs - half

			if common.Sqr(vcand[0])+common.Sqr(vcand[2]) > common.Sqr(vmax+cs/2) {
				continue
			}

			penalty := d.processSample(vcand, cs, pos, rad, vel, dvel, minPenalty, debug)
			ns++
			if penalty < minPenalty {
				minPenalty = penalty
				copy(nvel, vcand)
			}
		}
	}

	return ns
}

// vector normalization that ignores the y-component.
func dtNormalize2D(v []float32) {
	d := float32(math.Sqrt(float64(v[0]*v[0] + v[2]*v[2])))
	if d == 0 {
		return
	}

	d = 1.0 / d
	v[0] *= d
	v[2] *= d
}

// vector normalization that ignores the y-component.
func dtRorate2D(dest, v []float32, ang float32) {
	c := float32(math.Cos(float64(ang)))
	s := float32(math.Sin(float64(ang)))
	dest[0] = v[0]*c - v[2]*s
	dest[2] = v[0]*s + v[2]*c
	dest[1] = v[1]
}

func (d *dtObstacleAvoidanceQuery) sampleVelocityAdaptive(pos []float32, rad, vmax float32,
	vel, dvel, nvel []float32, params *DtObstacleAvoidanceParams,
	debug *DtObstacleAvoidanceDebugData) int32 {
	d.prepare(pos, dvel)
	d.m_params = *params
	d.m_invHorizTime = 1.0 / d.m_params.HorizTime
	d.m_vmax = vmax
	d.m_invVmax = math.MaxFloat32
	if vmax > 0 {
		d.m_invVmax = 1.0 / vmax
	}

	common.Vset(nvel, 0, 0, 0)

	if debug != nil {
		debug.reset()
	}

	// Build sampling pattern aligned to desired velocity.
	pat := make([]float32, (DT_MAX_PATTERN_DIVS*DT_MAX_PATTERN_RINGS+1)*2)
	npat := 0

	ndivs := d.m_params.AdaptiveDivs
	nrings := d.m_params.AdaptiveRings
	depth := d.m_params.AdaptiveDepth

	nd := common.Clamp(ndivs, 1, DT_MAX_PATTERN_DIVS)
	nr := common.Clamp(nrings, 1, DT_MAX_PATTERN_RINGS)
	da := (1.0 / float32(nd)) * DT_PI * 2
	ca := float32(math.Cos(float64(da)))
	sa := float32(math.Sin(float64(da)))

	// desired direction
	ddir := make([]float32, 6)
	copy(ddir, dvel)
	dtNormalize2D(ddir)
	dtRorate2D(ddir[3:], ddir, da*0.5) // rotated by da/2

	// Always add sample at zero
	pat[npat*2+0] = 0
	pat[npat*2+1] = 0
	npat++

	for j := 0; j < nr; j++ {
		r := float32(nr-j) / float32(nr)
		pat[npat*2+0] = ddir[(j%2)*3] * r
		pat[npat*2+1] = ddir[(j%2)*3+2] * r
		last1 := pat[npat*2:]
		last2 := last1
		npat++

		for i := 1; i < nd-1; i += 2 {
			// get next point on the "right" (rotate CW)
			pat[npat*2+0] = last1[0]*ca + last1[1]*sa
			pat[npat*2+1] = -last1[0]*sa + last1[1]*ca
			// get next point on the "left" (rotate CCW)
			pat[npat*2+2] = last2[0]*ca - last2[1]*sa
			pat[npat*2+3] = last2[0]*sa + last2[1]*ca

			last1 = pat[:npat*2]
			last2 = last1[2:]
			npat += 2
		}

		if (nd & 1) == 0 {
			pat[npat*2+2] = last2[0]*ca - last2[1]*sa
			pat[npat*2+3] = last2[0]*sa + last2[1]*ca
			npat++
		}
	}

	// Start sampling.
	cr := vmax * (1.0 - d.m_params.VelBias)
	res := make([]float32, 3)
	common.Vset(res, dvel[0]*d.m_params.VelBias, 0, dvel[2]*d.m_params.VelBias)
	ns := 0

	for k := 0; k < depth; k++ {
		minPenalty := float32(math.MaxFloat32)
		bvel := make([]float32, 3)
		common.Vset(bvel, 0, 0, 0)

		for i := 0; i < npat; i++ {
			vcand := make([]float32, 3)
			vcand[0] = res[0] + pat[i*2+0]*cr
			vcand[1] = 0
			vcand[2] = res[2] + pat[i*2+1]*cr

			if common.Sqr(vcand[0])+common.Sqr(vcand[2]) > common.Sqr(vmax+0.001) {
				continue
			}

			penalty := d.processSample(vcand, cr/10, pos, rad, vel, dvel, minPenalty, debug)
			ns++
			if penalty < minPenalty {
				minPenalty = penalty
				copy(bvel, vcand)
			}
		}

		copy(res, bvel)

		cr *= 0.5
	}

	copy(nvel, res)

	return int32(ns)
}
