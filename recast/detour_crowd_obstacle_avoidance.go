package recast

import "math"

type dtObstacleCircle struct {
	p      [3]float64 ///< Position of the obstacle
	vel    [3]float64 ///< Velocity of the obstacle
	dvel   [3]float64 ///< Velocity of the obstacle
	rad    float64    ///< Radius of the obstacle
	dp, np [3]float64 ///< Use for side selection during sampling.
}

type dtObstacleSegment struct {
	p, q  [3]float64 ///< End points of the obstacle segment
	touch bool
}

type dtObstacleAvoidanceDebugData struct {
	m_nsamples   int
	m_maxSamples int
	m_vel        []float64
	m_ssize      []float64
	m_pen        []float64
	m_vpen       []float64
	m_vcpen      []float64
	m_spen       []float64
	m_tpen       []float64
}

const DT_MAX_PATTERN_DIVS = 32 ///< Max numver of adaptive divs.
const DT_MAX_PATTERN_RINGS = 4 ///< Max number of adaptive rings.

type dtObstacleAvoidanceParams struct {
	velBias       float64
	weightDesVel  float64
	weightCurVel  float64
	weightSide    float64
	weightToi     float64
	horizTime     float64
	gridSize      int ///< grid
	adaptiveDivs  int ///< adaptive
	adaptiveRings int ///< adaptive
	adaptiveDepth int ///< adaptive
}

type dtObstacleAvoidanceQuery struct {
	m_params       dtObstacleAvoidanceParams
	m_invHorizTime float64
	m_vmax         float64
	m_invVmax      float64

	m_maxCircles int
	m_circles    []*dtObstacleCircle
	m_ncircles   int

	m_maxSegments int
	m_segments    []*dtObstacleSegment
	m_nsegments   int
}

const DT_PI = 3.14159265

func sweepCircleCircle(c0 []float64, r0 float64, v, c1 []float64, r1 float64, tmin, tmax *float64) int {
	EPS := 0.0001
	s := dtVsub(c1, c0)
	r := r0 + r1
	c := dtVdot2D(s, s) - r*r
	a := dtVdot2D(v, v)
	if a < EPS {
		return 0
	} // not moving

	// Overlap, calc time to exit.
	b := dtVdot2D(v, s)
	d := b*b - a*c
	if d < 0.0 {
		return 0
	} // no intersection.
	a = 1.0 / a
	rd := math.Sqrt(d)
	*tmin = (b - rd) * a
	*tmax = (b + rd) * a
	return 1
}

func isectRaySeg(ap, u []float64, bp, bq []float64, t *float64) int {
	v := dtVsub(bq, bp)
	w := dtVsub(ap, bp)
	d := dtVperp2D(u, v)
	if math.Abs(d) < 1e-6 {
		return 0
	}
	d = 1.0 / d
	*t = dtVperp2D(v, w) * d
	if *t < 0 || *t > 1 {
		return 0
	}
	s := dtVperp2D(u, w) * d
	if s < 0 || s > 1 {
		return 0
	}
	return 1
}
func newDtObstacleAvoidanceDebugData(maxSamples int) *dtObstacleAvoidanceDebugData {
	d := &dtObstacleAvoidanceDebugData{}
	d.m_maxSamples = maxSamples
	d.m_vel = make([]float64, 3*d.m_maxSamples)
	d.m_pen = make([]float64, d.m_maxSamples)
	d.m_ssize = make([]float64, d.m_maxSamples)
	d.m_vpen = make([]float64, d.m_maxSamples)
	d.m_vcpen = make([]float64, d.m_maxSamples)
	d.m_spen = make([]float64, d.m_maxSamples)
	d.m_tpen = make([]float64, d.m_maxSamples)
	return d
}
func (d *dtObstacleAvoidanceDebugData) reset() {
	d.m_nsamples = 0
}

func (d *dtObstacleAvoidanceDebugData) addSample(vel []float64, ssize float64, pen float64,
	vpen, vcpen, spen, tpen float64) {
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

func (d *dtObstacleAvoidanceDebugData) normalizeArray(arr []float64, n int) {
	// Normalize penaly range.
	minPen := math.MaxFloat64
	maxPen := math.SmallestNonzeroFloat64
	for i := 0; i < n; i++ {
		minPen = dtMin(minPen, arr[i])
		maxPen = dtMax(maxPen, arr[i])
	}
	penRange := maxPen - minPen
	s := float64(1)
	if penRange > 0.001 {
		s = 1.0 / penRange
	}
	for i := 0; i < n; i++ {
		arr[i] = dtClamp((arr[i]-minPen)*s, 0.0, 1.0)
	}

}

func (d *dtObstacleAvoidanceDebugData) normalizeSamples() {
	d.normalizeArray(d.m_pen, d.m_nsamples)
	d.normalizeArray(d.m_vpen, d.m_nsamples)
	d.normalizeArray(d.m_vcpen, d.m_nsamples)
	d.normalizeArray(d.m_spen, d.m_nsamples)
	d.normalizeArray(d.m_tpen, d.m_nsamples)
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

func (d *dtObstacleAvoidanceQuery) addCircle(pos []float64, rad float64,
	vel, dvel []float64) {
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

func (d *dtObstacleAvoidanceQuery) addSegment(p, q []float64) {
	if d.m_nsegments >= d.m_maxSegments {
		return
	}

	seg := d.m_segments[d.m_nsegments]
	d.m_nsegments++
	copy(seg.p[:], p)
	copy(seg.q[:], q)
}

func (d *dtObstacleAvoidanceQuery) prepare(pos, dvel []float64) {
	// Prepare obstacles
	for i := 0; i < d.m_ncircles; i++ {
		cir := d.m_circles[i]

		// Side
		pa := pos
		pb := cir.p

		orig := []float64{0, 0, 0}
		copy(cir.dp[:], dtVsub(pb[:], pa))
		dtVnormalize(cir.dp[:])
		dv := dtVsub(cir.dvel[:], dvel)

		a := dtTriArea2D(orig, cir.dp[:], dv)
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
		_, res := dtDistancePtSegSqr2D(pos, seg.p[:], seg.q[:])
		seg.touch = res < dtSqr(r)
	}
}

/* Calculate the collision penalty for a given velocity vector
 *
 * @param vcand sampled velocity
 * @param dvel desired velocity
 * @param minPenalty threshold penalty for early out
 */
func (d *dtObstacleAvoidanceQuery) processSample(vcand []float64, cs float64,
	pos []float64, rad float64,
	vel, dvel []float64,
	minPenalty float64,
	debug *dtObstacleAvoidanceDebugData) float64 {
	// penalty for straying away from the desired and current velocities
	vpen := d.m_params.weightDesVel * (dtVdist2D(vcand, dvel) * d.m_invVmax)
	vcpen := d.m_params.weightCurVel * (dtVdist2D(vcand, vel) * d.m_invVmax)

	// find the threshold hit time to bail out based on the early out penalty
	// (see how the penalty is calculated below to understand)
	minPen := minPenalty - vpen - vcpen
	tThresold := (d.m_params.weightToi/minPen - 0.1) * d.m_params.horizTime
	if tThresold-d.m_params.horizTime > math.MaxFloat32 {
		return minPenalty
	} // already too much
	// Find min time of impact and exit amongst all obstacles.
	tmin := d.m_params.horizTime
	side := float64(0)
	nside := 0

	for i := 0; i < d.m_ncircles; i++ {
		cir := d.m_circles[i]

		// RVO
		vab := make([]float64, 3)
		vab = dtVscale(vcand, 2)
		vab = dtVsub(vab, vel)
		vab = dtVsub(vab, cir.vel[:])

		// Side
		side += dtClamp(dtMin(dtVdot2D(cir.dp[:], vab)*0.5+0.5, dtVdot2D(cir.np[:], vab)*2), 0.0, 1.0)
		nside++

		var htmin = float64(0)
		var htmax = float64(0)
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
		htmin := float64(0)

		if seg.touch {
			// Special case when the agent is very close to the segment.
			sdir := make([]float64, 3)
			snorm := make([]float64, 3)
			sdir = dtVsub(seg.q[:], seg.p[:])
			snorm[0] = -sdir[2]
			snorm[2] = sdir[0]
			// If the velocity is pointing towards the segment, no collision.
			if dtVdot2D(snorm, vcand) < 0.0 {
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
		side /= float64(nside)
	}

	spen := d.m_params.weightSide * side
	tpen := d.m_params.weightToi * (1.0 / (0.1 + tmin*d.m_invHorizTime))

	penalty := vpen + vcpen + spen + tpen

	// Store different penalties for debug viewing
	if debug != nil {
		debug.addSample(vcand, cs, penalty, vpen, vcpen, spen, tpen)
	}

	return penalty
}

func (d *dtObstacleAvoidanceQuery) sampleVelocityGrid(pos []float64, rad float64, vmax float64, vel, dvel, nvel []float64,
	params *dtObstacleAvoidanceParams,
	debug *dtObstacleAvoidanceDebugData) int {
	d.prepare(pos, dvel)
	d.m_params = *params
	d.m_invHorizTime = 1.0 / d.m_params.horizTime
	d.m_vmax = vmax
	d.m_invVmax = math.MaxFloat64
	if vmax > 0 {
		d.m_invVmax = 1.0 / vmax
	}

	dtVset(nvel, 0, 0, 0)

	if debug != nil {
		debug.reset()
	}

	cvx := dvel[0] * d.m_params.velBias
	cvz := dvel[2] * d.m_params.velBias
	cs := vmax * 2 * (1 - d.m_params.velBias) / float64((d.m_params.gridSize - 1))
	half := float64(d.m_params.gridSize-1) * cs * 0.5

	minPenalty := math.MaxFloat64
	ns := 0

	for y := 0; y < d.m_params.gridSize; y++ {
		for x := 0; x < d.m_params.gridSize; x++ {
			vcand := make([]float64, 3)
			vcand[0] = cvx + float64(x)*cs - half
			vcand[1] = 0
			vcand[2] = cvz + float64(y)*cs - half

			if dtSqr(vcand[0])+dtSqr(vcand[2]) > dtSqr(vmax+cs/2) {
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
func dtNormalize2D(v []float64) {
	d := math.Sqrt(v[0]*v[0] + v[2]*v[2])
	if d == 0 {
		return
	}

	d = 1.0 / d
	v[0] *= d
	v[2] *= d
}

// vector normalization that ignores the y-component.
func dtRorate2D(dest, v []float64, ang float64) {
	c := math.Cos(ang)
	s := math.Sin(ang)
	dest[0] = v[0]*c - v[2]*s
	dest[2] = v[0]*s + v[2]*c
	dest[1] = v[1]
}

func (d *dtObstacleAvoidanceQuery) sampleVelocityAdaptive(pos []float64, rad, vmax float64,
	vel, dvel, nvel []float64, params *dtObstacleAvoidanceParams,
	debug *dtObstacleAvoidanceDebugData) int {
	d.prepare(pos, dvel)
	d.m_params = *params
	d.m_invHorizTime = 1.0 / d.m_params.horizTime
	d.m_vmax = vmax
	d.m_invVmax = math.MaxFloat64
	if vmax > 0 {
		d.m_invVmax = 1.0 / vmax
	}

	dtVset(nvel, 0, 0, 0)

	if debug != nil {
		debug.reset()
	}

	// Build sampling pattern aligned to desired velocity.
	pat := make([]float64, (DT_MAX_PATTERN_DIVS*DT_MAX_PATTERN_RINGS+1)*2)
	npat := 0

	ndivs := d.m_params.adaptiveDivs
	nrings := d.m_params.adaptiveRings
	depth := d.m_params.adaptiveDepth

	nd := dtClamp(ndivs, 1, DT_MAX_PATTERN_DIVS)
	nr := dtClamp(nrings, 1, DT_MAX_PATTERN_RINGS)
	da := (1.0 / float64(nd)) * DT_PI * 2
	ca := math.Cos(da)
	sa := math.Sin(da)

	// desired direction
	ddir := make([]float64, 6)
	copy(ddir, dvel)
	dtNormalize2D(ddir)
	dtRorate2D(ddir[3:], ddir, da*0.5) // rotated by da/2

	// Always add sample at zero
	pat[npat*2+0] = 0
	pat[npat*2+1] = 0
	npat++

	for j := 0; j < nr; j++ {
		r := float64((nr - j)) / float64(nr)
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
	cr := vmax * (1.0 - d.m_params.velBias)
	res := make([]float64, 3)
	dtVset(res, dvel[0]*d.m_params.velBias, 0, dvel[2]*d.m_params.velBias)
	ns := 0

	for k := 0; k < depth; k++ {
		minPenalty := math.MaxFloat64
		bvel := make([]float64, 3)
		dtVset(bvel, 0, 0, 0)

		for i := 0; i < npat; i++ {
			vcand := make([]float64, 3)
			vcand[0] = res[0] + pat[i*2+0]*cr
			vcand[1] = 0
			vcand[2] = res[2] + pat[i*2+1]*cr

			if dtSqr(vcand[0])+dtSqr(vcand[2]) > dtSqr(vmax+0.001) {
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

	return ns
}
