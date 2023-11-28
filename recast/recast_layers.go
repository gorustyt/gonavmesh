package recast

import "gonavamesh/common"

const RC_MAX_LAYERS_DEF = 63
const RC_MAX_NEIS_DEF = 16

// Keep type checking.
const RC_MAX_LAYERS = RC_MAX_LAYERS_DEF

const RC_MAX_NEIS = RC_MAX_NEIS_DEF

type rcLayerRegion struct {
	layers     [RC_MAX_LAYERS]int
	neis       [RC_MAX_NEIS]int
	ymin, ymax int
	layerId    int // Layer ID
	nlayers    int // Layer count
	nneis      int // Neighbour count
	base       int // Flag indicating if the region is the base of merged regions.
}

func contains(a []int, an int, v int) bool {
	n := an
	for i := 0; i < n; i++ {
		if a[i] == v {
			return true
		}

	}
	return false
}

func addUnique(a []int, an *int, anMax int, v int) bool {
	if contains(a, *an, v) {
		return true
	}

	if *an >= anMax {
		return false
	}

	a[*an] = v
	*an++
	return true
}

func LayersOverlapRange(amin, amax, bmin, bmax int) bool {
	if amin > bmax || amax < bmin {
		return false
	}
	return true
}

type rcLayerSweepSpan struct {
	ns  int // number samples
	id  int // region id
	nei int // neighbour id
}

// / Represents a set of heightfield layers.
// / @ingroup recast
// / @see rcAllocHeightfieldLayerSet, rcFreeHeightfieldLayerSet
type RcHeightfieldLayerSet struct {
	Layers  []*RcHeightfieldLayer ///< The layers in the set. [Size: #nlayers]
	Nlayers int                   ///< The number of layers in the set.
}

// / Represents a heightfield layer within a layer set.
// / @see RcHeightfieldLayerSet
type RcHeightfieldLayer struct {
	Bmin    [3]float64 ///< The minimum bounds in world space. [(x, y, z)]
	Bmax    [3]float64 ///< The maximum bounds in world space. [(x, y, z)]
	Cs      float64    ///< The size of each cell. (On the xz-plane.)
	Ch      float64    ///< The height of each cell. (The minimum increment along the y-axis.)
	Width   int        ///< The width of the heightfield. (Along the x-axis in cell units.)
	Height  int        ///< The height of the heightfield. (Along the z-axis in cell units.)
	Minx    int        ///< The minimum x-bounds of usable data.
	Maxx    int        ///< The maximum x-bounds of usable data.
	Miny    int        ///< The minimum y-bounds of usable data. (Along the z-axis.)
	Maxy    int        ///< The maximum y-bounds of usable data. (Along the z-axis.)
	Hmin    int        ///< The minimum height bounds of usable data. (Along the y-axis.)
	Hmax    int        ///< The maximum height bounds of usable data. (Along the y-axis.)
	Heights []int      ///< The heightfield. [Size: width * height]
	Areas   []int      ///< Area ids. [Size: Same as #heights]
	Cons    []int      ///< Packed neighbor connection information. [Size: Same as #heights]
}

// / @par
// /
// / See the #RcConfig documentation for more information on the configuration parameters.
// /
// / @see rcAllocHeightfieldLayerSet, RcCompactHeightfield, RcHeightfieldLayerSet, RcConfig
func RcBuildHeightfieldLayers(chf *RcCompactHeightfield, borderSize, walkableHeight int, lset *RcHeightfieldLayerSet) bool {

	w := chf.Width
	h := chf.Height
	srcReg := make([]int, chf.SpanCount)
	for i := range srcReg {
		srcReg[i] = 0xff
	}

	nsweeps := chf.Width
	sweeps := make([]*rcLayerSweepSpan, nsweeps)

	// Partition walkable area into monotone regions.
	prevCount := make([]int, 256)
	regId := 0

	for y := borderSize; y < h-borderSize; y++ {
		sweepId := 0

		for x := borderSize; x < w-borderSize; x++ {
			c := chf.Cells[x+y*w]

			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni; i++ {
				s := chf.Spans[i]
				if chf.Areas[i] == RC_NULL_AREA {
					continue
				}
				sid := 0xff

				// -x
				if rcGetCon(s, 0) != RC_NOT_CONNECTED {
					ax := x + common.GetDirOffsetX(0)
					ay := y + common.GetDirOffsetY(0)
					ai := chf.Cells[ax+ay*w].Index + rcGetCon(s, 0)
					if chf.Areas[ai] != RC_NULL_AREA && srcReg[ai] != 0xff {
						sid = srcReg[ai]
					}

				}

				if sid == 0xff {
					sid = sweepId
					sweepId++
					sweeps[sid].nei = 0xff
					sweeps[sid].ns = 0
				}

				// -y
				if rcGetCon(s, 3) != RC_NOT_CONNECTED {
					ax := x + common.GetDirOffsetX(3)
					ay := y + common.GetDirOffsetY(3)
					ai := chf.Cells[ax+ay*w].Index + rcGetCon(s, 3)
					nr := srcReg[ai]
					if nr != 0xff {
						// Set neighbour when first valid neighbour is encoutered.
						if sweeps[sid].ns == 0 {
							sweeps[sid].nei = nr
						}

						if sweeps[sid].nei == nr {
							// Update existing neighbour
							sweeps[sid].ns++
							prevCount[nr]++
						} else {
							// This is hit if there is nore than one neighbour.
							// Invalidate the neighbour.
							sweeps[sid].nei = 0xff
						}
					}
				}

				srcReg[i] = sid
			}
		}

		// Create unique ID.
		for i := 0; i < sweepId; i++ {
			// If the neighbour is set and there is only one continuous connection to it,
			// the sweep will be merged with the previous one, else new region is created.
			if sweeps[i].nei != 0xff && prevCount[sweeps[i].nei] == sweeps[i].ns {
				sweeps[i].id = sweeps[i].nei
			} else {
				if regId == 255 {
					return false
				}
				sweeps[i].id = regId
				regId++
			}
		}

		// Remap local sweep ids to region ids.
		for x := borderSize; x < w-borderSize; x++ {
			c := chf.Cells[x+y*w]
			i := c.Index
			ni := (int)(c.Index + c.Count)
			for ; i < ni; i++ {
				if srcReg[i] != 0xff {
					srcReg[i] = sweeps[srcReg[i]].id
				}

			}
		}
	}

	// Allocate and init layer regions.
	nregs := regId
	regs := make([]*rcLayerRegion, nregs)
	for i := range regs {
		regs[i] = &rcLayerRegion{}
	}
	for i := 0; i < nregs; i++ {
		regs[i].layerId = 0xff
		regs[i].ymin = 0xffff
		regs[i].ymax = 0
	}

	// Find region neighbours and overlapping regions.
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			c := chf.Cells[x+y*w]

			lregs := make([]int, RC_MAX_LAYERS)
			nlregs := 0
			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni; i++ {
				s := chf.Spans[i]
				ri := srcReg[i]
				if ri == 0xff {
					continue
				}

				regs[ri].ymin = common.Min(regs[ri].ymin, s.Y)
				regs[ri].ymax = common.Max(regs[ri].ymax, s.Y)

				// Collect all region layers.
				if nlregs < RC_MAX_LAYERS {
					lregs[nlregs] = ri
					nlregs++
				}
				// Update neighbours
				for dir := 0; dir < 4; dir++ {
					if rcGetCon(s, dir) != RC_NOT_CONNECTED {
						ax := x + common.GetDirOffsetX(dir)
						ay := y + common.GetDirOffsetY(dir)
						ai := chf.Cells[ax+ay*w].Index + rcGetCon(s, dir)
						rai := srcReg[ai]
						if rai != 0xff && rai != ri {
							// Don't check return value -- if we cannot add the neighbor
							// it will just cause a few more regions to be created, which
							// is fine.
							addUnique(regs[ri].neis[:], &regs[ri].nneis, RC_MAX_NEIS, rai)
						}
					}
				}
			}

			// Update overlapping regions.
			for i := 0; i < nlregs-1; i++ {
				for j := i + 1; j < nlregs; j++ {
					if lregs[i] != lregs[j] {
						ri := regs[lregs[i]]
						rj := regs[lregs[j]]

						if !addUnique(ri.layers[:], &ri.nlayers, RC_MAX_LAYERS, lregs[j]) ||
							!addUnique(rj.layers[:], &rj.nlayers, RC_MAX_LAYERS, lregs[i]) {
							return false
						}
					}
				}
			}
		}

	}

	// Create 2D layers from regions.
	layerId := 0

	MAX_STACK := 64
	stack := make([]int, MAX_STACK)
	nstack := 0

	for i := 0; i < nregs; i++ {
		root := regs[i]
		// Skip already visited.
		if root.layerId != 0xff {
			continue
		}

		// Start search.
		root.layerId = layerId
		root.base = 1

		nstack = 0
		stack[nstack] = i
		nstack++

		for nstack > 0 {
			// Pop front
			reg := regs[stack[0]]
			nstack--
			for j := 0; j < nstack; j++ {
				stack[j] = stack[j+1]
			}

			nneis := reg.nneis
			for j := 0; j < nneis; j++ {
				nei := reg.neis[j]
				regn := regs[nei]
				// Skip already visited.
				if regn.layerId != 0xff {
					continue
				}

				// Skip if the neighbour is overlapping root region.
				if contains(root.layers[:], root.nlayers, nei) {
					continue
				}

				// Skip if the height range would become too large.
				ymin := common.Min(root.ymin, regn.ymin)
				ymax := common.Max(root.ymax, regn.ymax)
				if (ymax - ymin) >= 255 {
					continue
				}

				if nstack < MAX_STACK {
					// Deepen
					stack[nstack] = nei
					nstack++

					// Mark layer id
					regn.layerId = layerId
					// Merge current layers to root.
					for k := 0; k < regn.nlayers; k++ {
						if !addUnique(root.layers[:], &root.nlayers, RC_MAX_LAYERS, regn.layers[k]) {
							return false
						}
					}
					root.ymin = common.Min(root.ymin, regn.ymin)
					root.ymax = common.Max(root.ymax, regn.ymax)
				}
			}
		}

		layerId++
	}

	// Merge non-overlapping regions that are close in height.
	mergeHeight := walkableHeight * 4

	for i := 0; i < nregs; i++ {
		ri := regs[i]
		if ri.base == 0 {
			continue
		}

		newId := ri.layerId

		for {
			oldId := 0xff

			for j := 0; j < nregs; j++ {
				if i == j {
					continue
				}
				rj := regs[j]
				if rj.base == 0 {
					continue
				}

				// Skip if the regions are not close to each other.
				if !LayersOverlapRange(ri.ymin, ri.ymax+mergeHeight, rj.ymin, rj.ymax+mergeHeight) {
					continue
				}

				// Skip if the height range would become too large.
				ymin := common.Min(ri.ymin, rj.ymin)
				ymax := common.Max(ri.ymax, rj.ymax)
				if (ymax - ymin) >= 255 {
					continue
				}

				// Make sure that there is no overlap when merging 'ri' and 'rj'.
				overlap := false
				// Iterate over all regions which have the same layerId as 'rj'
				for k := 0; k < nregs; k++ {
					if regs[k].layerId != rj.layerId {
						continue
					}

					// Check if region 'k' is overlapping region 'ri'
					// Index to 'regs' is the same as region id.
					if contains(ri.layers[:], ri.nlayers, k) {
						overlap = true
						break
					}
				}
				// Cannot merge of regions overlap.
				if overlap {
					continue
				}

				// Can merge i and j.
				oldId = rj.layerId
				break
			}

			// Could not find anything to merge with, stop.
			if oldId == 0xff {
				break
			}

			// Merge
			for j := 0; j < nregs; j++ {
				rj := regs[j]
				if rj.layerId == oldId {
					rj.base = 0
					// Remap layerIds.
					rj.layerId = newId
					// Add overlaid layers from 'rj' to 'ri'.
					for k := 0; k < rj.nlayers; k++ {
						if !addUnique(ri.layers[:], &ri.nlayers, RC_MAX_LAYERS, rj.layers[k]) {
							return false
						}
					}

					// Update height bounds.
					ri.ymin = common.Min(ri.ymin, rj.ymin)
					ri.ymax = common.Max(ri.ymax, rj.ymax)
				}
			}
		}
	}

	// Compact layerIds
	remap := make([]int, 256)

	// Find number of unique layers.
	layerId = 0
	for i := 0; i < nregs; i++ {
		remap[regs[i].layerId] = 1
	}

	for i := 0; i < 256; i++ {
		if remap[i] > 0 {
			remap[i] = layerId
			layerId++
		} else {
			remap[i] = 0xff
		}

	}
	// Remap ids.
	for i := 0; i < nregs; i++ {
		regs[i].layerId = remap[regs[i].layerId]
	}

	// No layers, return empty.
	if layerId == 0 {
		return true
	}

	// Create layers.
	lw := w - borderSize*2
	lh := h - borderSize*2

	// Build contracted bbox for layers.
	var bmin, bmax [3]float64
	copy(bmin[:], chf.Bmin[:])
	copy(bmax[:], chf.Bmax[:])
	bmin[0] += float64(borderSize) * chf.Cs
	bmin[2] += float64(borderSize) * chf.Cs
	bmax[0] -= float64(borderSize) * chf.Cs
	bmax[2] -= float64(borderSize) * chf.Cs

	lset.Nlayers = layerId

	lset.Layers = make([]*RcHeightfieldLayer, lset.Nlayers)
	for i := range lset.Layers {
		lset.Layers[i] = &RcHeightfieldLayer{}
	}

	// Store layers.
	for i := 0; i < lset.Nlayers; i++ {
		curId := i

		layer := lset.Layers[i]

		gridSize := lw * lh

		layer.Heights = make([]int, gridSize)
		for i := range layer.Heights {
			layer.Heights[i] = 0xff
		}

		layer.Areas = make([]int, gridSize)
		layer.Cons = make([]int, gridSize)
		// Find layer height bounds.
		hmin := 0
		hmax := 0
		for j := 0; j < nregs; j++ {
			if regs[j].base > 0 && regs[j].layerId == curId {
				hmin = regs[j].ymin
				hmax = regs[j].ymax
			}
		}

		layer.Width = lw
		layer.Height = lh
		layer.Cs = chf.Cs
		layer.Ch = chf.Ch

		// Adjust the bbox to fit the heightfield.
		copy(layer.Bmin[:], bmin[:])
		copy(layer.Bmax[:], bmax[:])
		layer.Bmin[1] = bmin[1] + float64(hmin)*chf.Ch
		layer.Bmax[1] = bmin[1] + float64(hmax)*chf.Ch
		layer.Hmin = hmin
		layer.Hmax = hmax

		// Update usable data region.
		layer.Minx = layer.Width
		layer.Maxx = 0
		layer.Miny = layer.Height
		layer.Maxy = 0

		// Copy height and area from compact heightfield.
		for y := 0; y < lh; y++ {
			for x := 0; x < lw; x++ {
				cx := borderSize + x
				cy := borderSize + y
				c := chf.Cells[cx+cy*w]
				j := c.Index
				nj := (c.Index + c.Count)
				for ; j < nj; j++ {
					s := chf.Spans[j]
					// Skip unassigned regions.
					if srcReg[j] == 0xff {
						continue
					}

					// Skip of does nto belong to current layer.
					lid := regs[srcReg[j]].layerId
					if lid != curId {
						continue
					}

					// Update data bounds.
					layer.Minx = common.Min(layer.Minx, x)
					layer.Maxx = common.Max(layer.Maxx, x)
					layer.Miny = common.Min(layer.Miny, y)
					layer.Maxy = common.Max(layer.Maxy, y)

					// Store height and area type.
					idx := x + y*lw
					layer.Heights[idx] = (s.Y - hmin)
					layer.Areas[idx] = chf.Areas[j]

					// Check connection.
					portal := 0
					con := 0
					for dir := 0; dir < 4; dir++ {
						if rcGetCon(s, dir) != RC_NOT_CONNECTED {
							ax := cx + common.GetDirOffsetX(dir)
							ay := cy + common.GetDirOffsetY(dir)
							ai := chf.Cells[ax+ay*w].Index + rcGetCon(s, dir)
							alid := 0xff
							if srcReg[ai] != 0xff {
								alid = regs[srcReg[ai]].layerId
							}
							// Portal mask
							if chf.Areas[ai] != RC_NULL_AREA && lid != alid {
								portal |= (1 << dir)
								// Update height so that it matches on both sides of the portal.
								as := chf.Spans[ai]
								if as.Y > hmin {
									layer.Heights[idx] = common.Max(layer.Heights[idx], (as.Y - hmin))
								}

							}
							// Valid connection mask
							if chf.Areas[ai] != RC_NULL_AREA && lid == alid {
								nx := ax - borderSize
								ny := ay - borderSize
								if nx >= 0 && ny >= 0 && nx < lw && ny < lh {
									con |= (1 << dir)
								}

							}
						}
					}

					layer.Cons[idx] = (portal << 4) | con
				}
			}
		}

		if layer.Minx > layer.Maxx {
			layer.Maxx = 0
			layer.Minx = layer.Maxx
		}

		if layer.Miny > layer.Maxy {
			layer.Maxy = 0
			layer.Miny = layer.Maxy
		}

	}

	return true
}
