package recast

import (
	"fmt"
	"reflect"
)

const (
	/// Heightfield border flag.
	/// If a heightfield region ID has this bit set, then the region is a border
	/// region and its spans are considered un-walkable.
	/// (Used during the region and contour build process.)
	/// @see rcCompactSpan::reg
	RC_BORDER_REG = 0x8000
)

func calculateDistanceField(chf *rcCompactHeightfield, src []int, maxDist *int) {
	w := chf.width
	h := chf.height

	// Init distance and points.
	for i := 0; i < chf.spanCount; i++ {
		src[i] = 0xffff
	}

	// Mark boundary cells.
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			c := chf.cells[x+y*w]
			i := c.index
			ni := (c.index + c.count)
			for ; i < ni; i++ {
				s := chf.spans[i]
				area := chf.areas[i]

				nc := 0
				for dir := 0; dir < 4; dir++ {
					if rcGetCon(s, dir) != RC_NOT_CONNECTED {
						ax := x + rcGetDirOffsetX(dir)
						ay := y + rcGetDirOffsetY(dir)
						ai := chf.cells[ax+ay*w].index + rcGetCon(s, dir)
						if area == chf.areas[ai] {
							nc++
						}

					}
				}
				if nc != 4 {
					src[i] = 0
				}

			}
		}
	}

	// Pass 1
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			c := chf.cells[x+y*w]
			i := c.index
			ni := (c.index + c.count)
			for ; i < ni; i++ {
				s := chf.spans[i]

				if rcGetCon(s, 0) != RC_NOT_CONNECTED {
					// (-1,0)
					ax := x + rcGetDirOffsetX(0)
					ay := y + rcGetDirOffsetY(0)
					ai := chf.cells[ax+ay*w].index + rcGetCon(s, 0)
					as := chf.spans[ai]
					if src[ai]+2 < src[i] {
						src[i] = src[ai] + 2
					}

					// (-1,-1)
					if rcGetCon(as, 3) != RC_NOT_CONNECTED {
						aax := ax + rcGetDirOffsetX(3)
						aay := ay + rcGetDirOffsetY(3)
						aai := chf.cells[aax+aay*w].index + rcGetCon(as, 3)
						if src[aai]+3 < src[i] {
							src[i] = src[aai] + 3
						}

					}
				}
				if rcGetCon(s, 3) != RC_NOT_CONNECTED {
					// (0,-1)
					ax := x + rcGetDirOffsetX(3)
					ay := y + rcGetDirOffsetY(3)
					ai := chf.cells[ax+ay*w].index + rcGetCon(s, 3)
					as := chf.spans[ai]
					if src[ai]+2 < src[i] {
						src[i] = src[ai] + 2
					}

					// (1,-1)
					if rcGetCon(as, 2) != RC_NOT_CONNECTED {
						aax := ax + rcGetDirOffsetX(2)
						aay := ay + rcGetDirOffsetY(2)
						aai := chf.cells[aax+aay*w].index + rcGetCon(as, 2)
						if src[aai]+3 < src[i] {
							src[i] = src[aai] + 3
						}

					}
				}
			}
		}
	}

	// Pass 2
	for y := h - 1; y >= 0; y-- {
		for x := w - 1; x >= 0; x-- {
			c := chf.cells[x+y*w]
			i := c.index
			ni := (c.index + c.count)
			for ; i < ni; i++ {
				s := chf.spans[i]

				if rcGetCon(s, 2) != RC_NOT_CONNECTED {
					// (1,0)
					ax := x + rcGetDirOffsetX(2)
					ay := y + rcGetDirOffsetY(2)
					ai := chf.cells[ax+ay*w].index + rcGetCon(s, 2)
					as := chf.spans[ai]
					if src[ai]+2 < src[i] {
						src[i] = src[ai] + 2
					}

					// (1,1)
					if rcGetCon(as, 1) != RC_NOT_CONNECTED {
						aax := ax + rcGetDirOffsetX(1)
						aay := ay + rcGetDirOffsetY(1)
						aai := chf.cells[aax+aay*w].index + rcGetCon(as, 1)
						if src[aai]+3 < src[i] {
							src[i] = src[aai] + 3
						}

					}
				}
				if rcGetCon(s, 1) != RC_NOT_CONNECTED {
					// (0,1)
					ax := x + rcGetDirOffsetX(1)
					ay := y + rcGetDirOffsetY(1)
					ai := chf.cells[ax+ay*w].index + rcGetCon(s, 1)
					as := chf.spans[ai]
					if src[ai]+2 < src[i] {
						src[i] = src[ai] + 2
					}

					// (-1,1)
					if rcGetCon(as, 0) != RC_NOT_CONNECTED {
						aax := ax + rcGetDirOffsetX(0)
						aay := ay + rcGetDirOffsetY(0)
						aai := chf.cells[aax+aay*w].index + rcGetCon(as, 0)
						if src[aai]+3 < src[i] {
							src[i] = src[aai] + 3
						}

					}
				}
			}
		}
	}

	*maxDist = 0
	for i := 0; i < chf.spanCount; i++ {
		*maxDist = rcMax(src[i], *maxDist)
	}

}

func boxBlur(chf *rcCompactHeightfield, thr int, src, dst []int) []int {
	w := chf.width
	h := chf.height

	thr *= 2

	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			c := chf.cells[x+y*w]
			i := c.index
			ni := (c.index + c.count)
			for ; i < ni; i++ {
				s := chf.spans[i]
				cd := src[i]
				if cd <= thr {
					dst[i] = cd
					continue
				}

				d := cd
				for dir := 0; dir < 4; dir++ {
					if rcGetCon(s, dir) != RC_NOT_CONNECTED {
						ax := x + rcGetDirOffsetX(dir)
						ay := y + rcGetDirOffsetY(dir)
						ai := chf.cells[ax+ay*w].index + rcGetCon(s, dir)
						d += src[ai]

						as := chf.spans[ai]
						dir2 := (dir + 1) & 0x3
						if rcGetCon(as, dir2) != RC_NOT_CONNECTED {
							ax2 := ax + rcGetDirOffsetX(dir2)
							ay2 := ay + rcGetDirOffsetY(dir2)
							ai2 := chf.cells[ax2+ay2*w].index + rcGetCon(as, dir2)
							d += src[ai2]
						} else {
							d += cd
						}
					} else {
						d += cd * 2
					}
				}
				dst[i] = ((d + 5) / 9)
			}
		}
	}
	return dst
}

type LevelStackEntry struct {
	x     int
	y     int
	index int
}

func newLevelStackEntry(x int, y int, index int) *LevelStackEntry {
	return &LevelStackEntry{
		x:     x,
		y:     y,
		index: index,
	}
}
func floodRegion(x, y, i int, level, r int,
	chf *rcCompactHeightfield, srcReg, srcDist []int, stack Stack[*LevelStackEntry]) bool {
	w := chf.width

	area := chf.areas[i]

	// Flood fill mark region.
	stack.Clear()
	stack.Push(newLevelStackEntry(x, y, i))
	srcReg[i] = r
	srcDist[i] = 0

	lev := 0
	if level >= 2 {
		lev = level - 2
	}
	count := 0

	for !stack.Empty() {
		back := stack.Pop()
		cx := back.x
		cy := back.y
		ci := back.index
		stack.Pop()

		cs := chf.spans[ci]

		// Check if any of the neighbours already have a valid region set.
		ar := 0
		for dir := 0; dir < 4; dir++ {
			// 8 connected
			if rcGetCon(cs, dir) != RC_NOT_CONNECTED {
				ax := cx + rcGetDirOffsetX(dir)
				ay := cy + rcGetDirOffsetY(dir)
				ai := chf.cells[ax+ay*w].index + rcGetCon(cs, dir)
				if chf.areas[ai] != area {
					continue
				}

				nr := srcReg[ai]
				if nr&RC_BORDER_REG > 0 {
					continue
				} // Do not take borders into account.

				if nr != 0 && nr != r {
					ar = nr
					break
				}

				as := chf.spans[ai]

				dir2 := (dir + 1) & 0x3
				if rcGetCon(as, dir2) != RC_NOT_CONNECTED {
					ax2 := ax + rcGetDirOffsetX(dir2)
					ay2 := ay + rcGetDirOffsetY(dir2)
					ai2 := chf.cells[ax2+ay2*w].index + rcGetCon(as, dir2)
					if chf.areas[ai2] != area {
						continue
					}

					nr2 := srcReg[ai2]
					if nr2 != 0 && nr2 != r {
						ar = nr2
						break
					}
				}
			}
		}
		if ar != 0 {
			srcReg[ci] = 0
			continue
		}

		count++

		// Expand neighbours.
		for dir := 0; dir < 4; dir++ {
			if rcGetCon(cs, dir) != RC_NOT_CONNECTED {
				ax := cx + rcGetDirOffsetX(dir)
				ay := cy + rcGetDirOffsetY(dir)
				ai := chf.cells[ax+ay*w].index + rcGetCon(cs, dir)
				if chf.areas[ai] != area {
					continue
				}

				if chf.dist[ai] >= lev && srcReg[ai] == 0 {
					srcReg[ai] = r
					srcDist[ai] = 0
					stack.Push(newLevelStackEntry(ax, ay, ai))
				}
			}
		}
	}

	return count > 0
}

// Struct to keep track of entries in the region table that have been changed.
type DirtyEntry struct {
	index     int
	region    int
	distance2 int
}

func newDirtyEntry(index int, region int, distance2 int) *DirtyEntry {
	return &DirtyEntry{index: index, region: region, distance2: distance2}
}

func expandRegions(maxIter, level int,
	chf *rcCompactHeightfield,
	srcReg, srcDist []int,
	stack Stack[*LevelStackEntry],
	fillStack bool) {
	w := chf.width
	h := chf.height

	if fillStack {
		// Find cells revealed by the raised level.
		stack.Clear()
		for y := 0; y < h; y++ {
			for x := 0; x < w; x++ {
				c := chf.cells[x+y*w]
				i := c.index
				ni := (c.index + c.count)
				for ; i < ni; i++ {
					if chf.dist[i] >= level && srcReg[i] == 0 && chf.areas[i] != RC_NULL_AREA {
						stack.Push(newLevelStackEntry(x, y, i))
					}
				}
			}
		}
	} else { // use cells in the input stack

		// mark all cells which already have a region
		for j := 0; j < stack.Len(); j++ {
			i := stack.Index(j).index
			if srcReg[i] != 0 {
				stack.Index(j).index = -1
			}

		}
	}

	var dirtyEntries = NewStack[*DirtyEntry](func() *DirtyEntry {
		return newDirtyEntry(0, 0, 0)
	})
	iter := 0
	for stack.Len() > 0 {
		failed := 0
		dirtyEntries.Clear()

		for j := 0; j < stack.Len(); j++ {
			x := stack.Index(j).x
			y := stack.Index(j).y
			i := stack.Index(j).index
			if i < 0 {
				failed++
				continue
			}

			r := srcReg[i]
			d2 := 0xffff
			area := chf.areas[i]
			s := chf.spans[i]
			for dir := 0; dir < 4; dir++ {
				if rcGetCon(s, dir) == RC_NOT_CONNECTED {
					continue
				}
				ax := x + rcGetDirOffsetX(dir)
				ay := y + rcGetDirOffsetY(dir)
				ai := chf.cells[ax+ay*w].index + rcGetCon(s, dir)
				if chf.areas[ai] != area {
					continue
				}
				if srcReg[ai] > 0 && (srcReg[ai]&RC_BORDER_REG) == 0 {
					if srcDist[ai]+2 < d2 {
						r = srcReg[ai]
						d2 = srcDist[ai] + 2
					}
				}
			}
			if r > 0 {
				stack.Index(j).index = -1 // mark as used
				dirtyEntries.Push(newDirtyEntry(i, r, d2))
			} else {
				failed++
			}
		}

		// Copy entries that differ between src and dst to keep them in sync.
		for i := 0; i < dirtyEntries.Len(); i++ {
			idx := dirtyEntries.Index(i).index
			srcReg[idx] = dirtyEntries.Index(i).region
			srcDist[idx] = dirtyEntries.Index(i).distance2
		}

		if failed == stack.Len() {
			break
		}

		if level > 0 {
			iter++
			if iter >= maxIter {
				break
			}

		}
	}
}

func sortCellsByLevel(startLevel int,
	chf *rcCompactHeightfield,
	srcReg []int,
	nbStacks int, stacks []Stack[*LevelStackEntry],
	loglevelsPerStack int) { // the levels per stack (2 in our case) as a bit shift

	w := chf.width
	h := chf.height
	startLevel = startLevel >> loglevelsPerStack

	for j := 0; j < nbStacks; j++ {
		stacks[j].Clear()
	}

	// put all cells in the level range into the appropriate stacks
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			c := chf.cells[x+y*w]
			i := c.index
			ni := (c.index + c.count)
			for ; i < ni; i++ {
				if chf.areas[i] == RC_NULL_AREA || srcReg[i] != 0 {
					continue
				}

				level := chf.dist[i] >> loglevelsPerStack
				sId := startLevel - level
				if sId >= nbStacks {
					continue
				}

				if sId < 0 {
					sId = 0
				}

				stacks[sId].Push(newLevelStackEntry(x, y, i))
			}
		}
	}
}

func appendStacks(srcStack Stack[*LevelStackEntry],
	dstStack Stack[*LevelStackEntry],
	srcReg []int) {
	for j := 0; j < srcStack.Len(); j++ {
		i := srcStack.Index(j).index
		if (i < 0) || (srcReg[i] != 0) {
			continue
		}

		dstStack.Push(srcStack.Index(j))
	}
}

type rcRegion[T any] struct {
	spanCount        int // Number of spans belonging to this region
	id               int // ID of the region
	areaType         int // Are type.
	remap            bool
	visited          bool
	overlap          bool
	connectsToBorder bool
	ymin, ymax       int
	connections      Stack[T]
	floors           Stack[T]
}

func newRcRegion[T any](i int) *rcRegion[T] {
	return &rcRegion[T]{
		id:   i,
		ymin: 0xffff,
	}
}

func removeAdjacentNeighbours(reg *rcRegion[int]) {
	// Remove adjacent duplicates.
	for i := 0; i < reg.connections.Len() && reg.connections.Len() > 1; {
		ni := (i + 1) % reg.connections.Len()
		if reg.connections.Index(i) == reg.connections.Index(ni) {
			// Remove duplicate
			for j := i; j < reg.connections.Len()-1; j++ {
				reg.connections.SetByIndex(j, reg.connections.Index(j+1))
			}

			reg.connections.Pop()
		} else {
			i++
		}

	}
}

func replaceNeighbour(reg *rcRegion[int], oldId, newId int) {
	neiChanged := false
	for i := 0; i < reg.connections.Len(); i++ {
		if reg.connections.Index(i) == oldId {
			reg.connections.SetByIndex(i, newId)
			neiChanged = true
		}
	}
	for i := 0; i < reg.floors.Len(); i++ {
		if reg.floors.Index(i) == oldId {
			reg.floors.SetByIndex(i, newId)
		}
	}
	if neiChanged {
		removeAdjacentNeighbours(reg)
	}

}

func canMergeWithRegion(rega *rcRegion[int], regb *rcRegion[int]) bool {
	if rega.areaType != regb.areaType {
		return false
	}

	n := 0
	for i := 0; i < rega.connections.Len(); i++ {
		if rega.connections.Index(i) == regb.id {
			n++
		}

	}
	if n > 1 {
		return false
	}

	for i := 0; i < rega.floors.Len(); i++ {
		if rega.floors.Index(i) == regb.id {
			return false
		}

	}
	return true
}

func addUniqueFloorRegion(reg *rcRegion[int], n int) {
	for i := 0; i < reg.floors.Len(); i++ {
		if reg.floors.Index(i) == n {
			return
		}

	}

	reg.floors.Push(n)
}

func mergeRegions(rega, regb *rcRegion[int]) bool {
	aid := rega.id
	bid := regb.id

	// Duplicate current neighbourhood.
	acon := NewStack[int](func() int {
		return 0
	})
	acon.Resize(rega.connections.Len())
	for i := 0; i < rega.connections.Len(); i++ {
		acon.SetByIndex(i, rega.connections.Index(i))
	}

	bcon := regb.connections

	// Find insertion point on A.
	insa := -1
	for i := 0; i < acon.Len(); i++ {
		if acon.Index(i) == bid {
			insa = i
			break
		}
	}
	if insa == -1 {
		return false
	}

	// Find insertion point on B.
	insb := -1
	for i := 0; i < bcon.Len(); i++ {
		if bcon.Index(i) == aid {
			insb = i
			break
		}
	}
	if insb == -1 {
		return false
	}

	// Merge neighbours.
	rega.connections.Clear()
	i := 0
	ni := acon.Len()
	for ; i < ni-1; i++ {
		rega.connections.Push(acon.Index((insa + 1 + i) % ni))
	}

	i = 0
	ni = bcon.Len()
	for ; i < ni-1; i++ {
		rega.connections.Push(bcon.Index((insb + 1 + i) % ni))
	}

	removeAdjacentNeighbours(rega)

	for j := 0; j < regb.floors.Len(); j++ {
		addUniqueFloorRegion(rega, regb.floors.Index(j))
	}

	rega.spanCount += regb.spanCount
	regb.spanCount = 0
	regb.connections.Resize(0)

	return true
}

func isRegionConnectedToBorder(reg *rcRegion[int]) bool {
	// Region is connected to border if
	// one of the neighbours is null id.
	for i := 0; i < reg.connections.Len(); i++ {
		if reg.connections.Index(i) == 0 {
			return true
		}

	}
	return false
}

func isSolidEdge(chf *rcCompactHeightfield, srcReg []int, x, y, i, dir int) bool {
	s := chf.spans[i]
	r := 0
	if rcGetCon(s, dir) != RC_NOT_CONNECTED {
		ax := x + rcGetDirOffsetX(dir)
		ay := y + rcGetDirOffsetY(dir)
		ai := chf.cells[ax+ay*chf.width].index + rcGetCon(s, dir)
		r = srcReg[ai]
	}
	if r == srcReg[i] {
		return false
	}

	return true
}

func regionWalkContour(x, y, i, dir int, chf *rcCompactHeightfield, srcReg []int,
	cont Stack[int]) {
	startDir := dir
	starti := i

	ss := chf.spans[i]
	curReg := 0
	if rcGetCon(ss, dir) != RC_NOT_CONNECTED {
		ax := x + rcGetDirOffsetX(dir)
		ay := y + rcGetDirOffsetY(dir)
		ai := chf.cells[ax+ay*chf.width].index + rcGetCon(ss, dir)
		curReg = srcReg[ai]
	}
	cont.Push(curReg)

	iter := 0
	for {
		iter++
		if iter >= 40000 {
			break
		}
		s := chf.spans[i]

		if isSolidEdge(chf, srcReg, x, y, i, dir) {
			// Choose the edge corner
			r := 0
			if rcGetCon(s, dir) != RC_NOT_CONNECTED {
				ax := x + rcGetDirOffsetX(dir)
				ay := y + rcGetDirOffsetY(dir)
				ai := chf.cells[ax+ay*chf.width].index + rcGetCon(s, dir)
				r = srcReg[ai]
			}
			if r != curReg {
				curReg = r
				cont.Push(curReg)
			}

			dir = (dir + 1) & 0x3 // Rotate CW
		} else {
			ni := -1
			nx := x + rcGetDirOffsetX(dir)
			ny := y + rcGetDirOffsetY(dir)
			if rcGetCon(s, dir) != RC_NOT_CONNECTED {
				nc := chf.cells[nx+ny*chf.width]
				ni = nc.index + rcGetCon(s, dir)
			}
			if ni == -1 {
				// Should not happen.
				return
			}
			x = nx
			y = ny
			i = ni
			dir = (dir + 3) & 0x3 // Rotate CCW
		}

		if starti == i && startDir == dir {
			break
		}
	}

	// Remove adjacent duplicates.
	if cont.Len() > 1 {
		for j := 0; j < cont.Len(); {
			nj := (j + 1) % cont.Len()
			if cont.Index(j) == cont.Index(nj) {
				for k := j; k < cont.Len()-1; k++ {
					cont.SetByIndex(k, cont.Index(k+1))
				}

				cont.Pop()
			} else {
				j++
			}

		}
	}
}

func addUniqueConnection(reg *rcRegion[int], n int) {
	for i := 0; i < reg.connections.Len(); i++ {
		if reg.connections.Index(i) == n {
			return
		}

	}

	reg.connections.Push(n)
}
func mergeAndFilterLayerRegions(minRegionArea int, maxRegionId *int, chf *rcCompactHeightfield,
	srcReg []int) bool {
	w := chf.width
	h := chf.height

	nreg := *maxRegionId + 1
	regions := NewStack[*rcRegion[int]](func() *rcRegion[int] {
		return newRcRegion[int](0)
	})

	// Construct regions
	regions.Reserve(nreg)
	for i := 0; i < nreg; i++ {
		regions.Push(newRcRegion[int](i))
	}

	// Find region neighbours and overlapping regions.
	lregs := NewStackArray(func() int { return 0 }, 32)
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			c := chf.cells[x+y*w]

			lregs.Clear()
			i := c.index
			ni := c.index + c.count
			for ; i < ni; i++ {
				s := chf.spans[i]
				ri := srcReg[i]
				if ri == 0 || ri >= nreg {
					continue
				}
				reg := regions.Index(ri)

				reg.spanCount++

				reg.ymin = rcMin(reg.ymin, s.y)
				reg.ymax = rcMax(reg.ymax, s.y)

				// Collect all region layers.
				lregs.Push(ri)

				// Update neighbours
				for dir := 0; dir < 4; dir++ {
					if rcGetCon(s, dir) != RC_NOT_CONNECTED {
						ax := x + rcGetDirOffsetX(dir)
						ay := y + rcGetDirOffsetY(dir)
						ai := chf.cells[ax+ay*w].index + rcGetCon(s, dir)
						rai := srcReg[ai]
						if rai > 0 && rai < nreg && rai != ri {
							addUniqueConnection(reg, rai)
						}

						if rai&RC_BORDER_REG > 0 {
							reg.connectsToBorder = true
						}

					}
				}

			}

			// Update overlapping regions.
			for i = 0; i < lregs.Len()-1; i++ {
				for j := i + 1; j < lregs.Len(); j++ {
					if lregs.Index(i) != lregs.Index(j) {
						ri := regions.Index(lregs.Index(i))
						rj := regions.Index(lregs.Index(j))
						addUniqueFloorRegion(ri, lregs.Index(j))
						addUniqueFloorRegion(rj, lregs.Index(i))
					}
				}
			}

		}
	}

	// Create 2D layers from regions.
	layerId := 1

	for i := 0; i < nreg; i++ {
		regions.Index(i).id = 0
	}

	// Merge montone regions to create non-overlapping areas.
	stack := NewStackArray(func() int {
		return 0
	}, 32)
	for i := 1; i < nreg; i++ {
		root := regions.Index(i)
		// Skip already visited.
		if root.id != 0 {
			continue
		}

		// Start search.
		root.id = layerId

		stack.Clear()
		stack.Push(i)

		for stack.Len() > 0 {
			// Pop front
			reg := regions.Index(stack.Index(0))
			for j := 0; j < stack.Len()-1; j++ {
				stack.SetByIndex(j, stack.Index(j+1))
			}

			stack.Resize(stack.Len() - 1)

			ncons := reg.connections.Len()
			for j := 0; j < ncons; j++ {
				nei := reg.connections.Index(j)
				regn := regions.Index(nei)
				// Skip already visited.
				if regn.id != 0 {
					continue
				}

				// Skip if the neighbour is overlapping root region.
				overlap := false
				for k := 0; k < root.floors.Len(); k++ {
					if root.floors.Index(k) == nei {
						overlap = true
						break
					}
				}
				if overlap {
					continue
				}

				// Deepen
				stack.Push(nei)

				// Mark layer id
				regn.id = layerId
				// Merge current layers to root.
				for k := 0; k < regn.floors.Len(); k++ {
					addUniqueFloorRegion(root, regn.floors.Index(k))
				}

				root.ymin = rcMin(root.ymin, regn.ymin)
				root.ymax = rcMax(root.ymax, regn.ymax)
				root.spanCount += regn.spanCount
				regn.spanCount = 0
				root.connectsToBorder = root.connectsToBorder || regn.connectsToBorder
			}
		}

		layerId++
	}

	// Remove small regions
	for i := 0; i < nreg; i++ {
		if regions.Index(i).spanCount > 0 && regions.Index(i).spanCount < minRegionArea && !regions.Index(i).connectsToBorder {
			reg := regions.Index(i).id
			for j := 0; j < nreg; j++ {
				if regions.Index(j).id == reg {
					regions.Index(j).id = 0
				}

			}

		}
	}

	// Compress region Ids.
	for i := 0; i < nreg; i++ {
		regions.Index(i).remap = false
		if regions.Index(i).id == 0 {
			continue
		} // Skip nil regions.
		if regions.Index(i).id&RC_BORDER_REG > 0 {
			continue
		} // Skip external regions.
		regions.Index(i).remap = true
	}

	regIdGen := 0
	for i := 0; i < nreg; i++ {
		if !regions.Index(i).remap {
			continue
		}

		oldId := regions.Index(i).id
		regIdGen++
		newId := regIdGen
		for j := i; j < nreg; j++ {
			if regions.Index(j).id == oldId {
				regions.Index(j).id = newId
				regions.Index(j).remap = false
			}
		}
	}
	*maxRegionId = regIdGen

	// Remap regions.
	for i := 0; i < chf.spanCount; i++ {
		if (srcReg[i] & RC_BORDER_REG) == 0 {
			srcReg[i] = regions.Index(srcReg[i]).id
		}

	}

	return true
}

func mergeAndFilterRegions(minRegionArea, mergeRegionSize int,
	maxRegionId *int,
	chf *rcCompactHeightfield,
	srcReg []int, overlaps Stack[int]) bool {
	w := chf.width
	h := chf.height

	nreg := *maxRegionId + 1
	regions := NewStack[*rcRegion[int]](func() *rcRegion[int] {
		return newRcRegion[int](0)
	})
	regions.Reserve(nreg)
	// Construct regions
	for i := 0; i < nreg; i++ {
		regions.Push(newRcRegion[int](i))
	}

	// Find edge of a region and find connections around the contour.
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			c := chf.cells[x+y*w]
			i := c.index
			ni := (c.index + c.count)
			for ; i < ni; i++ {
				r := srcReg[i]
				if r == 0 || r >= nreg {
					continue
				}

				reg := regions.Index(r)
				reg.spanCount++

				// Update floors.
				for j := c.index; j < ni; j++ {
					if i == j {
						continue
					}
					floorId := srcReg[j]
					if floorId == 0 || floorId >= nreg {
						continue
					}

					if floorId == r {
						reg.overlap = true
					}
					addUniqueFloorRegion(reg, floorId)
				}

				// Have found contour
				if reg.connections.Len() > 0 {
					continue
				}

				reg.areaType = chf.areas[i]

				// Check if this cell is next to a border.
				ndir := -1
				for dir := 0; dir < 4; dir++ {
					if isSolidEdge(chf, srcReg, x, y, i, dir) {
						ndir = dir
						break
					}
				}

				if ndir != -1 {
					// The cell is at border.
					// Walk around the contour to find all the neighbours.
					regionWalkContour(x, y, i, ndir, chf, srcReg, reg.connections)
				}
			}
		}
	}

	// Remove too small regions.
	stack := NewStackArray[int](func() int { return 0 }, 32)
	trace := NewStackArray[int](func() int {
		return 0
	}, 32)
	for i := 0; i < nreg; i++ {
		reg := regions.Index(i)
		if reg.id == 0 || (reg.id&RC_BORDER_REG > 0) {
			continue
		}

		if reg.spanCount == 0 {
			continue
		}

		if reg.visited {
			continue
		}

		// Count the total size of all the connected regions.
		// Also keep track of the regions connects to a tile border.
		connectsToBorder := false
		spanCount := 0
		stack.Clear()
		trace.Clear()

		reg.visited = true
		stack.Push(i)

		for !stack.Empty() {
			// Pop
			ri := stack.Pop()

			creg := regions.Index(ri)

			spanCount += creg.spanCount
			trace.Push(ri)

			for j := 0; j < creg.connections.Len(); j++ {
				if creg.connections.Index(j)&RC_BORDER_REG > 0 {
					connectsToBorder = true
					continue
				}
				neireg := regions.Index(creg.connections.Index(j))
				if neireg.visited {
					continue
				}

				if neireg.id == 0 || (neireg.id&RC_BORDER_REG) > 0 {
					continue
				}

				// Visit
				stack.Push(neireg.id)
				neireg.visited = true
			}
		}

		// If the accumulated regions size is too small, remove it.
		// Do not remove areas which connect to tile borders
		// as their size cannot be estimated correctly and removing them
		// can potentially remove necessary areas.
		if spanCount < minRegionArea && !connectsToBorder {
			// Kill all visited regions.
			for j := 0; j < trace.Len(); j++ {
				regions.Index(trace.Index(j)).spanCount = 0
				regions.Index(trace.Index(j)).id = 0
			}
		}
	}

	// Merge too small regions to neighbour regions.
	mergeCount := 0
BEGIN:
	mergeCount = 0
	for i := 0; i < nreg; i++ {
		reg := regions.Index(i)
		if reg.id == 0 || (reg.id&RC_BORDER_REG > 0) {
			continue
		}

		if reg.overlap {
			continue
		}

		if reg.spanCount == 0 {
			continue
		}

		// Check to see if the region should be merged.
		if reg.spanCount > mergeRegionSize && isRegionConnectedToBorder(reg) {
			continue
		}

		// Small region with more than 1 connection.
		// Or region which is not connected to a border at all.
		// Find smallest neighbour region that connects to this one.
		smallest := 0xfffffff
		mergeId := reg.id
		for j := 0; j < reg.connections.Len(); j++ {
			if reg.connections.Index(j)&RC_BORDER_REG > 0 {
				continue
			}
			mreg := regions.Index(reg.connections.Index(j))
			if mreg.id == 0 || (mreg.id&RC_BORDER_REG > 0) || mreg.overlap {
				continue
			}
			if mreg.spanCount < smallest &&
				canMergeWithRegion(reg, mreg) &&
				canMergeWithRegion(mreg, reg) {
				smallest = mreg.spanCount
				mergeId = mreg.id
			}
		}
		// Found new id.
		if mergeId != reg.id {
			oldId := reg.id
			target := regions.Index(mergeId)

			// Merge neighbours.
			if mergeRegions(target, reg) {
				// Fixup regions pointing to current region.
				for j := 0; j < nreg; j++ {
					if regions.Index(j).id == 0 || (regions.Index(j).id&RC_BORDER_REG) > 0 {
						continue
					}
					// If another region was already merged into current region
					// change the nid of the previous region too.
					if regions.Index(j).id == oldId {
						regions.Index(j).id = mergeId
					}

					// Replace the current region with the new one if the
					// current regions is neighbour.
					replaceNeighbour(regions.Index(j), oldId, mergeId)
				}
				mergeCount++
			}
		}
	}

	if mergeCount > 0 {
		goto BEGIN
	}

	// Compress region Ids.
	for i := 0; i < nreg; i++ {
		regions.Index(i).remap = false
		if regions.Index(i).id == 0 {
			continue
		} // Skip nil regions.
		if regions.Index(i).id&RC_BORDER_REG > 0 {
			continue
		} // Skip external regions.
		regions.Index(i).remap = true
	}

	regIdGen := 0
	for i := 0; i < nreg; i++ {
		if !regions.Index(i).remap {
			continue
		}

		oldId := regions.Index(i).id
		regIdGen++
		newId := regIdGen
		for j := i; j < nreg; j++ {
			if regions.Index(j).id == oldId {
				regions.Index(j).id = newId
				regions.Index(j).remap = false
			}
		}
	}
	*maxRegionId = regIdGen

	// Remap regions.
	for i := 0; i < chf.spanCount; i++ {
		if (srcReg[i] & RC_BORDER_REG) == 0 {
			srcReg[i] = regions.Index(srcReg[i]).id
		}

	}

	// Return regions that we found to be overlapping.
	for i := 0; i < nreg; i++ {
		if regions.Index(i).overlap {
			overlaps.Push(regions.Index(i).id)
		}
	}

	return true
}

// / @par
// /
// / This is usually the second to the last step in creating a fully built
// / compact heightfield.  This step is required before regions are built
// / using #rcBuildRegions or #rcBuildRegionsMonotone.
// /
// / After this step, the distance data is available via the rcCompactHeightfield::maxDistance
// / and rcCompactHeightfield::dist fields.
// /
// / @see rcCompactHeightfield, rcBuildRegions, rcBuildRegionsMonotone
func rcBuildDistanceField(chf *rcCompactHeightfield) bool {
	if chf.dist != nil {
		chf.dist = nil
		chf.dist = []int{}
	}
	src := make([]int, chf.spanCount)
	dst := make([]int, chf.spanCount)

	maxDist := 0
	calculateDistanceField(chf, src, &maxDist)
	chf.maxDistance = maxDist

	// Blur
	if !reflect.DeepEqual(boxBlur(chf, 1, src, dst), src) {
		dst, src = src, dst
	}
	// Store distance.
	chf.dist = src

	return true
}
func paintRectRegion(minx, maxx, miny, maxy, regId int,
	chf *rcCompactHeightfield, srcReg []int) {
	w := chf.width
	for y := miny; y < maxy; y++ {
		for x := minx; x < maxx; x++ {
			c := chf.cells[x+y*w]
			i := c.index
			ni := c.index + c.count
			for ; i < ni; i++ {
				if chf.areas[i] != RC_NULL_AREA {
					srcReg[i] = regId
				}

			}
		}
	}
}

const RC_NULL_NEI = 0xffff

type rcSweepSpan struct {
	rid int // row id
	id  int // region id
	ns  int // number samples
	nei int // neighbour id
}

// / @par
// /
// / Non-null regions will consist of connected, non-overlapping walkable spans that form a single contour.
// / Contours will form simple polygons.
// /
// / If multiple regions form an area that is smaller than @p minRegionArea, then all spans will be
// / re-assigned to the zero (null) region.
// /
// / Partitioning can result in smaller than necessary regions. @p mergeRegionArea helps
// / reduce unnecessarily small regions.
// /
// / See the #rcConfig documentation for more information on the configuration parameters.
// /
// / The region data will be available via the rcCompactHeightfield::maxRegions
// / and rcCompactSpan::reg fields.
// /
// / @warning The distance field must be created using #rcBuildDistanceField before attempting to build regions.
// /
// / @see rcCompactHeightfield, rcCompactSpan, rcBuildDistanceField, rcBuildRegionsMonotone, rcConfig
func rcBuildRegionsMonotone(chf *rcCompactHeightfield,
	borderSize, minRegionArea, mergeRegionArea int) bool {

	w := chf.width
	h := chf.height
	id := 1
	srcReg := make([]int, chf.spanCount)

	nsweeps := rcMax(chf.width, chf.height)
	sweeps := make([]*rcSweepSpan, nsweeps)

	// Mark border regions.
	if borderSize > 0 {
		// Make sure border will not overflow.
		bw := rcMin(w, borderSize)
		bh := rcMin(h, borderSize)
		// Paint regions
		paintRectRegion(0, bw, 0, h, id|RC_BORDER_REG, chf, srcReg)
		id++
		paintRectRegion(w-bw, w, 0, h, id|RC_BORDER_REG, chf, srcReg)
		id++
		paintRectRegion(0, w, 0, bh, id|RC_BORDER_REG, chf, srcReg)
		id++
		paintRectRegion(0, w, h-bh, h, id|RC_BORDER_REG, chf, srcReg)
		id++
	}

	chf.borderSize = borderSize

	prev := NewStackArray(func() int {
		return 0
	}, 256)

	// Sweep one line at a time.
	for y := borderSize; y < h-borderSize; y++ {
		// Collect spans from this row.
		prev.Resize(id + 1)
		prev.Resize(id, 0)
		rid := 1

		for x := borderSize; x < w-borderSize; x++ {
			c := chf.cells[x+y*w]

			i := c.index
			ni := (c.index + c.count)
			for ; i < ni; i++ {
				s := chf.spans[i]
				if chf.areas[i] == RC_NULL_AREA {
					continue
				}

				// -x
				previd := 0
				if rcGetCon(s, 0) != RC_NOT_CONNECTED {
					ax := x + rcGetDirOffsetX(0)
					ay := y + rcGetDirOffsetY(0)
					ai := chf.cells[ax+ay*w].index + rcGetCon(s, 0)
					if (srcReg[ai]&RC_BORDER_REG) == 0 && chf.areas[i] == chf.areas[ai] {
						previd = srcReg[ai]
					}

				}

				if previd == 0 {
					previd = rid
					rid++
					sweeps[previd].rid = previd
					sweeps[previd].ns = 0
					sweeps[previd].nei = 0
				}

				// -y
				if rcGetCon(s, 3) != RC_NOT_CONNECTED {
					ax := x + rcGetDirOffsetX(3)
					ay := y + rcGetDirOffsetY(3)
					ai := chf.cells[ax+ay*w].index + rcGetCon(s, 3)
					if srcReg[ai] > 0 && (srcReg[ai]&RC_BORDER_REG) == 0 && chf.areas[i] == chf.areas[ai] {
						nr := srcReg[ai]
						if sweeps[previd].nei == 0 || sweeps[previd].nei == nr {
							sweeps[previd].nei = nr
							sweeps[previd].ns++
							prev.SetByIndex(nr, prev.Index(nr)+1)
						} else {
							sweeps[previd].nei = RC_NULL_NEI
						}
					}
				}

				srcReg[i] = previd
			}
		}

		// Create unique ID.
		for i := 1; i < rid; i++ {
			if sweeps[i].nei != RC_NULL_NEI && sweeps[i].nei != 0 && prev.Index(sweeps[i].nei) == sweeps[i].ns {
				sweeps[i].id = sweeps[i].nei
			} else {
				sweeps[i].id = id
				id++
			}
		}

		// Remap IDs
		for x := borderSize; x < w-borderSize; x++ {
			c := chf.cells[x+y*w]
			i := c.index
			ni := (c.index + c.count)
			for ; i < ni; i++ {
				if srcReg[i] > 0 && srcReg[i] < rid {
					srcReg[i] = sweeps[srcReg[i]].id
				}

			}
		}
	}

	// Merge regions and filter out small regions.
	overlaps := NewStack(func() int {
		return 0
	})
	chf.maxRegions = id
	if !mergeAndFilterRegions(minRegionArea, mergeRegionArea, &chf.maxRegions, chf, srcReg, overlaps) {
		return false
	}

	// Monotone partitioning does not generate overlapping regions.

	// Store the result out.
	for i := 0; i < chf.spanCount; i++ {
		chf.spans[i].reg = srcReg[i]
	}

	return true
}

// / @par
// /
// / Non-null regions will consist of connected, non-overlapping walkable spans that form a single contour.
// / Contours will form simple polygons.
// /
// / If multiple regions form an area that is smaller than @p minRegionArea, then all spans will be
// / re-assigned to the zero (null) region.
// /
// / Watershed partitioning can result in smaller than necessary regions, especially in diagonal corridors.
// / @p mergeRegionArea helps reduce unnecessarily small regions.
// /
// / See the #rcConfig documentation for more information on the configuration parameters.
// /
// / The region data will be available via the rcCompactHeightfield::maxRegions
// / and rcCompactSpan::reg fields.
// /
// / @warning The distance field must be created using #rcBuildDistanceField before attempting to build regions.
// /
// / @see rcCompactHeightfield, rcCompactSpan, rcBuildDistanceField, rcBuildRegionsMonotone, rcConfig
func rcBuildRegions(chf *rcCompactHeightfield, borderSize, minRegionArea, mergeRegionArea int) bool {

	w := chf.width
	h := chf.height
	buf := make([]int, chf.spanCount*2)
	LOG_NB_STACKS := 3
	NB_STACKS := 1 << LOG_NB_STACKS
	lvlStacks := make([]Stack[*LevelStackEntry], NB_STACKS)
	for i := 0; i < NB_STACKS; i++ {
		lvlStacks[i].Reserve(256)
	}

	stack := NewStack(func() *LevelStackEntry {
		return newLevelStackEntry(0, 0, 0)
	})
	stack.Reserve(256)

	srcReg := buf
	srcDist := buf[chf.spanCount:]
	regionId := 1
	level := (chf.maxDistance + 1) & ^1

	// TODO: Figure better formula, expandIters defines how much the
	// watershed "overflows" and simplifies the regions. Tying it to
	// agent radius was usually good indication how greedy it could be.
	//	const int expandIters = 4 + walkableRadius * 2;
	expandIters := 8

	if borderSize > 0 {
		// Make sure border will not overflow.
		bw := rcMin(w, borderSize)
		bh := rcMin(h, borderSize)

		// Paint regions
		paintRectRegion(0, bw, 0, h, regionId|RC_BORDER_REG, chf, srcReg)
		regionId++
		paintRectRegion(w-bw, w, 0, h, regionId|RC_BORDER_REG, chf, srcReg)
		regionId++
		paintRectRegion(0, w, 0, bh, regionId|RC_BORDER_REG, chf, srcReg)
		regionId++
		paintRectRegion(0, w, h-bh, h, regionId|RC_BORDER_REG, chf, srcReg)
		regionId++
	}

	chf.borderSize = borderSize

	sId := -1
	for level > 0 {
		level = 0
		if level >= 2 {
			level = level - 2
		}
		sId = (sId + 1) & (NB_STACKS - 1)
		if sId == 0 {
			sortCellsByLevel(level, chf, srcReg, NB_STACKS, lvlStacks, 1)
		} else {
			appendStacks(lvlStacks[sId-1], lvlStacks[sId], srcReg) // copy left overs from last level}

			// Expand current regions until no empty connected cells found.
			expandRegions(expandIters, level, chf, srcReg, srcDist, lvlStacks[sId], false)
			// Mark new regions with IDs.
			for j := 0; j < lvlStacks[sId].Len(); j++ {
				current := lvlStacks[sId].Index(j)
				x := current.x
				y := current.y
				i := current.index
				if i >= 0 && srcReg[i] == 0 {
					if floodRegion(x, y, i, level, regionId, chf, srcReg, srcDist, stack) {
						if regionId == 0xFFFF {
							return false
						}

						regionId++
					}
				}
			}
		}
	}
	// Expand current regions until no empty connected cells found.
	expandRegions(expandIters*8, 0, chf, srcReg, srcDist, stack, true)

	// Merge regions and filter out small regions.
	overlaps := NewStack(func() int {
		return 0
	})
	chf.maxRegions = regionId
	if !mergeAndFilterRegions(minRegionArea, mergeRegionArea, &chf.maxRegions, chf, srcReg, overlaps) {
		return false
	}

	// If overlapping regions were found during merging, split those regions.
	if overlaps.Len() > 0 {
		fmt.Printf("rcBuildRegions: %d overlapping regions.", overlaps.Len())
	}

	// Write the result out.
	for i := 0; i < chf.spanCount; i++ {
		chf.spans[i].reg = srcReg[i]
	}

	return true
}

func rcBuildLayerRegions(chf *rcCompactHeightfield, borderSize, minRegionArea int) bool {

	w := chf.width
	h := chf.height
	id := 1
	srcReg := make([]int, chf.spanCount)

	nsweeps := rcMax(chf.width, chf.height)
	sweeps := make([]*rcSweepSpan, nsweeps)

	// Mark border regions.
	if borderSize > 0 {
		// Make sure border will not overflow.
		bw := rcMin(w, borderSize)
		bh := rcMin(h, borderSize)
		// Paint regions
		paintRectRegion(0, bw, 0, h, id|RC_BORDER_REG, chf, srcReg)
		id++
		paintRectRegion(w-bw, w, 0, h, id|RC_BORDER_REG, chf, srcReg)
		id++
		paintRectRegion(0, w, 0, bh, id|RC_BORDER_REG, chf, srcReg)
		id++
		paintRectRegion(0, w, h-bh, h, id|RC_BORDER_REG, chf, srcReg)
		id++
	}

	chf.borderSize = borderSize

	prev := NewStackArray(func() int {
		return 0
	}, 256)

	// Sweep one line at a time.
	for y := borderSize; y < h-borderSize; y++ {
		// Collect spans from this row.
		prev.Resize(id + 1)
		rid := 1

		for x := borderSize; x < w-borderSize; x++ {
			c := chf.cells[x+y*w]

			i := c.index
			ni := (c.index + c.count)
			for ; i < ni; i++ {
				s := chf.spans[i]
				if chf.areas[i] == RC_NULL_AREA {
					continue
				}

				// -x
				previd := 0
				if rcGetCon(s, 0) != RC_NOT_CONNECTED {
					ax := x + rcGetDirOffsetX(0)
					ay := y + rcGetDirOffsetY(0)
					ai := chf.cells[ax+ay*w].index + rcGetCon(s, 0)
					if (srcReg[ai]&RC_BORDER_REG) == 0 && chf.areas[i] == chf.areas[ai] {
						previd = srcReg[ai]
					}

				}

				if previd == 0 {
					previd = rid
					rid++
					sweeps[previd].rid = previd
					sweeps[previd].ns = 0
					sweeps[previd].nei = 0
				}

				// -y
				if rcGetCon(s, 3) != RC_NOT_CONNECTED {
					ax := x + rcGetDirOffsetX(3)
					ay := y + rcGetDirOffsetY(3)
					ai := chf.cells[ax+ay*w].index + rcGetCon(s, 3)
					if srcReg[ai] > 0 && (srcReg[ai]&RC_BORDER_REG) == 0 && chf.areas[i] == chf.areas[ai] {
						nr := srcReg[ai]
						if sweeps[previd].nei == 0 || sweeps[previd].nei == nr {
							sweeps[previd].nei = nr
							sweeps[previd].ns++
							prev.SetByIndex(nr, prev.Index(nr)+1)
						} else {
							sweeps[previd].nei = RC_NULL_NEI
						}
					}
				}

				srcReg[i] = previd
			}
		}

		// Create unique ID.
		for i := 1; i < rid; i++ {
			if sweeps[i].nei != RC_NULL_NEI && sweeps[i].nei != 0 && prev.Index(sweeps[i].nei) == sweeps[i].ns {
				sweeps[i].id = sweeps[i].nei
			} else {
				sweeps[i].id = id
				id++
			}
		}

		// Remap IDs
		for x := borderSize; x < w-borderSize; x++ {
			c := chf.cells[x+y*w]

			i := c.index
			ni := (int)(c.index + c.count)
			for ; i < ni; i++ {
				if srcReg[i] > 0 && srcReg[i] < rid {
					srcReg[i] = sweeps[srcReg[i]].id
				}

			}
		}
	}

	// Merge monotone regions to layers and remove small regions.
	chf.maxRegions = id
	if !mergeAndFilterLayerRegions(minRegionArea, &chf.maxRegions, chf, srcReg) {
		return false
	}

	// Store the result out.
	for i := 0; i < chf.spanCount; i++ {
		chf.spans[i].reg = srcReg[i]
	}

	return true
}
