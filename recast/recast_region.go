package recast

import (
	"fmt"
	"github.com/gorustyt/gonavmesh/common"
	"reflect"
)

const (
	/// Heightfield border flag.
	/// If a heightfield region ID has this bit set, then the region is a border
	/// region and its spans are considered un-walkable.
	/// (Used during the region and contour build process.)
	/// @see RcCompactSpan::reg
	RC_BORDER_REG = 0x8000
)

func calculateDistanceField(chf *RcCompactHeightfield, src []uint16, maxDist *uint16) {
	w := chf.Width
	h := chf.Height

	// Init distance and points.
	for i := int32(0); i < chf.SpanCount; i++ {
		src[i] = 0xffff
	}

	// Mark boundary cells.
	for y := int32(0); y < h; y++ {
		for x := int32(0); x < w; x++ {
			c := chf.Cells[x+y*w]
			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni; i++ {
				s := chf.Spans[i]
				area := chf.Areas[i]

				nc := 0
				for dir := int32(0); dir < 4; dir++ {
					if rcGetCon(s, dir) != RC_NOT_CONNECTED {
						ax := x + common.GetDirOffsetX(dir)
						ay := y + common.GetDirOffsetY(dir)
						ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(s, dir)
						if area == chf.Areas[ai] {
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
	for y := int32(0); y < h; y++ {
		for x := int32(0); x < w; x++ {
			c := chf.Cells[x+y*w]
			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni; i++ {
				s := chf.Spans[i]

				if rcGetCon(s, 0) != RC_NOT_CONNECTED {
					// (-1,0)
					ax := x + common.GetDirOffsetX(0)
					ay := y + common.GetDirOffsetY(0)
					ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(s, 0)
					as := chf.Spans[ai]
					if src[ai]+2 < src[i] {
						src[i] = src[ai] + 2
					}

					// (-1,-1)
					if rcGetCon(as, 3) != RC_NOT_CONNECTED {
						aax := ax + common.GetDirOffsetX(3)
						aay := ay + common.GetDirOffsetY(3)
						aai := int32(chf.Cells[aax+aay*w].Index) + rcGetCon(as, 3)
						if src[aai]+3 < src[i] {
							src[i] = src[aai] + 3
						}

					}
				}
				if rcGetCon(s, 3) != RC_NOT_CONNECTED {
					// (0,-1)
					ax := x + common.GetDirOffsetX(3)
					ay := y + common.GetDirOffsetY(3)
					ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(s, 3)
					as := chf.Spans[ai]
					if src[ai]+2 < src[i] {
						src[i] = src[ai] + 2
					}

					// (1,-1)
					if rcGetCon(as, 2) != RC_NOT_CONNECTED {
						aax := ax + common.GetDirOffsetX(2)
						aay := ay + common.GetDirOffsetY(2)
						aai := int32(chf.Cells[aax+aay*w].Index) + rcGetCon(as, 2)
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
			c := chf.Cells[x+y*w]
			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni; i++ {
				s := chf.Spans[i]

				if rcGetCon(s, 2) != RC_NOT_CONNECTED {
					// (1,0)
					ax := x + common.GetDirOffsetX(2)
					ay := y + common.GetDirOffsetY(2)
					ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(s, 2)
					as := chf.Spans[ai]
					if src[ai]+2 < src[i] {
						src[i] = src[ai] + 2
					}

					// (1,1)
					if rcGetCon(as, 1) != RC_NOT_CONNECTED {
						aax := ax + common.GetDirOffsetX(1)
						aay := ay + common.GetDirOffsetY(1)
						aai := int32(chf.Cells[aax+aay*w].Index) + rcGetCon(as, 1)
						if src[aai]+3 < src[i] {
							src[i] = src[aai] + 3
						}

					}
				}
				if rcGetCon(s, 1) != RC_NOT_CONNECTED {
					// (0,1)
					ax := x + common.GetDirOffsetX(1)
					ay := y + common.GetDirOffsetY(1)
					ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(s, 1)
					as := chf.Spans[ai]
					if src[ai]+2 < src[i] {
						src[i] = src[ai] + 2
					}

					// (-1,1)
					if rcGetCon(as, 0) != RC_NOT_CONNECTED {
						aax := ax + common.GetDirOffsetX(0)
						aay := ay + common.GetDirOffsetY(0)
						aai := int32(chf.Cells[aax+aay*w].Index) + rcGetCon(as, 0)
						if src[aai]+3 < src[i] {
							src[i] = src[aai] + 3
						}

					}
				}
			}
		}
	}

	*maxDist = 0
	for i := int32(0); i < chf.SpanCount; i++ {
		*maxDist = common.Max(src[i], *maxDist)
	}

}

func boxBlur(chf *RcCompactHeightfield, thr int32, src, dst []uint16) []uint16 {
	w := chf.Width
	h := chf.Height

	thr *= 2

	for y := int32(0); y < h; y++ {
		for x := int32(0); x < w; x++ {
			c := chf.Cells[x+y*w]
			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni; i++ {
				s := chf.Spans[i]
				cd := src[i]
				if int32(cd) <= thr {
					dst[i] = cd
					continue
				}

				d := cd
				for dir := int32(0); dir < 4; dir++ {
					if rcGetCon(s, dir) != RC_NOT_CONNECTED {
						ax := x + common.GetDirOffsetX(dir)
						ay := y + common.GetDirOffsetY(dir)
						ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(s, dir)
						d += src[ai]

						as := chf.Spans[ai]
						dir2 := (dir + 1) & 0x3
						if rcGetCon(as, dir2) != RC_NOT_CONNECTED {
							ax2 := ax + common.GetDirOffsetX(dir2)
							ay2 := ay + common.GetDirOffsetY(dir2)
							ai2 := int32(chf.Cells[ax2+ay2*w].Index) + rcGetCon(as, dir2)
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
	x     int32
	y     int32
	index int32
}

func newLevelStackEntry(x int32, y int32, index int32) *LevelStackEntry {
	return &LevelStackEntry{
		x:     x,
		y:     y,
		index: index,
	}
}
func floodRegion(x, y, i int32, level, r uint16,
	chf *RcCompactHeightfield, srcReg, srcDist []uint16, stack Stack[*LevelStackEntry]) bool {
	w := chf.Width

	area := chf.Areas[i]

	// Flood fill mark region.
	stack.Clear()
	stack.Push(newLevelStackEntry(x, y, i))
	srcReg[i] = r
	srcDist[i] = 0

	lev := uint16(0)
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

		cs := chf.Spans[ci]

		// Check if any of the neighbours already have a valid region set.
		ar := uint16(0)
		for dir := int32(0); dir < 4; dir++ {
			// 8 connected
			if rcGetCon(cs, dir) != RC_NOT_CONNECTED {
				ax := cx + common.GetDirOffsetX(dir)
				ay := cy + common.GetDirOffsetY(dir)
				ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(cs, dir)
				if chf.Areas[ai] != area {
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

				as := chf.Spans[ai]

				dir2 := (dir + 1) & 0x3
				if rcGetCon(as, dir2) != RC_NOT_CONNECTED {
					ax2 := ax + common.GetDirOffsetX(dir2)
					ay2 := ay + common.GetDirOffsetY(dir2)
					ai2 := int32(chf.Cells[ax2+ay2*w].Index) + rcGetCon(as, dir2)
					if chf.Areas[ai2] != area {
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
		for dir := int32(0); dir < 4; dir++ {
			if rcGetCon(cs, dir) != RC_NOT_CONNECTED {
				ax := cx + common.GetDirOffsetX(dir)
				ay := cy + common.GetDirOffsetY(dir)
				ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(cs, dir)
				if chf.Areas[ai] != area {
					continue
				}

				if chf.Dist[ai] >= lev && srcReg[ai] == 0 {
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
	index     int32
	region    uint16
	distance2 uint16
}

func newDirtyEntry(index int32, region uint16, distance2 uint16) *DirtyEntry {
	return &DirtyEntry{index: index, region: region, distance2: distance2}
}

func expandRegions(maxIter, level int32,
	chf *RcCompactHeightfield,
	srcReg, srcDist []uint16,
	stack Stack[*LevelStackEntry],
	fillStack bool) {
	w := chf.Width
	h := chf.Height

	if fillStack {
		// Find cells revealed by the raised level.
		stack.Clear()
		for y := int32(0); y < h; y++ {
			for x := int32(0); x < w; x++ {
				c := chf.Cells[x+y*w]
				i := c.Index
				ni := (c.Index + c.Count)
				for ; i < ni; i++ {
					if int32(chf.Dist[i]) >= level && srcReg[i] == 0 && chf.Areas[i] != RC_NULL_AREA {
						stack.Push(newLevelStackEntry(x, y, int32(i)))
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
	iter := int32(0)
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
			d2 := uint16(0xffff)
			area := chf.Areas[i]
			s := chf.Spans[i]
			for dir := int32(0); dir < 4; dir++ {
				if rcGetCon(s, dir) == RC_NOT_CONNECTED {
					continue
				}
				ax := x + common.GetDirOffsetX(dir)
				ay := y + common.GetDirOffsetY(dir)
				ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(s, dir)
				if chf.Areas[ai] != area {
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

func sortCellsByLevel(startLevel uint16,
	chf *RcCompactHeightfield,
	srcReg []uint16,
	nbStacks uint32, stacks []Stack[*LevelStackEntry],
	loglevelsPerStack uint16) { // the levels per stack (2 in our case) as a bit shift

	w := chf.Width
	h := chf.Height
	startLevel = startLevel >> loglevelsPerStack

	for j := uint32(0); j < nbStacks; j++ {
		stacks[j].Clear()
	}

	// put all cells in the level range into the appropriate stacks
	for y := int32(0); y < h; y++ {
		for x := int32(0); x < w; x++ {
			c := chf.Cells[x+y*w]
			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni; i++ {
				if chf.Areas[i] == RC_NULL_AREA || srcReg[i] != 0 {
					continue
				}

				level := chf.Dist[i] >> loglevelsPerStack
				sId := startLevel - level
				if uint32(sId) >= nbStacks {
					continue
				}

				if sId < 0 {
					sId = 0
				}

				stacks[sId].Push(newLevelStackEntry(x, y, int32(i)))
			}
		}
	}
}

func appendStacks(srcStack Stack[*LevelStackEntry],
	dstStack Stack[*LevelStackEntry],
	srcReg []uint16) {
	for j := 0; j < srcStack.Len(); j++ {
		i := srcStack.Index(j).index
		if (i < 0) || (srcReg[i] != 0) {
			continue
		}

		dstStack.Push(srcStack.Index(j))
	}
}

type rcRegion[T any] struct {
	spanCount        int32  // Number of spans belonging to this region
	id               uint16 // ID of the region
	areaType         uint8  // Are type.
	remap            bool
	visited          bool
	overlap          bool
	connectsToBorder bool
	ymin, ymax       uint16
	connections      Stack[T]
	floors           Stack[T]
}

func newRcRegion[T any](i uint16) *rcRegion[T] {
	return &rcRegion[T]{
		id:   i,
		ymin: 0xffff,
	}
}

func removeAdjacentNeighbours(reg *rcRegion[int32]) {
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

func replaceNeighbour(reg *rcRegion[int32], oldId, newId uint16) {
	neiChanged := false
	for i := 0; i < reg.connections.Len(); i++ {
		if reg.connections.Index(i) == int32(oldId) {
			reg.connections.SetByIndex(i, int32(newId))
			neiChanged = true
		}
	}
	for i := 0; i < reg.floors.Len(); i++ {
		if reg.floors.Index(i) == int32(oldId) {
			reg.floors.SetByIndex(i, int32(newId))
		}
	}
	if neiChanged {
		removeAdjacentNeighbours(reg)
	}

}

func canMergeWithRegion(rega *rcRegion[int32], regb *rcRegion[int32]) bool {
	if rega.areaType != regb.areaType {
		return false
	}

	n := 0
	for i := 0; i < rega.connections.Len(); i++ {
		if rega.connections.Index(i) == int32(regb.id) {
			n++
		}

	}
	if n > 1 {
		return false
	}

	for i := 0; i < rega.floors.Len(); i++ {
		if rega.floors.Index(i) == int32(regb.id) {
			return false
		}

	}
	return true
}

func addUniqueFloorRegion(reg *rcRegion[int32], n int32) {
	for i := 0; i < reg.floors.Len(); i++ {
		if reg.floors.Index(i) == n {
			return
		}

	}

	reg.floors.Push(n)
}

func mergeRegions(rega, regb *rcRegion[int32]) bool {
	aid := rega.id
	bid := regb.id

	// Duplicate current neighbourhood.
	acon := NewStack[int32](func() int32 {
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
		if acon.Index(i) == int32(bid) {
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
		if bcon.Index(i) == int32(aid) {
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

func isRegionConnectedToBorder(reg *rcRegion[int32]) bool {
	// Region is connected to border if
	// one of the neighbours is null id.
	for i := 0; i < reg.connections.Len(); i++ {
		if reg.connections.Index(i) == 0 {
			return true
		}

	}
	return false
}

func isSolidEdge(chf *RcCompactHeightfield, srcReg []uint16, x, y, i, dir int32) bool {
	s := chf.Spans[i]
	r := uint16(0)
	if rcGetCon(s, dir) != RC_NOT_CONNECTED {
		ax := x + common.GetDirOffsetX(dir)
		ay := y + common.GetDirOffsetY(dir)
		ai := int32(chf.Cells[ax+ay*chf.Width].Index) + rcGetCon(s, dir)
		r = srcReg[ai]
	}
	if r == srcReg[i] {
		return false
	}

	return true
}

func regionWalkContour(x, y, i, dir int32, chf *RcCompactHeightfield, srcReg []uint16,
	cont Stack[int32]) {
	startDir := dir
	starti := i

	ss := chf.Spans[i]
	curReg := uint16(0)
	if rcGetCon(ss, dir) != RC_NOT_CONNECTED {
		ax := x + common.GetDirOffsetX(dir)
		ay := y + common.GetDirOffsetY(dir)
		ai := int32(chf.Cells[ax+ay*chf.Width].Index) + rcGetCon(ss, dir)
		curReg = srcReg[ai]
	}
	cont.Push(int32(curReg))

	iter := 0
	for {
		iter++
		if iter >= 40000 {
			break
		}
		s := chf.Spans[i]

		if isSolidEdge(chf, srcReg, x, y, i, dir) {
			// Choose the edge corner
			r := uint16(0)
			if rcGetCon(s, dir) != RC_NOT_CONNECTED {
				ax := x + common.GetDirOffsetX(dir)
				ay := y + common.GetDirOffsetY(dir)
				ai := int32(chf.Cells[ax+ay*chf.Width].Index) + rcGetCon(s, dir)
				r = srcReg[ai]
			}
			if r != curReg {
				curReg = r
				cont.Push(int32(curReg))
			}

			dir = (dir + 1) & 0x3 // Rotate CW
		} else {
			ni := int32(-1)
			nx := x + common.GetDirOffsetX(dir)
			ny := y + common.GetDirOffsetY(dir)
			if rcGetCon(s, dir) != RC_NOT_CONNECTED {
				nc := chf.Cells[nx+ny*chf.Width]
				ni = int32(nc.Index) + rcGetCon(s, dir)
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

func addUniqueConnection(reg *rcRegion[int32], n int32) {
	for i := 0; i < reg.connections.Len(); i++ {
		if reg.connections.Index(i) == n {
			return
		}

	}

	reg.connections.Push(n)
}
func mergeAndFilterLayerRegions(minRegionArea int32, maxRegionId *uint16, chf *RcCompactHeightfield,
	srcReg []uint16) bool {
	w := chf.Width
	h := chf.Height

	nreg := *maxRegionId + 1
	regions := NewStack[*rcRegion[int32]](func() *rcRegion[int32] {
		return newRcRegion[int32](0)
	})

	// Construct regions
	regions.Reserve(int(nreg))
	for i := 0; i < int(nreg); i++ {
		regions.Push(newRcRegion[int32](uint16(i)))
	}

	// Find region neighbours and overlapping regions.
	lregs := NewStackArray(func() int { return 0 }, 32)
	for y := int32(0); y < h; y++ {
		for x := int32(0); x < w; x++ {
			c := chf.Cells[x+y*w]

			lregs.Clear()
			i := c.Index
			ni := c.Index + c.Count
			for ; i < ni; i++ {
				s := chf.Spans[i]
				ri := srcReg[i]
				if ri == 0 || ri >= nreg {
					continue
				}
				reg := regions.Index(int(ri))

				reg.spanCount++

				reg.ymin = common.Min(reg.ymin, s.Y)
				reg.ymax = common.Max(reg.ymax, s.Y)

				// Collect all region layers.
				lregs.Push(int(ri))

				// Update neighbours
				for dir := int32(0); dir < 4; dir++ {
					if rcGetCon(s, dir) != RC_NOT_CONNECTED {
						ax := x + common.GetDirOffsetX(dir)
						ay := y + common.GetDirOffsetY(dir)
						ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(s, dir)
						rai := srcReg[ai]
						if rai > 0 && rai < nreg && rai != ri {
							addUniqueConnection(reg, int32(rai))
						}

						if rai&RC_BORDER_REG > 0 {
							reg.connectsToBorder = true
						}

					}
				}

			}

			// Update overlapping regions.
			for i := 0; i < lregs.Len()-1; i++ {
				for j := i + 1; j < lregs.Len(); j++ {
					if lregs.Index(i) != lregs.Index(j) {
						ri := regions.Index(lregs.Index(i))
						rj := regions.Index(lregs.Index(j))
						addUniqueFloorRegion(ri, int32(lregs.Index(j)))
						addUniqueFloorRegion(rj, int32(lregs.Index(i)))
					}
				}
			}

		}
	}

	// Create 2D layers from regions.
	layerId := 1

	for i := 0; i < int(nreg); i++ {
		regions.Index(i).id = 0
	}

	// Merge montone regions to create non-overlapping areas.
	stack := NewStackArray(func() int {
		return 0
	}, 32)
	for i := 1; i < int(nreg); i++ {
		root := regions.Index(i)
		// Skip already visited.
		if root.id != 0 {
			continue
		}

		// Start search.
		root.id = uint16(layerId)

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
				regn := regions.Index(int(nei))
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
				stack.Push(int(nei))

				// Mark layer id
				regn.id = uint16(layerId)
				// Merge current layers to root.
				for k := 0; k < regn.floors.Len(); k++ {
					addUniqueFloorRegion(root, regn.floors.Index(k))
				}

				root.ymin = common.Min(root.ymin, regn.ymin)
				root.ymax = common.Max(root.ymax, regn.ymax)
				root.spanCount += regn.spanCount
				regn.spanCount = 0
				root.connectsToBorder = root.connectsToBorder || regn.connectsToBorder
			}
		}

		layerId++
	}

	// Remove small regions
	for i := 0; i < int(nreg); i++ {
		if regions.Index(i).spanCount > 0 && regions.Index(i).spanCount < minRegionArea && !regions.Index(i).connectsToBorder {
			reg := regions.Index(i).id
			for j := 0; j < int(nreg); j++ {
				if regions.Index(j).id == reg {
					regions.Index(j).id = 0
				}

			}

		}
	}

	// Compress region Ids.
	for i := 0; i < int(nreg); i++ {
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
	for i := 0; i < int(nreg); i++ {
		if !regions.Index(i).remap {
			continue
		}

		oldId := regions.Index(i).id
		regIdGen++
		newId := regIdGen
		for j := i; j < int(nreg); j++ {
			if regions.Index(j).id == oldId {
				regions.Index(j).id = uint16(newId)
				regions.Index(j).remap = false
			}
		}
	}
	*maxRegionId = uint16(regIdGen)

	// Remap regions.
	for i := 0; i < int(chf.SpanCount); i++ {
		if (srcReg[i] & RC_BORDER_REG) == 0 {
			srcReg[i] = regions.Index(int(srcReg[i])).id
		}

	}

	return true
}

func mergeAndFilterRegions(minRegionArea, mergeRegionSize int32,
	maxRegionId *uint16,
	chf *RcCompactHeightfield,
	srcReg []uint16, overlaps Stack[int32]) bool {
	w := chf.Width
	h := chf.Height

	nreg := *maxRegionId + 1
	regions := NewStack[*rcRegion[int32]](func() *rcRegion[int32] {
		return newRcRegion[int32](0)
	})
	regions.Reserve(int(nreg))
	// Construct regions
	for i := uint16(0); i < nreg; i++ {
		regions.Push(newRcRegion[int32](i))
	}

	// Find edge of a region and find connections around the contour.
	for y := int32(0); y < h; y++ {
		for x := int32(0); x < w; x++ {
			c := chf.Cells[x+y*w]
			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni; i++ {
				r := srcReg[i]
				if r == 0 || r >= nreg {
					continue
				}

				reg := regions.Index(int(r))
				reg.spanCount++

				// Update floors.
				for j := c.Index; j < ni; j++ {
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
					addUniqueFloorRegion(reg, int32(floorId))
				}

				// Have found contour
				if reg.connections.Len() > 0 {
					continue
				}

				reg.areaType = chf.Areas[i]

				// Check if this cell is next to a border.
				ndir := int32(-1)
				for dir := int32(0); dir < 4; dir++ {
					if isSolidEdge(chf, srcReg, x, y, int32(i), dir) {
						ndir = dir
						break
					}
				}

				if ndir != -1 {
					// The cell is at border.
					// Walk around the contour to find all the neighbours.
					regionWalkContour(x, y, int32(i), ndir, chf, srcReg, reg.connections)
				}
			}
		}
	}

	// Remove too small regions.
	stack := NewStackArray[int32](func() int32 { return 0 }, 32)
	trace := NewStackArray[int32](func() int32 {
		return 0
	}, 32)
	for i := 0; i < int(nreg); i++ {
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
		spanCount := int32(0)
		stack.Clear()
		trace.Clear()

		reg.visited = true
		stack.Push(int32(i))

		for !stack.Empty() {
			// Pop
			ri := stack.Pop()

			creg := regions.Index(int(ri))

			spanCount += creg.spanCount
			trace.Push(ri)

			for j := 0; j < creg.connections.Len(); j++ {
				if creg.connections.Index(j)&RC_BORDER_REG > 0 {
					connectsToBorder = true
					continue
				}
				neireg := regions.Index(int(creg.connections.Index(j)))
				if neireg.visited {
					continue
				}

				if neireg.id == 0 || (neireg.id&RC_BORDER_REG) > 0 {
					continue
				}

				// Visit
				stack.Push(int32(neireg.id))
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
				regions.Index(int(trace.Index(j))).spanCount = 0
				regions.Index(int(trace.Index(j))).id = 0
			}
		}
	}

	// Merge too small regions to neighbour regions.
	mergeCount := 0
	common.DoWhile(func() (stop bool) {
		for i := 0; i < int(nreg); i++ {
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
			smallest := int32(0xfffffff)
			mergeId := reg.id
			for j := 0; j < reg.connections.Len(); j++ {
				if reg.connections.Index(j)&RC_BORDER_REG > 0 {
					continue
				}
				mreg := regions.Index(int(reg.connections.Index(j)))
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
				target := regions.Index(int(mergeId))

				// Merge neighbours.
				if mergeRegions(target, reg) {
					// Fixup regions pointing to current region.
					for j := 0; j < int(nreg); j++ {
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
		return false
	}, func() bool {
		return mergeCount > 0
	})

	// Compress region Ids.
	for i := 0; i < int(nreg); i++ {
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
	for i := 0; i < int(nreg); i++ {
		if !regions.Index(i).remap {
			continue
		}

		oldId := regions.Index(i).id
		regIdGen++
		newId := regIdGen
		for j := i; j < int(nreg); j++ {
			if regions.Index(j).id == oldId {
				regions.Index(j).id = uint16(newId)
				regions.Index(j).remap = false
			}
		}
	}
	*maxRegionId = uint16(regIdGen)

	// Remap regions.
	for i := 0; i < int(chf.SpanCount); i++ {
		if (srcReg[i] & RC_BORDER_REG) == 0 {
			srcReg[i] = regions.Index(int(srcReg[i])).id
		}

	}

	// Return regions that we found to be overlapping.
	for i := 0; i < int(nreg); i++ {
		if regions.Index(i).overlap {
			overlaps.Push(int32(regions.Index(i).id))
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
// / After this step, the distance data is available via the RcCompactHeightfield::maxDistance
// / and RcCompactHeightfield::dist fields.
// /
// / @see RcCompactHeightfield, rcBuildRegions, rcBuildRegionsMonotone
func RcBuildDistanceField(chf *RcCompactHeightfield) bool {
	if chf.Dist != nil {
		chf.Dist = nil
		chf.Dist = []uint16{}
	}
	src := make([]uint16, chf.SpanCount)
	dst := make([]uint16, chf.SpanCount)

	maxDist := uint16(0)
	calculateDistanceField(chf, src, &maxDist)
	chf.MaxDistance = maxDist

	// Blur
	if !reflect.DeepEqual(boxBlur(chf, 1, src, dst), src) {
		dst, src = src, dst
	}
	// Store distance.
	chf.Dist = src

	return true
}
func paintRectRegion(minx, maxx, miny, maxy int32, regId uint16,
	chf *RcCompactHeightfield, srcReg []uint16) {
	w := chf.Width
	for y := miny; y < maxy; y++ {
		for x := minx; x < maxx; x++ {
			c := chf.Cells[x+y*w]
			i := c.Index
			ni := c.Index + c.Count
			for ; i < ni; i++ {
				if chf.Areas[i] != RC_NULL_AREA {
					srcReg[i] = regId
				}

			}
		}
	}
}

const RC_NULL_NEI = 0xffff

type rcSweepSpan struct {
	rid uint16 // row id
	id  uint16 // region id
	ns  uint16 // number samples
	nei uint16 // neighbour id
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
// / See the #RcConfig documentation for more information on the configuration parameters.
// /
// / The region data will be available via the RcCompactHeightfield::maxRegions
// / and RcCompactSpan::reg fields.
// /
// / @warning The distance field must be created using #rcBuildDistanceField before attempting to build regions.
// /
// / @see RcCompactHeightfield, RcCompactSpan, rcBuildDistanceField, rcBuildRegionsMonotone, RcConfig
func RcBuildRegionsMonotone(chf *RcCompactHeightfield,
	borderSize, minRegionArea, mergeRegionArea int32) bool {

	w := chf.Width
	h := chf.Height
	id := uint16(1)
	srcReg := make([]uint16, chf.SpanCount)

	nsweeps := common.Max(chf.Width, chf.Height)
	sweeps := make([]*rcSweepSpan, nsweeps)

	// Mark border regions.
	if borderSize > 0 {
		// Make sure border will not overflow.
		bw := common.Min(w, borderSize)
		bh := common.Min(h, borderSize)
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

	chf.BorderSize = borderSize

	prev := NewStackArray(func() int {
		return 0
	}, 256)

	// Sweep one line at a time.
	for y := borderSize; y < h-borderSize; y++ {
		// Collect spans from this row.
		prev.Resize(int(id + 1))
		prev.Resize(int(id), 0)
		rid := uint16(1)

		for x := borderSize; x < w-borderSize; x++ {
			c := chf.Cells[x+y*w]

			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni; i++ {
				s := chf.Spans[i]
				if chf.Areas[i] == RC_NULL_AREA {
					continue
				}

				// -x
				previd := uint16(0)
				if rcGetCon(s, 0) != RC_NOT_CONNECTED {
					ax := x + common.GetDirOffsetX(0)
					ay := y + common.GetDirOffsetY(0)
					ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(s, 0)
					if (srcReg[ai]&RC_BORDER_REG) == 0 && chf.Areas[i] == chf.Areas[ai] {
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
					ax := x + common.GetDirOffsetX(3)
					ay := y + common.GetDirOffsetY(3)
					ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(s, 3)
					if srcReg[ai] > 0 && (srcReg[ai]&RC_BORDER_REG) == 0 && chf.Areas[i] == chf.Areas[ai] {
						nr := srcReg[ai]
						if sweeps[previd].nei == 0 || sweeps[previd].nei == nr {
							sweeps[previd].nei = nr
							sweeps[previd].ns++
							prev.SetByIndex(int(nr), prev.Index(int(nr))+1)
						} else {
							sweeps[previd].nei = RC_NULL_NEI
						}
					}
				}

				srcReg[i] = uint16(previd)
			}
		}

		// Create unique ID.
		for i := uint16(1); i < rid; i++ {
			if sweeps[i].nei != RC_NULL_NEI && sweeps[i].nei != 0 && prev.Index(int(sweeps[i].nei)) == int(sweeps[i].ns) {
				sweeps[i].id = sweeps[i].nei
			} else {
				sweeps[i].id = id
				id++
			}
		}

		// Remap IDs
		for x := borderSize; x < w-borderSize; x++ {
			c := chf.Cells[x+y*w]
			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni; i++ {
				if srcReg[i] > 0 && srcReg[i] < rid {
					srcReg[i] = uint16(sweeps[srcReg[i]].id)
				}

			}
		}
	}

	// Merge regions and filter out small regions.
	overlaps := NewStack(func() int32 {
		return 0
	})
	chf.MaxRegions = id
	if !mergeAndFilterRegions(minRegionArea, mergeRegionArea, &chf.MaxRegions, chf, srcReg, overlaps) {
		return false
	}

	// Monotone partitioning does not generate overlapping regions.

	// Store the result out.
	for i := int32(0); i < chf.SpanCount; i++ {
		chf.Spans[i].Reg = srcReg[i]
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
// / See the #RcConfig documentation for more information on the configuration parameters.
// /
// / The region data will be available via the RcCompactHeightfield::maxRegions
// / and RcCompactSpan::reg fields.
// /
// / @warning The distance field must be created using #rcBuildDistanceField before attempting to build regions.
// /
// / @see RcCompactHeightfield, RcCompactSpan, rcBuildDistanceField, rcBuildRegionsMonotone, RcConfig
func RcBuildRegions(chf *RcCompactHeightfield, borderSize, minRegionArea, mergeRegionArea int32) bool {

	w := chf.Width
	h := chf.Height
	buf := make([]uint16, chf.SpanCount*2)
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
	srcDist := buf[chf.SpanCount:]
	regionId := uint16(1)
	level := (chf.MaxDistance + 1) & (^uint16(1))

	// TODO: Figure better formula, expandIters defines how much the
	// watershed "overflows" and simplifies the regions. Tying it to
	// agent radius was usually good indication how greedy it could be.
	//	const int expandIters = 4 + walkableRadius * 2;
	expandIters := int32(8)

	if borderSize > 0 {
		// Make sure border will not overflow.
		bw := common.Min(w, borderSize)
		bh := common.Min(h, borderSize)

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

	chf.BorderSize = borderSize

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
			expandRegions(expandIters, int32(level), chf, srcReg, srcDist, lvlStacks[sId], false)
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
	overlaps := NewStack(func() int32 {
		return 0
	})
	chf.MaxRegions = regionId
	if !mergeAndFilterRegions(minRegionArea, mergeRegionArea, &chf.MaxRegions, chf, srcReg, overlaps) {
		return false
	}

	// If overlapping regions were found during merging, split those regions.
	if overlaps.Len() > 0 {
		fmt.Printf("rcBuildRegions: %d overlapping regions.", overlaps.Len())
	}

	// Write the result out.
	for i := int32(0); i < chf.SpanCount; i++ {
		chf.Spans[i].Reg = srcReg[i]
	}

	return true
}

func RcBuildLayerRegions(chf *RcCompactHeightfield, borderSize, minRegionArea int32) bool {

	w := chf.Width
	h := chf.Height
	id := uint16(1)
	srcReg := make([]uint16, chf.SpanCount)

	nsweeps := common.Max(chf.Width, chf.Height)
	sweeps := make([]*rcSweepSpan, nsweeps)

	// Mark border regions.
	if borderSize > 0 {
		// Make sure border will not overflow.
		bw := common.Min(w, borderSize)
		bh := common.Min(h, borderSize)
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

	chf.BorderSize = borderSize

	prev := NewStackArray(func() int {
		return 0
	}, 256)

	// Sweep one line at a time.
	for y := borderSize; y < h-borderSize; y++ {
		// Collect spans from this row.
		prev.Resize(int(id + 1))
		rid := uint16(1)

		for x := borderSize; x < w-borderSize; x++ {
			c := chf.Cells[x+y*w]

			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni; i++ {
				s := chf.Spans[i]
				if chf.Areas[i] == RC_NULL_AREA {
					continue
				}

				// -x
				previd := uint16(0)
				if rcGetCon(s, 0) != RC_NOT_CONNECTED {
					ax := x + common.GetDirOffsetX(0)
					ay := y + common.GetDirOffsetY(0)
					ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(s, 0)
					if (srcReg[ai]&RC_BORDER_REG) == 0 && chf.Areas[i] == chf.Areas[ai] {
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
					ax := x + common.GetDirOffsetX(3)
					ay := y + common.GetDirOffsetY(3)
					ai := int32(chf.Cells[ax+ay*w].Index) + rcGetCon(s, 3)
					if srcReg[ai] > 0 && (srcReg[ai]&RC_BORDER_REG) == 0 && chf.Areas[i] == chf.Areas[ai] {
						nr := srcReg[ai]
						if sweeps[previd].nei == 0 || sweeps[previd].nei == nr {
							sweeps[previd].nei = nr
							sweeps[previd].ns++
							prev.SetByIndex(int(nr), prev.Index(int(nr))+1)
						} else {
							sweeps[previd].nei = RC_NULL_NEI
						}
					}
				}

				srcReg[i] = previd
			}
		}

		// Create unique ID.
		for i := uint16(1); i < rid; i++ {
			if sweeps[i].nei != RC_NULL_NEI && sweeps[i].nei != 0 && prev.Index(int(sweeps[i].nei)) == int(sweeps[i].ns) {
				sweeps[i].id = sweeps[i].nei
			} else {
				sweeps[i].id = id
				id++
			}
		}

		// Remap IDs
		for x := borderSize; x < w-borderSize; x++ {
			c := chf.Cells[x+y*w]

			i := c.Index
			ni := int32(c.Index + c.Count)
			for ; int32(i) < ni; i++ {
				if srcReg[i] > 0 && srcReg[i] < rid {
					srcReg[i] = sweeps[srcReg[i]].id
				}

			}
		}
	}

	// Merge monotone regions to layers and remove small regions.
	chf.MaxRegions = id
	if !mergeAndFilterLayerRegions(minRegionArea, &chf.MaxRegions, chf, srcReg) {
		return false
	}

	// Store the result out.
	for i := int32(0); i < chf.SpanCount; i++ {
		chf.Spans[i].Reg = srcReg[i]
	}

	return true
}
