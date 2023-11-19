package recast

import "math"

type Item struct {
	id   int
	x, y int
	next int
}

type dtProximityGrid struct {
	m_cellSize    float64
	m_invCellSize float64
	m_pool        []*Item
	m_poolHead    int
	m_poolSize    int

	m_buckets     []int
	m_bucketsSize int

	m_bounds [4]int
}

func hashPos2(x, y, n int) int {
	return ((x * 73856093) ^ (y * 19349663)) & (n - 1)
}

func newDtProximityGrid(poolSize int, cellSize float64) (d *dtProximityGrid) {
	d = &dtProximityGrid{}
	if poolSize <= 0 {
		panic("")
	}
	if cellSize <= 0.0 {
		panic("")
	}
	d.m_cellSize = cellSize
	d.m_invCellSize = 1.0 / d.m_cellSize

	// Allocate hashs buckets
	d.m_bucketsSize = dtNextPow2(poolSize)
	d.m_buckets = make([]int, d.m_bucketsSize)
	// Allocate pool of items.
	d.m_poolSize = poolSize
	d.m_poolHead = 0
	d.m_pool = make([]*Item, d.m_poolSize)

	d.Clear()

	return d
}

func (d *dtProximityGrid) Clear() {
	if d.m_buckets == nil {
		d.m_buckets = make([]int, d.m_bucketsSize)
	}
	for i := range d.m_buckets {
		d.m_buckets[i] = 0xff
	}
	d.m_poolHead = 0
	d.m_bounds[0] = 0xffff
	d.m_bounds[1] = 0xffff
	d.m_bounds[2] = -0xffff
	d.m_bounds[3] = -0xffff
}

func (d *dtProximityGrid) addItem(id int,
	minx, miny,
	maxx, maxy float64) {
	iminx := int(math.Floor(minx * d.m_invCellSize))
	iminy := int(math.Floor(miny * d.m_invCellSize))
	imaxx := int(math.Floor(maxx * d.m_invCellSize))
	imaxy := int(math.Floor(maxy * d.m_invCellSize))

	d.m_bounds[0] = dtMin(d.m_bounds[0], iminx)
	d.m_bounds[1] = dtMin(d.m_bounds[1], iminy)
	d.m_bounds[2] = dtMax(d.m_bounds[2], imaxx)
	d.m_bounds[3] = dtMax(d.m_bounds[3], imaxy)

	for y := iminy; y <= imaxy; y++ {
		for x := iminx; x <= imaxx; x++ {
			if d.m_poolHead < d.m_poolSize {
				h := hashPos2(x, y, d.m_bucketsSize)
				idx := d.m_poolHead
				d.m_poolHead++
				item := d.m_pool[idx]
				item.x = x
				item.y = y
				item.id = id
				item.next = d.m_buckets[h]
				d.m_buckets[h] = idx
			}
		}
	}
}

func (d *dtProximityGrid) queryItems(minx, miny,
	maxx, maxy float64, ids []int, maxIds int) int {
	iminx := int(math.Floor(minx * d.m_invCellSize))
	iminy := int(math.Floor(miny * d.m_invCellSize))
	imaxx := int(math.Floor(maxx * d.m_invCellSize))
	imaxy := int(math.Floor(maxy * d.m_invCellSize))

	n := 0

	for y := iminy; y <= imaxy; y++ {
		for x := iminx; x <= imaxx; x++ {
			h := hashPos2(x, y, d.m_bucketsSize)
			idx := d.m_buckets[h]
			for idx != 0xffff {
				item := d.m_pool[idx]
				if item.x == x && item.y == y {
					// Check if the id exists already.
					end := n
					i := 0
					for i != end && i != item.id {
						i++
					}

					// Item not found, add it.
					if i == end {
						if n >= maxIds {
							return n
						}

						ids[n] = item.id
						n++
					}
				}
				idx = item.next
			}
		}
	}

	return n
}

func (d *dtProximityGrid) getItemCountAt(x, y int) int {
	n := 0

	h := hashPos2(x, y, d.m_bucketsSize)
	idx := d.m_buckets[h]
	for idx != 0xffff {
		item := d.m_pool[idx]
		if item.x == x && item.y == y {
			n++
		}

		idx = item.next
	}

	return n
}
