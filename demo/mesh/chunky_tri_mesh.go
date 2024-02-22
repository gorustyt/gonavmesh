package mesh

import (
	"math"
	"sort"
)

type rcChunkyTriMeshNode struct {
	bmin [2]float32
	bmax [2]float32
	i    int
	n    int
}

type rcChunkyTriMesh struct {
	nodes           []*rcChunkyTriMeshNode
	nnodes          int
	tris            []int
	ntris           int
	maxTrisPerChunk int
}

func newRcChunkyTriMesh() *rcChunkyTriMesh {
	return &rcChunkyTriMesh{}
}

type BoundsItem struct {
	bmin [2]float32
	bmax [2]float32
	i    int
}

func compareItemX(a, b *BoundsItem) int {

	if a.bmin[0] < b.bmin[0] {
		return -1
	}

	if a.bmin[0] > b.bmin[0] {
		return 1
	}

	return 0
}

func compareItemY(a, b *BoundsItem) int {

	if a.bmin[1] < b.bmin[1] {
		return -1
	}

	if a.bmin[1] > b.bmin[1] {
		return 1
	}

	return 0
}

func calcExtends(items []*BoundsItem, nitems int, imin, imax int, bmin, bmax []float32) {
	bmin[0] = items[imin].bmin[0]
	bmin[1] = items[imin].bmin[1]

	bmax[0] = items[imin].bmax[0]
	bmax[1] = items[imin].bmax[1]

	for i := imin + 1; i < imax; i++ {
		it := items[i]
		if it.bmin[0] < bmin[0] {
			bmin[0] = it.bmin[0]
		}
		if it.bmin[1] < bmin[1] {
			bmin[1] = it.bmin[1]
		}

		if it.bmax[0] > bmax[0] {
			bmax[0] = it.bmax[0]
		}
		if it.bmax[1] > bmax[1] {
			bmax[1] = it.bmax[1]
		}
	}
}

func subdivide(items []*BoundsItem, nitems, imin, imax, trisPerChunk int,
	curNode *int, nodes []*rcChunkyTriMeshNode, maxNodes int,
	curTri *int, outTris []int, inTris []int) {
	inum := imax - imin
	icur := curNode

	if *curNode >= maxNodes {
		return
	}

	node := nodes[*curNode]
	*curNode++
	if inum <= trisPerChunk {
		// Leaf
		calcExtends(items, nitems, imin, imax, node.bmin[:], node.bmax[:])

		// Copy triangles.
		node.i = *curTri
		node.n = inum

		for i := imin; i < imax; i++ {
			src := inTris[items[i].i*3:]
			dst := outTris[*curTri*3:]
			*curTri++
			dst[0] = src[0]
			dst[1] = src[1]
			dst[2] = src[2]
		}
	} else {
		// Split
		calcExtends(items, nitems, imin, imax, node.bmin[:], node.bmax[:])

		axis := longestAxis(node.bmax[0]-node.bmin[0],
			node.bmax[1]-node.bmin[1])

		if axis == 0 {
			// Sort along x-axis
			sort.Slice(items[imin:imin+inum], func(i, j int) bool {
				return compareItemX(items[imin+i], items[imin+j]) > 0
			})

		} else if axis == 1 {
			// Sort along y-axis
			sort.Slice(items[imin:imin+inum], func(i, j int) bool {
				return compareItemY(items[imin+i], items[imin+j]) > 0
			})
		}

		isplit := imin + inum/2

		// Left
		subdivide(items, nitems, imin, isplit, trisPerChunk, curNode, nodes, maxNodes, curTri, outTris, inTris)
		// Right
		subdivide(items, nitems, isplit, imax, trisPerChunk, curNode, nodes, maxNodes, curTri, outTris, inTris)

		iescape := *curNode - *icur
		// Negative index means escape.
		node.i = -iescape
	}
}
func longestAxis(x, y float32) int {
	if y > x {
		return 1
	}
	return 0
}

func rcCreateChunkyTriMesh(verts []float32, tris []int, ntris int,
	trisPerChunk int, cm *rcChunkyTriMesh) bool {
	nchunks := (ntris + trisPerChunk - 1) / trisPerChunk

	cm.nodes = make([]*rcChunkyTriMeshNode, nchunks*4)
	cm.tris = make([]int, ntris*3)
	cm.ntris = ntris

	// Build tree
	items := make([]*BoundsItem, ntris)
	for i := 0; i < ntris; i++ {
		t := tris[i*3:]
		it := items[i]
		it.i = i
		// Calc triangle XZ bounds.
		it.bmax[0] = verts[t[0]*3+0]
		it.bmin[0] = it.bmax[0]
		it.bmax[1] = verts[t[0]*3+2]
		it.bmin[1] = it.bmax[1]
		for j := 1; j < 3; j++ {
			v := verts[t[j]*3:]
			if v[0] < it.bmin[0] {
				it.bmin[0] = v[0]
			}
			if v[2] < it.bmin[1] {
				it.bmin[1] = v[2]
			}

			if v[0] > it.bmax[0] {
				it.bmax[0] = v[0]
			}
			if v[2] > it.bmax[1] {
				it.bmax[1] = v[2]
			}
		}
	}

	curTri := 0
	curNode := 0
	subdivide(items, ntris, 0, ntris, trisPerChunk, &curNode, cm.nodes, nchunks*4, &curTri, cm.tris, tris)

	cm.nnodes = curNode

	// Calc max tris per node.
	cm.maxTrisPerChunk = 0
	for i := 0; i < cm.nnodes; i++ {
		node := cm.nodes[i]
		isLeaf := node.i >= 0
		if !isLeaf {
			continue
		}
		if node.n > cm.maxTrisPerChunk {
			cm.maxTrisPerChunk = node.n
		}

	}

	return true
}

func checkOverlapRect(amin, amax,
	bmin, bmax [2]float32) bool {
	overlap := true
	if amin[0] > bmax[0] || amax[0] < bmin[0] {
		overlap = false
	}
	if amin[1] > bmax[1] || amax[1] < bmin[1] {
		overlap = false
	}
	return overlap
}

func rcGetChunksOverlappingRect(cm *rcChunkyTriMesh,
	bmin, bmax [2]float32,
	ids []int, maxIds int) int {
	// Traverse tree
	i := 0
	n := 0
	for i < cm.nnodes {
		node := cm.nodes[i]
		overlap := checkOverlapRect(bmin, bmax, node.bmin, node.bmax)
		isLeafNode := node.i >= 0

		if isLeafNode && overlap {
			if n < maxIds {
				ids[n] = i
				n++
			}
		}

		if overlap || isLeafNode {
			i++
		} else {
			escapeIndex := -node.i
			i += escapeIndex
		}
	}

	return n
}

func checkOverlapSegment(p, q, bmin, bmax [2]float32) bool {
	EPSILON := 1e-6

	var tmin = 0.0
	var tmax = 1.0
	var d [2]float32
	d[0] = q[0] - p[0]
	d[1] = q[1] - p[1]

	for i := 0; i < 2; i++ {
		if math.Abs(d[i]) < EPSILON {
			// Ray is parallel to slab. No hit if origin not within slab
			if p[i] < bmin[i] || p[i] > bmax[i] {
				return false
			}

		} else {
			// Compute intersection t value of ray with near and far plane of slab
			ood := 1.0 / d[i]
			t1 := (bmin[i] - p[i]) * ood
			t2 := (bmax[i] - p[i]) * ood
			if t1 > t2 {
				tmp := t1
				t1 = t2
				t2 = tmp
			}
			if t1 > tmin {
				tmin = t1
			}
			if t2 < tmax {
				tmax = t2
			}
			if tmin > tmax {
				return false
			}
		}
	}
	return true
}

func rcGetChunksOverlappingSegment(cm *rcChunkyTriMesh,
	p, q [2]float32,
	ids []int, maxIds int) int {
	// Traverse tree
	i := 0
	n := 0
	for i < cm.nnodes {
		node := cm.nodes[i]
		overlap := checkOverlapSegment(p, q, node.bmin, node.bmax)
		isLeafNode := node.i >= 0

		if isLeafNode && overlap {
			if n < maxIds {
				ids[n] = i
				n++
			}
		}

		if overlap || isLeafNode {
			i++
		} else {
			escapeIndex := -node.i
			i += escapeIndex
		}
	}

	return n
}
