package recast

import "math"

type name struct {
}

func prev(i, n int) int {
	if i-1 >= 0 {
		return i - 1
	}
	return n - 1
}

func next(i, n int) int {
	if i+1 < n {
		return i + 1
	}
	return 0
}

// 多边形三角化
func triangulate(verts []*Point, indices []int) [][]int {
	n := len(verts)
	var tris [][]int //三角形
	for i := 0; i < n; i++ {
		i1 := next(i, n)
		i2 := next(i1, n)
		if diagonal(i, i2, n, verts, indices) { //是一个凸点
			indices[i1] |= 0x80000000
		}
	}
	for n > 3 {
		minLen := -1
		mini := -1
		for i := 0; i < n; i++ { //找出长度最小的两个凸点
			i1 := next(i, n)
			if indices[i1]&0x80000000 > 0 {
				p0 := verts[(indices[i]&0x0fffffff)*4]
				p2 := verts[(indices[next(i1, n)]&0x0fffffff)*4]

				dx := p2.X - p0.X
				dy := p2.Y - p0.Y
				length := dx*dx + dy*dy
				if minLen < 0 || length < minLen {
					minLen = length
					mini = i
				}
			}
		}
		if mini == -1 {
			// We might get here because the contour has overlapping segments, like this:
			//
			//  A o-o=====o---o B
			//   /  |C   D|    \.
			//  o   o     o     o
			//  :   :     :     :
			// We'll try to recover by loosing up the inCone test a bit so that a diagonal
			// like A-B or C-D can be found and we can continue.
			minLen = -1
			mini = -1
			for i := 0; i < n; i++ {
				i1 := next(i, n)
				i2 := next(i1, n)
				if diagonalLoose(i, i2, n, verts, indices) {
					p0 := verts[(indices[i]&0x0fffffff)*4]
					p2 := verts[(indices[next(i2, n)]&0x0fffffff)*4]
					dx := p2.X - p0.X
					dy := p2.Y - p0.Y
					length := dx*dx + dy*dy

					if minLen < 0 || length < minLen {
						minLen = length
						mini = i
					}
				}
			}
			if mini == -1 {
				// The contour is messed up. This sometimes happens
				// if the contour simplification is too aggressive.
				return tris
			}
		}

		i := mini
		i1 := next(i, n)
		i2 := next(i1, n)
		var tri []int
		tri = append(tri, indices[i]&0x0fffffff)
		tri = append(tri, indices[i1]&0x0fffffff)
		tri = append(tri, indices[i2]&0x0fffffff)
		tris = append(tris, tri)
		// Removes P[i1] by copying P[i+1]...P[n-1] left one index.
		n--
		for k := i1; k < n; k++ {
			indices[k] = indices[k+1]
		}
		if i1 >= n {
			i1 = 0
		}
		i = prev(i1, n)
		// Update diagonal flags.
		if diagonal(prev(i, n), i1, n, verts, indices) {
			indices[i] |= 0x80000000
		} else {
			indices[i] &= 0x0fffffff
		}

		if diagonal(i, next(i1, n), n, verts, indices) {
			indices[i1] |= 0x80000000
		} else {
			indices[i1] &= 0x0fffffff
		}

	}

	//剩下的三个点
	var tri []int
	tri = append(tri, indices[0]&0x0fffffff)
	tri = append(tri, indices[1]&0x0fffffff)
	tri = append(tri, indices[2]&0x0fffffff)
	tris = append(tris, tri)
	return tris
}

type rcEdge struct {
	vert     [2]*Point //边的两个点
	polyEdge [2]int    //邻接的两个多边形的边的索引
	poly     [2]int    //邻接的两个多边形的索引
}

func buildMeshAdjacency(polys [][]*Point) map[*Point][2]int {
	var (
		edges     []*rcEdge
		firstEdge = map[*Point]int{}
		nextEdge  = map[int]int{}
		adjacency = map[*Point][2]int{}
	)
	for i := 0; i < len(polys); i++ {
		v := polys[i]
		for j := 0; j < len(v); j++ {
			if v[j] == nil {
				break
			}
			v0 := v[j]
			v1 := v[j+1]
			if v1 == nil {
				v1 = v[0]
			}
			edges = append(edges, &rcEdge{
				vert:     [2]*Point{v0, v1},
				polyEdge: [2]int{j, 0},
				poly:     [2]int{i, i},
			})
			nextEdge[len(edges)-1] = firstEdge[v0] //a<-b<-c<-d
			firstEdge[v0] = len(edges) - 1         //v0的一条边
		}
	}
	for i := 0; i < len(polys); i++ {
		v := polys[i]
		for j := 0; j < len(v); j++ {
			if v[j] == nil {
				break
			}
			v0 := v[j]
			v1 := v[j+1]
			if v1 == nil {
				v1 = v[0]
			}
			for e, ok := firstEdge[v1]; ok; e = nextEdge[e] {
				edge := edges[e]
				if edge.vert[1] == v0 && edge.poly[0] == edge.poly[1] {
					edge.poly[1] = i
					edge.polyEdge[1] = j
					break
				}
			}
		}
	}

	for i := 0; i < len(edges); i++ {
		e := edges[i]
		if e.poly[0] != e.poly[1] {

		}
	}
	return adjacency
}

// 构建多边形
func rcBuildPolyMesh(cset *rcContourSet) *rcPolyMesh {
	var (
		mesh rcPolyMesh
		vflags=map[int]struct{}{} //应该移除的顶点,key==>index
	)
	mesh.bmin = cset.bmin
	mesh.bmax = cset.bmax
	mesh.cs = cset.cs
	mesh.ch = cset.ch
	mesh.borderSize = cset.borderSize
	mesh.maxEdgeError = cset.maxError
	maxVertices := 0
	maxTris := 0
	maxVertsPerCont := 0
	for i := 0; i < len(cset.conts); i++ {
		nverts := len(cset.conts[i].verts) //顶点个数
		// Skip null contours.
		if nverts < 3 { //跳过空的顶点,多边形需要三个点
			continue
		}
		maxVertices += nverts
		maxTris += nverts - 2
		maxVertsPerCont = int(math.Max(float64(maxVertsPerCont), float64(nverts)))
	}
	if maxVertices >= 0xfffe { //顶点数量太多
		return nil
	}

	for i := 0; i < len(cset.conts); i++ {
		cont := cset.conts[i]
		nverts := len(cont.verts)
		// Skip null contours.
		if nverts < 3 {
			continue
		}
		var (
			indices = make([]int, nverts) //顶点索引
		)
		// Triangulate contour
		for j := 0; j < nverts; j++ {
			indices[j] = j
		}

		tris := triangulate(cont.verts, indices)
		ntris := len(tris)
		if ntris <= 0 {
			panic("not expect ")
		}
		for j := 0; j <  nverts; j++{
		v := cont.verts[j*4]
		indices[j] = addVertex((unsigned short)v[0], (unsigned short)v[1], (unsigned short)v[2],
		mesh.verts, firstVert, nextVert, mesh.nverts);
		if v[3] & RC_BORDER_VERTEX {
		// This vertex should be removed.
		vflags[indices[j]] = 1;
		}
		}
	}
	return &mesh
}
