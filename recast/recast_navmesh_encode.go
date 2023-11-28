package recast

type NavMeshData struct {
	Header      *DtMeshHeader
	NavVerts    []float64
	NavPolys    []*DtPoly
	navDMeshes  []*DtPolyDetail
	NavDVerts   []float64
	NavBvtree   []*DtBVNode
	NavDTris    []int
	OffMeshCons []*DtOffMeshConnection
	Links       []*DtLink // Ignore links; just leave enough space for them. They'll be created on load.// Ignore links; just leave enough space for them. They'll be created on load.
}
