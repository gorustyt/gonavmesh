package recast

import (
	"github.com/gorustyt/gonavmesh/common/rw"
	"github.com/gorustyt/gonavmesh/detour"
	"unsafe"
)

type NavMeshData struct {
	Header      *detour.DtMeshHeader
	NavVerts    []float32
	NavPolys    []*detour.DtPoly
	NavDMeshes  []*detour.DtPolyDetail
	NavDVerts   []float32
	NavBvtree   []*detour.DtBVNode
	NavDTris    []uint8
	OffMeshCons []*detour.DtOffMeshConnection
	Links       []*detour.DtLink // Ignore links; just leave enough space for them. They'll be created on load.// Ignore links; just leave enough space for them. They'll be created on load.
}

func dtAlign4(x int) int { return (x + 3) & ^3 }

func (d *NavMeshData) ToProto() {

}

func (d *NavMeshData) FromProto() {

}

func (d *NavMeshData) ToBin() (res []byte) {
	w := rw.NewNavMeshDataBinWriter()
	d.Header.ToBin(w)
	// Patch header pointers.
	headerSize = dtAlign4(int(unsafe.Sizeof(dtMeshHeader)));
	vertsSize = dtAlign4(sizeof(float)*3*header->vertCount);
	polysSize = dtAlign4(sizeof(dtPoly)*header->polyCount);
	linksSize = dtAlign4(sizeof(dtLink)*(header->maxLinkCount));
	detailMeshesSize = dtAlign4(sizeof(dtPolyDetail)*header->detailMeshCount);
	detailVertsSize = dtAlign4(sizeof(float)*3*header->detailVertCount);
	detailTrisSize = dtAlign4(sizeof(unsigned char)*4*header->detailTriCount);
	bvtreeSize = dtAlign4(sizeof(dtBVNode)*header->bvNodeCount);
	offMeshLinksSize = dtAlign4(sizeof(dtOffMeshConnection)*header->offMeshConCount);
	return w.GetWriteBytes()
}

func (d *NavMeshData) FromBin(data []byte) {
	r := rw.NewNavMeshDataBinReader(data)
	d.Header.FromBin(r)
}
