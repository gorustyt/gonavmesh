package detour

import (
	"github.com/gorustyt/gonavmesh/common/rw"
	"unsafe"
)

type NavMeshData struct {
	Header      *DtMeshHeader
	NavVerts    []float32
	NavPolys    []*DtPoly
	NavDMeshes  []*DtPolyDetail
	NavDVerts   []float32
	NavBvtree   []*DtBVNode
	NavDTris    []uint8
	OffMeshCons []*DtOffMeshConnection
	Links       []*DtLink // Ignore links; just leave enough space for them. They'll be created on load.// Ignore links; just leave enough space for them. They'll be created on load.
}

func dtAlign4(x int) int { return (x + 3) & ^3 }

func getStructAlignOffset[T any](v T) int {
	old := unsafe.Sizeof(v)
	return getAlignOffset(int(old))
}

func getAlignOffset(old int) int {
	return dtAlign4(old) - old
}
func (d *NavMeshData) ToProto() {

}

func (d *NavMeshData) FromProto() {

}

func (d *NavMeshData) ToBin() (res []byte) {
	w := rw.NewNavMeshDataBinWriter()
	d.Header.ToBin(w)
	w.PadZero(getStructAlignOffset(*d.Header))
	w.WriteFloat32s(d.NavVerts)
	w.PadZero(getAlignOffset(4 * len(d.NavVerts)))
	for _, v := range d.NavPolys {
		v.ToBin(w)
	}
	w.PadZero(getAlignOffset(int(unsafe.Sizeof(DtPoly{})) * len(d.NavPolys)))

	for _, v := range d.Links {
		v.ToBin(w)
	}
	w.PadZero(getAlignOffset(int(unsafe.Sizeof(DtLink{})) * len(d.Links)))

	for _, v := range d.NavDMeshes {
		v.ToBin(w)
	}
	w.PadZero(getAlignOffset(int(unsafe.Sizeof(DtPolyDetail{})) * len(d.NavDMeshes)))

	w.WriteFloat32(d.NavDVerts)
	w.PadZero(getAlignOffset(4 * len(d.NavDVerts)))

	w.WriteInt8s(d.NavDTris)
	w.PadZero(getAlignOffset(1 * len(d.NavDTris)))

	for _, v := range d.NavBvtree {
		v.ToBin(w)
	}
	w.PadZero(getAlignOffset(int(unsafe.Sizeof(DtBVNode{})) * len(d.NavBvtree)))

	for _, v := range d.OffMeshCons {
		v.ToBin(w)
	}
	w.PadZero(getAlignOffset(int(unsafe.Sizeof(DtOffMeshConnection{})) * len(d.OffMeshCons)))

	return w.GetWriteBytes()
}

func (d *NavMeshData) FromBin(data []byte) {
	r := rw.NewNavMeshDataBinReader(data)
	// Patch header pointers.
	d.Header = (&DtMeshHeader{}).FromBin(r)
	r.Skip(getStructAlignOffset(*d.Header))
	d.NavVerts = make([]float32, d.Header.VertCount*3)
	r.ReadFloat32s(d.NavVerts)
	r.Skip(getAlignOffset(4 * len(d.NavVerts)))
	d.NavPolys = make([]*DtPoly, d.Header.PolyCount)
	for i := range d.NavPolys {
		d.NavPolys[i] = (&DtPoly{}).FromBin(r)
	}
	r.Skip(getAlignOffset(int(unsafe.Sizeof(DtPoly{})) * len(d.NavPolys)))
	d.Links = make([]*DtLink, d.Header.MaxLinkCount)
	for i := range d.Links {
		d.Links[i] = (&DtLink{}).FromBin(r)
	}
	r.Skip(getAlignOffset(int(unsafe.Sizeof(DtLink{})) * len(d.Links)))
	d.NavDMeshes = make([]*DtPolyDetail, d.Header.DetailMeshCount)
	for i := range d.NavDMeshes {
		d.NavDMeshes[i] = (&DtPolyDetail{}).FromBin(r)
	}
	r.Skip(getAlignOffset(int(unsafe.Sizeof(DtPolyDetail{})) * len(d.NavDMeshes)))
	d.NavDVerts = make([]float32, 3*d.Header.DetailVertCount)
	r.ReadFloat32s(d.NavDVerts)
	r.Skip(getAlignOffset(4 * len(d.NavDVerts)))
	d.NavDTris = make([]uint8, 4*d.Header.DetailTriCount)
	r.ReadUInt8s(d.NavDTris)
	r.Skip(getAlignOffset(1 * len(d.NavDTris)))
	d.NavBvtree = make([]*DtBVNode, d.Header.BvNodeCount)
	for i := range d.NavBvtree {
		d.NavBvtree[i] = (&DtBVNode{}).FromBin(r)
	}
	r.Skip(getAlignOffset(int(unsafe.Sizeof(DtBVNode{})) * len(d.NavBvtree)))

	d.OffMeshCons = make([]*DtOffMeshConnection, d.Header.OffMeshConCount)
	for i := range d.OffMeshCons {
		d.OffMeshCons[i] = (&DtOffMeshConnection{}).FromBin(r)
	}
	r.Skip(getAlignOffset(int(unsafe.Sizeof(DtOffMeshConnection{})) * len(d.OffMeshCons)))

}
