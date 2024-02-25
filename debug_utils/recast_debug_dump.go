package debug_utils

import (
	"fmt"
	"github.com/gorustyt/gonavmesh/common/rw"
	"github.com/gorustyt/gonavmesh/recast"
	"log/slog"
)

func DuDumpPolyMeshToObj(pmesh *recast.RcPolyMesh, w *rw.ReaderWriter) bool {
	if w == nil {
		slog.Info("duDumpPolyMeshToObj: input IO is null.\n")
		return false
	}

	nvp := pmesh.Nvp
	cs := pmesh.Cs
	ch := pmesh.Ch
	orig := pmesh.Bmin

	w.WriteString("# Recast Navmesh\n")
	w.WriteString("o NavMesh\n")

	w.WriteString("\n")

	for i := int32(0); i < pmesh.Nverts; i++ {
		v := pmesh.Verts[i*3:]
		x := orig[0] + float32(v[0])*cs
		y := orig[1] + float32(v[1]+1)*ch + 0.1
		z := orig[2] + float32(v[2])*cs
		w.WriteString(fmt.Sprintf("v %f %f %f\n", x, y, z))
	}

	w.WriteString("\n")

	for i := int32(0); i < pmesh.Npolys; i++ {
		p := pmesh.Polys[i*nvp*2:]
		for j := int32(2); j < nvp; j++ {
			if p[j] == recast.RC_MESH_NULL_IDX {
				break
			}
			w.WriteString(fmt.Sprintf("f %d %d %d\n", p[0]+1, p[j-1]+1, p[j]+1))
		}
	}

	return true
}

func DuDumpPolyMeshDetailToObj(dmesh recast.RcPolyMeshDetail, w *rw.ReaderWriter) bool {
	if w == nil {
		slog.Info("duDumpPolyMeshDetailToObj: input IO is null.\n")
		return false
	}
	w.WriteString("# Recast Navmesh\n")
	w.WriteString("o NavMesh\n")

	w.WriteString("\n")

	for i := int32(0); i < dmesh.Nverts; i++ {
		v := dmesh.Verts[i*3:]
		w.WriteString(fmt.Sprintf("v %f %f %f\n", v[0], v[1], v[2]))
	}

	w.WriteString("\n")

	for i := int32(0); i < dmesh.Nmeshes; i++ {
		m := dmesh.Meshes[i*4:]
		bverts := m[0]
		btris := m[2]
		ntris := m[3]
		tris := dmesh.Tris[btris*4:]
		for j := uint32(0); j < ntris; j++ {
			w.WriteString(fmt.Sprintf("f %d %d %d\n",
				(int(bverts)+int(tris[j*4+0]))+1,
				(int(bverts)+int(tris[j*4+1]))+1,
				(int(bverts)+int(tris[j*4+2]))+1))
		}
	}

	return true
}

const CSET_MAGIC = ('c' << 24) | ('s' << 16) | ('e' << 8) | 't'

const CSET_VERSION = 2

func DuDumpContourSet(cset *recast.RcContourSet, w *rw.ReaderWriter) bool {
	if w == nil {
		slog.Info("duDumpContourSet: input IO is null.\n")
		return false
	}

	w.WriteInt32(CSET_MAGIC)
	w.WriteInt32(CSET_VERSION)
	w.WriteInt32(cset.Nconts)
	w.WriteFloat32s(cset.Bmin)
	w.WriteFloat32s(cset.Bmax)

	w.WriteFloat32(cset.Cs)
	w.WriteFloat32(cset.Ch)

	w.WriteInt32(cset.Width)
	w.WriteInt32(cset.Height)
	w.WriteInt32(cset.BorderSize)
	for i := int32(0); i < cset.Nconts; i++ {
		cont := cset.Conts[i]
		w.WriteInt32(cont.Nverts)
		w.WriteInt32(cont.Nrverts)

		w.WriteInt16(cont.Reg)
		w.WriteInt8(cont.Area)
		w.WriteInt32s(cont.Verts)
		w.WriteInt32s(cont.Rverts)
	}

	return true
}

func DuReadContourSet(cset recast.RcContourSet, r *rw.ReaderWriter) bool {
	if r == nil {
		slog.Info("duReadContourSet: input IO is null.\n")
		return false
	}
	magic := r.ReadInt32()
	version := r.ReadInt32()
	if magic != CSET_MAGIC {
		slog.Info("duReadContourSet: Bad voodoo.\n")
		return false
	}
	if version != CSET_VERSION {
		slog.Info("duReadContourSet: Bad version.\n")
		return false
	}

	cset.Nconts = r.ReadInt32()

	cset.Conts = make([]*recast.RcContour, cset.Nconts)
	for i := range cset.Conts {
		cset.Conts[i] = &recast.RcContour{}
	}
	r.ReadFloat32s(cset.Bmin)
	r.ReadFloat32s(cset.Bmax)

	cset.Cs = r.ReadFloat32()
	cset.Ch = r.ReadFloat32()
	cset.Width = r.ReadInt32()
	cset.Height = r.ReadInt32()
	cset.BorderSize = r.ReadInt32()
	for i := int32(0); i < cset.Nconts; i++ {
		cont := cset.Conts[i]
		cont.Nverts = r.ReadInt32()
		cont.Nrverts = r.ReadInt32()
		cont.Reg = r.ReadUInt16()
		cont.Area = r.ReadUInt8()
		cont.Verts = make([]int32, 4*cont.Nverts)
		cont.Rverts = make([]int32, 4*cont.Nrverts)
		r.ReadInt32s(cont.Verts)
		r.ReadInt32s(cont.Rverts)
	}

	return true
}

const CHF_MAGIC = ('r' << 24) | ('c' << 16) | ('h' << 8) | 'f'

const CHF_VERSION = 3

func DuDumpCompactHeightfield(chf *recast.RcCompactHeightfield, w *rw.ReaderWriter) bool {
	if w == nil {
		slog.Info("duDumpCompactHeightfield: input IO is null.\n")
		return false
	}
	w.WriteInt32(CHF_MAGIC)
	w.WriteInt32(CHF_VERSION)
	w.WriteInt32(chf.Width)
	w.WriteInt32(chf.Height)
	w.WriteInt32(chf.SpanCount)
	w.WriteInt32(chf.WalkableHeight)
	w.WriteInt32(chf.WalkableClimb)
	w.WriteInt32(chf.BorderSize)
	w.WriteInt16(chf.MaxDistance)
	w.WriteInt16(chf.MaxRegions)
	w.WriteFloat32s(chf.Bmin[:])
	w.WriteFloat32s(chf.Bmax[:])
	w.WriteFloat32(chf.Cs)
	w.WriteFloat32(chf.Ch)
	tmp := 0
	if len(chf.Cells) != 0 {
		tmp |= 1
	}
	if len(chf.Spans) != 0 {
		tmp |= 2
	}
	if len(chf.Dist) != 0 {
		tmp |= 4
	}
	if len(chf.Areas) != 0 {
		tmp |= 8
	}

	w.WriteInt32(tmp)

	if len(chf.Cells) != 0 {
		for _, v := range chf.Cells {
			v.ToBin(w)
		}
	}

	if len(chf.Spans) != 0 {
		for _, v := range chf.Spans {
			v.ToBin(w)
		}
	}
	if len(chf.Dist) != 0 {
		w.WriteInt16s(chf.Dist)
	}
	if len(chf.Areas) != 0 {
		w.WriteInt8s(chf.Areas)
	}
	return true
}

func DuReadCompactHeightfield(chf *recast.RcCompactHeightfield, r *rw.ReaderWriter) bool {
	if r == nil {
		slog.Info("duReadCompactHeightfield: input IO is null.\n")
		return false
	}

	magic := r.ReadInt32()
	version := r.ReadInt32()
	if magic != CHF_MAGIC {
		slog.Info("duReadCompactHeightfield: Bad voodoo.\n")
		return false
	}
	if version != CHF_VERSION {
		slog.Info("duReadCompactHeightfield: Bad version.\n")
		return false
	}
	chf.Width = r.ReadInt32()
	chf.Height = r.ReadInt32()
	chf.SpanCount = r.ReadInt32()
	chf.WalkableHeight = r.ReadInt32()
	chf.WalkableClimb = r.ReadInt32()
	chf.BorderSize = r.ReadInt32()

	chf.MaxDistance = r.ReadUInt16()
	chf.MaxRegions = r.ReadUInt16()
	r.ReadFloat32s(chf.Bmin[:])
	r.ReadFloat32s(chf.Bmax[:])
	chf.Cs = r.ReadFloat32()
	chf.Ch = r.ReadFloat32()
	tmp := r.ReadInt32()
	if tmp&1 != 0 {
		chf.Cells = make([]*recast.RcCompactCell, chf.Width*chf.Height)
		for i := range chf.Cells {
			chf.Cells[i] = &recast.RcCompactCell{}
			chf.Cells[i].FromBin(r)
		}
	}
	if tmp&2 != 0 {
		chf.Spans = make([]*recast.RcCompactSpan, chf.SpanCount)
		for i := range chf.Spans {
			chf.Spans[i] = &recast.RcCompactSpan{}
			chf.Spans[i].FromBin(r)
		}
	}
	if tmp&4 != 0 {
		chf.Dist = make([]uint16, chf.SpanCount)
		r.ReadUInt16s(chf.Dist)
	}
	if tmp&8 != 0 {
		chf.Areas = make([]uint8, chf.SpanCount)
		r.ReadUInt8s(chf.Areas)
	}

	return true
}

func logLine(name string, pc float32) {
	t := float32(0.)
	slog.Info("%s:\t%.2fms\t(%.1f%%)", name, t/1000.0, t*pc)
}

func DuLogBuildTimes(totalTimeUsec int) {
	pc := 100.0 / float32(totalTimeUsec)

	slog.Info("Build Times")
	logLine("- Rasterize", pc)
	logLine("- Build Compact", pc)
	logLine("- Filter Border", pc)
	logLine("- Filter Walkable", pc)
	logLine("- Erode Area", pc)
	logLine("- Median Area", pc)
	logLine("- Mark Box Area", pc)
	logLine("- Mark Convex Area", pc)
	logLine("- Mark Cylinder Area", pc)
	logLine("- Build Distance Field", pc)
	logLine("    - Distance", pc)
	logLine("    - Blur", pc)
	logLine("- Build Regions", pc)
	logLine("    - Watershed", pc)
	logLine("      - Expand", pc)
	logLine("      - Find Basins", pc)
	logLine("    - Filter", pc)
	logLine("- Build Layers", pc)
	logLine("- Build Contours", pc)
	logLine("    - Trace", pc)
	logLine("    - Simplify", pc)
	logLine("- Build Polymesh", pc)
	logLine("- Build Polymesh Detail", pc)
	logLine("- Merge Polymeshes", pc)
	logLine("- Merge Polymesh Details", pc)
	slog.Info("=== TOTAL:\t%.2fms", totalTimeUsec/1000.0)
}
