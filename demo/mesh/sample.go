package mesh

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/gonavmesh/common/rw"
	"github.com/gorustyt/gonavmesh/detour"
	"io"
	"log/slog"
)

const NAVMESHSET_MAGIC = 'M'<<24 | 'S'<<16 | 'E'<<8 | 'T' //'MSET';
const NAVMESHSET_VERSION = 1

type Sample struct {
}

func NewSample() *Sample {
	return &Sample{}
}

type NavMeshSetHeader struct {
	Magic    int
	Version  int
	NumTiles int
	Params   *detour.NavMeshParams
}

func NewNavMeshSetHeader() *NavMeshSetHeader {
	return &NavMeshSetHeader{
		Params: &detour.NavMeshParams{},
	}
}

func (n *NavMeshSetHeader) Encode(w *rw.ReaderWriter) {
	w.WriteInt32(n.Magic)
	w.WriteInt32(n.Version)
	w.WriteInt32(n.NumTiles)
	n.Params.ToBin(w)
}

func (n *NavMeshSetHeader) Decode(r *rw.ReaderWriter) {
	n.Magic = int(r.ReadInt32())
	n.Version = int(r.ReadInt32())
	n.NumTiles = int(r.ReadInt32())
	n.Params.FromBin(r)
}

type NavMeshTileHeader struct {
	TileRef  detour.DtTileRef
	DataSize int
}

func (n *NavMeshTileHeader) Decode(r *rw.ReaderWriter) {
	n.TileRef = detour.DtTileRef(r.ReadInt32())
	n.DataSize = int(r.ReadInt32())
}
func (n *NavMeshTileHeader) Encode(w *rw.ReaderWriter) {
	w.WriteInt32(int32(n.TileRef))
	w.WriteInt32(n.DataSize)
}

func (s *Sample) saveAll(writer fyne.URIWriteCloser, mesh *detour.DtNavMesh) {
	if mesh == nil {
		return
	}
	w := rw.NewNavMeshDataBinWriter()
	// Store header.
	header := NewNavMeshSetHeader()
	header.Magic = NAVMESHSET_MAGIC
	header.Version = NAVMESHSET_VERSION
	header.NumTiles = 0
	for i := int32(0); i < mesh.GetMaxTiles(); i++ {
		tile := mesh.GetTile(int(i))
		if tile == nil || tile.Header == nil || tile.Data == nil {
			continue
		}
		header.NumTiles++
	}
	header.Params = mesh.GetParams()
	header.Encode(w)

	// Store tiles.
	for i := int32(0); i < mesh.GetMaxTiles(); i++ {
		tile := mesh.GetTile(int(i))
		if tile == nil || tile.Header == nil || tile.Data == nil {
			continue
		}

		var tileHeader NavMeshTileHeader
		tileHeader.TileRef = mesh.GetTileRef(tile)
		dataSize := tile.Data.ToBin(w)
		tileHeader.DataSize = dataSize
		tileHeader.Encode(w)
		_, err := writer.Write(w.GetWriteBytes())
		if err != nil {
			panic(err)
		}
	}
	writer.Close()
}

func (s *Sample) Load(reader fyne.URIReadCloser) {
	s.loadAll(reader)
}

func (s *Sample) loadAll(reader fyne.URIReadCloser) detour.IDtNavMesh {
	data, err := io.ReadAll(reader)
	if err != nil {
		panic(err)
	}
	r := rw.NewNavMeshDataBinReader(data)
	// Read header.
	header := NewNavMeshSetHeader()
	header.Decode(r)
	if header.Magic != NAVMESHSET_MAGIC {
		slog.Error("header.Magic != NAVMESHSET_MAGIC")
		return nil
	}
	if header.Version != NAVMESHSET_VERSION {
		slog.Error("header.Version != NAVMESHSET_VERSION ")
		return nil
	}

	mesh, _ := detour.NewDtNavMeshWithParams(header.Params)

	// Read tiles.
	for i := 0; i < header.NumTiles; i++ {
		var tileHeader NavMeshTileHeader
		tileHeader.Decode(r)

		if tileHeader.TileRef == 0 || tileHeader.DataSize == 0 {
			break
		}

		meshData := &detour.NavMeshData{}
		err = meshData.FromBin(r)
		if err != nil {
			panic(err)
		}
		mesh.AddTile(meshData, detour.DT_TILE_FREE_DATA, tileHeader.TileRef)
	}
	reader.Close()
	return mesh
}
