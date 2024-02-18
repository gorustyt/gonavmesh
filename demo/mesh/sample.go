package mesh

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/gonavmesh/detour"
	"github.com/gorustyt/gonavmesh/recast"
	"log/slog"
)

const NAVMESHSET_MAGIC = 'M'<<24 | 'S'<<16 | 'E'<<8 | 'T' //'MSET';
const NAVMESHSET_VERSION = 1

type Sample struct {
}
type NavMeshSetHeader struct {
	Magic    int
	Version  int
	NumTiles int
	Params   detour.NavMeshParams
}

func (n *NavMeshSetHeader) Encode(writer fyne.URIWriteCloser) {

}

func (n *NavMeshSetHeader) Decode(reader fyne.URIReadCloser) {

}

type NavMeshTileHeader struct {
	TileRef  detour.DtTileRef
	DataSize int
}

func (n *NavMeshTileHeader) Decode(reader fyne.URIReadCloser) {

}
func (n *NavMeshTileHeader) Encode(writer fyne.URIWriteCloser) {

}

func (s *Sample) saveAll(writer fyne.URIWriteCloser, mesh *detour.DtNavMesh) {
	if mesh == nil {
		return
	}
	// Store header.
	var header NavMeshSetHeader
	header.Magic = NAVMESHSET_MAGIC
	header.Version = NAVMESHSET_VERSION
	header.NumTiles = 0
	for i := 0; i < mesh.GetMaxTiles(); i++ {
		tile := mesh.GetTile(i)
		if tile == nil || tile.Header == nil || tile.DataSize == 0 {
			continue
		}
		header.NumTiles++
	}
	header.Params = *mesh.GetParams()
	header.Encode(writer)

	// Store tiles.
	for i := 0; i < mesh.GetMaxTiles(); i++ {
		tile := mesh.GetTile(i)
		if tile == nil || tile.Header == nil || tile.DataSize == 0 {
			continue
		}

		var tileHeader NavMeshTileHeader
		tileHeader.TileRef = mesh.GetTileRef(tile)
		tileHeader.DataSize = tile.DataSize
		tileHeader.Encode(writer)

		fwrite(tile.Data, tile.dataSize, 1, fp)
	}
	writer.Close()
}

func (s *Sample)loadAll(reader fyne.URIReadCloser) *detour.DtNavMesh {

// Read header.
var header NavMeshSetHeader ;
header.Encode(reader)
if header.Magic != NAVMESHSET_MAGIC {
	slog.Error("header.Magic != NAVMESHSET_MAGIC")
	return nil
}
if header.Version != NAVMESHSET_VERSION {
	slog.Error("header.Version != NAVMESHSET_VERSION ")
	return nil
}

var mesh *detour.DtNavMesh

 status := mesh.Init(header.Params);
if (recast.DtStatusFailed(status)) {
fclose(fp);
return 0;
}

// Read tiles.
for i := 0; i < header.NumTiles; i++{
var tileHeader NavMeshTileHeader ;
tileHeader.Decode(reader)

if (!tileHeader.TileRef || !tileHeader.dataSize){
	break;
}


unsigned char* data = (unsigned char*)dtAlloc(tileHeader.dataSize, DT_ALLOC_PERM);
if (!data) break;
memset(data, 0, tileHeader.dataSize);
readLen = fread(data, tileHeader.dataSize, 1, fp);
if (readLen != 1) {
dtFree(data);
fclose(fp);
return 0;
}

mesh.AddTile(data, tileHeader.dataSize, detour.DT_TILE_FREE_DATA, tileHeader.TileRef, 0);
}



return mesh;
}