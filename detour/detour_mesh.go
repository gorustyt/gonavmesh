package detour

import (
	"github.com/gorustyt/gonavmesh/common/rw"
	"math"
	"unsafe"
)

const (
	/// The maximum number of vertices per navigation polygon.
	/// @ingroup detour
	DT_VERTS_PER_POLYGON = 6
	DT_NULL_LINK         = 0xffffffff

	/// A flag that indicates that an entity links to an external entity.
	/// (E.g. A polygon edge is a portal that links to another polygon.)
	DT_EXT_LINK                   = 0x8000
	DT_RAY_CAST_LIMIT_PROPORTIONS = 50.0
	/// A flag that indicates that an off-mesh connection can be traversed in both directions. (Is bidirectional.)
	DT_OFFMESH_CON_BIDIR = 1
	/// @{
	/// @name Tile Serialization Constants
	/// These constants are used to detect whether a navigation tile's data
	/// and state format is compatible with the current build.
	///

	/// A magic number used to detect compatibility of navigation tile data.
	DT_NAVMESH_MAGIC = 'D'<<24 | 'N'<<16 | 'A'<<8 | 'V'

	/// A version number used to detect compatibility of navigation tile data.
	DT_NAVMESH_VERSION = 7

	/// A magic number used to detect the compatibility of navigation tile states.
	DT_NAVMESH_STATE_MAGIC = 'D'<<24 | 'N'<<16 | 'M'<<8 | 'S'

	/// A version number used to detect compatibility of navigation tile states.
	DT_NAVMESH_STATE_VERSION = 1
	/// The maximum number of user defined area ids.
	/// @ingroup detour
	DT_MAX_AREAS = 64
)

const (
	/// Flags representing the type of a navigation mesh polygon.
	/// The polygon is a standard convex polygon that is part of the surface of the mesh.
	DT_POLYTYPE_GROUND = 0
	/// The polygon is an off-mesh connection consisting of two vertices.
	DT_POLYTYPE_OFFMESH_CONNECTION = 1
)

const (
	DT_DETAIL_EDGE_BOUNDARY = 0x01 ///< Detail triangle edge is part of the poly boundary

)

const (
	DT_RAYCAST_USE_COSTS = 0x01 ///< Raycast should calculate movement cost along the ray and fill RaycastHit::cost

)
const (

	/// Tile flags used for various functions and fields.
	/// For an example, see DtNavMesh::addTile().

	/// The navigation mesh owns the tile memory and is responsible for freeing it.
	DT_TILE_FREE_DATA = 0x01
)

var (
	// Undefine (or define in a build config) the following line to use 64bit polyref.
	// Generally not needed, useful for very large worlds.
	// Note: tiles build using 32bit refs are not compatible with 64bit refs!
	DT_POLYREF64 = 1
	DT_SALT_BITS = 16
	DT_TILE_BITS = 28
	DT_POLY_BITS = 20
)

func InitMesh(dT_POLYREF64 int) {
	DT_POLYREF64 = dT_POLYREF64

}

type DtPolyRef int
type DtTileRef int

// / Defines a polygon within a DtMeshTile object.
// / @ingroup detour
type DtPoly struct {
	/// Index to first link in linked list. (Or #DT_NULL_LINK if there is no link.)
	FirstLink uint32

	/// The indices of the polygon's vertices.
	/// The actual vertices are located in DtMeshTile::verts.
	Verts [DT_VERTS_PER_POLYGON]uint16

	/// Packed data representing neighbor polygons references and flags for each edge.
	Neis [DT_VERTS_PER_POLYGON]uint16

	/// The user defined polygon flags.
	Flags uint16

	/// The number of vertices in the polygon.
	VertCount uint8

	/// The bit packed area id and polygon type.
	/// @note Use the structure's set and get methods to access this value.
	AreaAndtype uint8
}

func (d *DtPoly) ToBin(w *rw.ReaderWriter) {
	w.WriteInt32(d.FirstLink)
	w.WriteInt16s(d.Verts[:])
	w.WriteInt16s(d.Neis[:])
	w.WriteInt16(d.Flags)
	w.WriteInt8(d.VertCount)
	w.WriteInt8(d.AreaAndtype)
}

func (d *DtPoly) FromBin(w *rw.ReaderWriter) *DtPoly {
	d.FirstLink = w.ReadUInt32()
	d.Flags = w.ReadUInt16()
	w.ReadUInt16s(d.Verts[:])
	w.ReadUInt16s(d.Neis[:])
	d.VertCount = w.ReadUInt8()
	d.AreaAndtype = w.ReadUInt8()
	return d
}

// / Sets the user defined area id. [Limit: < #DT_MAX_AREAS]
func (p *DtPoly) SetArea(a uint8) { p.AreaAndtype = (p.AreaAndtype & 0xc0) | (a & 0x3f) }

// / Sets the polygon type. (See: #dtPolyTypes.)
func (p *DtPoly) SetType(t uint8) { p.AreaAndtype = (p.AreaAndtype & 0x3f) | (t << 6) }

// / Gets the user defined area id.
func (p *DtPoly) GetArea() uint8 { return p.AreaAndtype & 0x3f }

// / Gets the polygon type. (See: #dtPolyTypes)
func (p *DtPoly) GetType() uint8 { return p.AreaAndtype >> 6 }

// / Defines the location of detail sub-mesh data within a DtMeshTile.
type DtPolyDetail struct {
	VertBase  uint32 ///< The offset of the vertices in the DtMeshTile::detailVerts array.
	TriBase   uint32 ///< The offset of the triangles in the DtMeshTile::detailTris array.
	VertCount uint8  ///< The number of vertices in the sub-mesh.
	TriCount  uint8  ///< The number of triangles in the sub-mesh.
}

func (d *DtPolyDetail) ToBin(w *rw.ReaderWriter) {
	w.WriteInt32(d.VertBase)
	w.WriteInt32(d.TriBase)
	w.WriteInt8(d.VertCount)
	w.WriteInt8(d.TriCount)
}

func (d *DtPolyDetail) FromBin(w *rw.ReaderWriter) *DtPolyDetail {
	d.VertBase = w.ReadUInt32()
	d.TriBase = w.ReadUInt32()
	d.VertCount = w.ReadUInt8()
	d.TriCount = w.ReadUInt8()
	return d
}

// Defines a link between polygons.
// / @note This structure is rarely if ever used by the end user.
// / @see DtMeshTile
type DtLink struct {
	Ref  DtPolyRef ///< Neighbour reference. (The neighbor that is linked to.)
	Next uint32    ///< Index of the next link.
	Edge uint8     ///< Index of the polygon edge that owns this link.
	Side uint8     ///< If a boundary link, defines on which side the link is.
	Bmin uint8     ///< If a boundary link, defines the minimum sub-edge area.
	Bmax uint8     ///< If a boundary link, defines the maximum sub-edge area.
}

func (d *DtLink) ToBin(w *rw.ReaderWriter) {
	w.WriteInt32(int32(d.Ref))
	w.WriteInt8(d.Next)
	w.WriteInt8(d.Edge)
	w.WriteInt8(d.Side)
	w.WriteInt8(d.Bmin)
	w.WriteInt8(d.Bmax)
}

func (d *DtLink) FromBin(w *rw.ReaderWriter) *DtLink {
	d.Ref = DtPolyRef(w.ReadInt32())
	d.Next = w.ReadUInt32()
	d.Edge = w.ReadUInt8()
	d.Side = w.ReadUInt8()
	d.Bmin = w.ReadUInt8()
	d.Bmax = w.ReadUInt8()
	return d
}

// / Bounding volume node.
// / @note This structure is rarely if ever used by the end user.
// / @see DtMeshTile
type DtBVNode struct {
	Bmin [3]uint16 ///< Minimum bounds of the node's AABB. [(x, y, z)]
	Bmax [3]uint16 ///< Maximum bounds of the node's AABB. [(x, y, z)]
	I    int32     ///< The node's index. (Negative for escape sequence.)
}

func (d *DtBVNode) ToBin(w *rw.ReaderWriter) {
	w.WriteInt32(d.Bmin[:])
	w.WriteInt16s(d.Bmax[:])
	w.WriteInt32(d.I)
}

func (d *DtBVNode) FromBin(w *rw.ReaderWriter) *DtBVNode {
	w.ReadUInt16s(d.Bmin[:])
	w.ReadUInt16s(d.Bmax[:])
	d.I = w.ReadInt32()
	return d
}

// / Defines a navigation mesh tile.
// / @ingroup detour
type DtMeshTile struct {
	salt uint32 ///< Counter describing modifications to the tile.

	linksFreeList uint32          ///< Index to the next free link.
	Header        *DtMeshHeader   ///< The tile header.
	Polys         []*DtPoly       ///< The tile polygons. [Size: DtMeshHeader::polyCount]
	Verts         []float32       ///< The tile vertices. [(x, y, z) * DtMeshHeader::vertCount]
	Links         []*DtLink       ///< The tile links. [Size: DtMeshHeader::maxLinkCount]
	DetailMeshes  []*DtPolyDetail ///< The tile's detail sub-meshes. [Size: DtMeshHeader::detailMeshCount]

	/// The detail mesh's unique vertices. [(x, y, z) * DtMeshHeader::detailVertCount]
	DetailVerts []float32

	/// The detail mesh's triangles. [(vertA, vertB, vertC, triFlags) * DtMeshHeader::detailTriCount].
	/// See dtDetailTriEdgeFlags and dtGetDetailTriEdgeFlags.
	DetailTris []uint8

	/// The tile bounding volume nodes. [Size: DtMeshHeader::bvNodeCount]
	/// (Will be null if bounding volumes are disabled.)
	BvTree []*DtBVNode

	OffMeshCons []*DtOffMeshConnection ///< The tile off-mesh connections. [Size: DtMeshHeader::offMeshConCount]
	Flags       int32                  ///< Tile flags. (See: #dtTileFlags)
	Next        *DtMeshTile            ///< The next free tile, or the next tile in the spatial grid.
	//存储的数据
	Data *NavMeshData
}

func (d *DtMeshTile) GetIndexByPloy(ploy *DtPoly) int {
	for i, v := range d.Polys {
		if v == ploy {
			return i
		}
	}
	return -1
}

// / Provides high level information related to a DtMeshTile object.
// / @ingroup detour
type DtMeshHeader struct {
	Magic           int32  ///< Tile magic number. (Used to identify the data format.)
	Version         int32  ///< Tile data format version number.
	X               int32  ///< The x-position of the tile within the DtNavMesh tile grid. (x, y, layer)
	Y               int32  ///< The y-position of the tile within the DtNavMesh tile grid. (x, y, layer)
	Layer           int32  ///< The layer of the tile within the DtNavMesh tile grid. (x, y, layer)
	UserId          uint32 ///< The user defined id of the tile.
	PolyCount       int32  ///< The number of polygons in the tile.
	VertCount       int32  ///< The number of vertices in the tile.
	MaxLinkCount    int32  ///< The number of allocated links.
	DetailMeshCount int32  ///< The number of sub-meshes in the detail mesh.

	/// The number of unique vertices in the detail mesh. (In addition to the polygon vertices.)
	DetailVertCount int32

	DetailTriCount  int32      ///< The number of triangles in the detail mesh.
	BvNodeCount     int32      ///< The number of bounding volume nodes. (Zero if bounding volumes are disabled.)
	OffMeshConCount int32      ///< The number of off-mesh connections.
	OffMeshBase     int32      ///< The index of the first polygon which is an off-mesh connection.
	WalkableHeight  float32    ///< The height of the agents using the tile.
	WalkableRadius  float32    ///< The radius of the agents using the tile.
	WalkableClimb   float32    ///< The maximum climb height of the agents using the tile.
	Bmin            [3]float32 ///< The minimum bounds of the tile's AABB. [(x, y, z)]
	Bmax            [3]float32 ///< The maximum bounds of the tile's AABB. [(x, y, z)]

	/// The bounding volume quantization factor.
	BvQuantFactor float32
}

func (d *DtMeshHeader) Size() int {
	return int(unsafe.Sizeof(d))
}

func (d *DtMeshHeader) ToBin(w *rw.ReaderWriter) {
	w.WriteInt32(d.Magic)
	w.WriteInt32(d.Version)
	w.WriteInt32(d.X)
	w.WriteInt32(d.Y)
	w.WriteInt32(d.Layer)
	w.WriteInt32(d.UserId)
	w.WriteInt32(d.PolyCount)
	w.WriteInt32(d.VertCount)
	w.WriteInt32(d.MaxLinkCount)
	w.WriteInt32(d.DetailMeshCount)
	w.WriteInt32(d.DetailVertCount)
	w.WriteInt32(d.DetailTriCount)
	w.WriteInt32(d.BvNodeCount)
	w.WriteInt32(d.OffMeshConCount)
	w.WriteInt32(d.OffMeshBase)
	w.WriteFloat32(d.WalkableHeight)
	w.WriteFloat32(d.WalkableRadius)
	w.WriteFloat32(d.WalkableClimb)
	w.WriteFloat32(d.Bmin[:])
	w.WriteFloat32(d.Bmax[:])
	w.WriteFloat32(d.BvQuantFactor)
}

func (d *DtMeshHeader) FromBin(w *rw.ReaderWriter) *DtMeshHeader {
	d.Magic = w.ReadInt32()
	d.Version = w.ReadInt32()
	d.X = w.ReadInt32()
	d.Y = w.ReadInt32()
	d.Layer = w.ReadInt32()
	d.UserId = w.ReadUInt32()
	d.PolyCount = w.ReadInt32()
	d.VertCount = w.ReadInt32()
	d.MaxLinkCount = w.ReadInt32()
	d.DetailMeshCount = w.ReadInt32()
	d.DetailVertCount = w.ReadInt32()
	d.DetailTriCount = w.ReadInt32()
	d.BvNodeCount = w.ReadInt32()
	d.OffMeshConCount = w.ReadInt32()
	d.OffMeshBase = w.ReadInt32()
	d.WalkableHeight = w.ReadFloat32()
	d.WalkableRadius = w.ReadFloat32()
	d.WalkableClimb = w.ReadFloat32()
	w.ReadFloat32s(d.Bmin[:])
	w.ReadFloat32s(d.Bmax[:])
	d.BvQuantFactor = w.ReadFloat32()
	return d
}

// / Defines an navigation mesh off-mesh connection within a DtMeshTile object.
// / An off-mesh connection is a user defined traversable connection made up to two vertices.
type DtOffMeshConnection struct {
	/// The endpoints of the connection. [(ax, ay, az, bx, by, bz)]
	Pos [6]float32

	/// The radius of the endpoints. [Limit: >= 0]
	Rad float32

	/// The polygon reference of the connection within the tile.
	Poly uint16

	/// Link flags.
	/// @note These are not the connection's user defined flags. Those are assigned via the
	/// connection's DtPoly definition. These are link flags used for internal purposes.
	Flags uint8

	/// End point side.
	Side uint8

	/// The id of the offmesh connection. (User assigned when the navigation mesh is built.)
	UserId uint32
}

func (d *DtOffMeshConnection) ToBin(w *rw.ReaderWriter) {
	w.WriteFloat32s(d.Pos)
	w.WriteFloat32(d.Rad)
	w.WriteInt16(d.Poly)
	w.WriteInt8(d.Flags)
	w.WriteInt8(d.Side)
	w.WriteInt32(d.UserId)
}

func (d *DtOffMeshConnection) FromBin(w *rw.ReaderWriter) *DtOffMeshConnection {
	w.ReadFloat32s(d.Pos[:])
	d.Rad = w.ReadFloat32()
	d.Poly = w.ReadUInt16()
	d.Flags = w.ReadUInt8()
	d.Side = w.ReadUInt8()
	d.UserId = w.ReadUInt32()
	return d
}

// / Configuration parameters used to define multi-tile navigation meshes.
// / The values are used to allocate space during the initialization of a navigation mesh.
// / @see DtNavMesh::init()
// / @ingroup detour
type NavMeshParams struct {
	Orig       [3]float32 ///< The world space origin of the navigation mesh's tile space. [(x, y, z)]
	TileWidth  float32    ///< The width of each tile. (Along the x-axis.)
	TileHeight float32    ///< The height of each tile. (Along the z-axis.)
	MaxTiles   int32      ///< The maximum number of tiles the navigation mesh can contain. This and maxPolys are used to calculate how many bits are needed to identify tiles and polygons uniquely.
	MaxPolys   int32      ///< The maximum number of polygons each tile can contain. This and maxTiles are used to calculate how many bits are needed to identify tiles and polygons uniquely.
}

func (d *NavMeshParams) FromBin(r *rw.ReaderWriter) {
	r.ReadFloat32s(d.Orig[:])
	d.TileWidth = r.ReadFloat32()
	d.TileHeight = r.ReadFloat32()
	d.MaxTiles = r.ReadInt32()
	d.MaxPolys = r.ReadInt32()
}

func (d *NavMeshParams) ToBin(w *rw.ReaderWriter) {
	w.WriteFloat32s(d.Orig[:])
	w.WriteFloat32(d.TileWidth)
	w.WriteFloat32(d.TileHeight)
	w.WriteInt32(d.MaxTiles)
	w.WriteInt32(d.MaxPolys)
}

type IDtNavMesh interface {
	GetParams() *NavMeshParams
	/// Adds a tile to the navigation mesh.
	///  @param[in]		data		Data for the new tile mesh. (See: #dtCreateNavMeshData)
	///  @param[in]		dataSize	Data size of the new tile mesh.
	///  @param[in]		flags		Tile flags. (See: #dtTileFlags)
	///  @param[in]		lastRef		The desired reference for the tile. (When reloading a tile.) [opt] [Default: 0]
	///  @param[out]	result		The tile reference. (If the tile was succesfully added.) [opt]
	/// @return The status flags for the operation.
	AddTile(data *NavMeshData, flags int32, lastRef DtTileRef) (result DtTileRef, status DtStatus)
	/// Removes the specified tile from the navigation mesh.
	///  @param[in]		ref			The reference of the tile to remove.
	/// @return The status flags for the operation.
	RemoveTile(ref DtTileRef) (data *NavMeshData, status DtStatus)
	/// @name Query Functions

	/// Calculates the tile grid location for the specified world position.
	///  @param[in]	pos  The world position for the query. [(x, y, z)]
	///  @param[out]	tx		The tile's x-location. (x, y)
	///  @param[out]	ty		The tile's y-location. (x, y)
	CalcTileLoc(pos []float32) (tx, ty int32)
	/// Gets the tile at the specified grid location.
	///  @param[in]	x		The tile's x-location. (x, y, layer)
	///  @param[in]	y		The tile's y-location. (x, y, layer)
	///  @param[in]	layer	The tile's layer. (x, y, layer)
	/// @return The tile, or null if the tile does not exist.
	GetTileAt(x, y, layer int32) *DtMeshTile
	/// Gets all tiles at the specified grid location. (All layers.)
	///  @param[in]		x			The tile's x-location. (x, y)
	///  @param[in]		y			The tile's y-location. (x, y)
	///  @param[out]	tiles		A pointer to an array of tiles that will hold the result.
	///  @param[in]		maxTiles	The maximum tiles the tiles parameter can hold.
	/// @return The number of tiles returned in the tiles array.
	GetTilesAt(x, y int32, maxTiles int32) ([]*DtMeshTile, int32)
	/// Gets the tile reference for the tile at specified grid location.
	///  @param[in]	x		The tile's x-location. (x, y, layer)
	///  @param[in]	y		The tile's y-location. (x, y, layer)
	///  @param[in]	layer	The tile's layer. (x, y, layer)
	/// @return The tile reference of the tile, or 0 if there is none.
	GetTileRefAt(x, y, layer int32) DtTileRef
	/// Gets the tile reference for the specified tile.
	///  @param[in]	tile	The tile.
	/// @return The tile reference of the tile.
	GetTileRef(tile *DtMeshTile) DtTileRef
	/// Gets the tile for the specified tile reference.
	///  @param[in]	ref		The tile reference of the tile to retrieve.
	/// @return The tile for the specified reference, or null if the
	///		reference is invalid.
	GetTileByRef(ref DtTileRef) *DtMeshTile
	/// The maximum number of tiles supported by the navigation mesh.
	/// @return The maximum number of tiles supported by the navigation mesh.
	GetMaxTiles() int32
	/// Gets the tile at the specified index.
	///  @param[in]	i		The tile index. [Limit: 0 >= index < #getMaxTiles()]
	/// @return The tile at the specified index.
	GetTile(i int) *DtMeshTile
	/// Gets the tile and polygon for the specified polygon reference.
	///  @param[in]		ref		The reference for the a polygon.
	///  @param[out]	tile	The tile containing the polygon.
	///  @param[out]	poly	The polygon.
	/// @return The status flags for the operation.
	GetTileAndPolyByRef(ref DtPolyRef) (tile *DtMeshTile, poly *DtPoly, status DtStatus)
	/// Returns the tile and polygon for the specified polygon reference.
	///  @param[in]		ref		A known valid reference for a polygon.
	///  @param[out]	tile	The tile containing the polygon.
	///  @param[out]	poly	The polygon.
	GetTileAndPolyByRefUnsafe(ref DtPolyRef) (tile *DtMeshTile, poly *DtPoly)
	/// Checks the validity of a polygon reference.
	///  @param[in]	ref		The polygon reference to check.
	/// @return True if polygon reference is valid for the navigation mesh.
	IsValidPolyRef(ref DtPolyRef) bool
	/// Gets the polygon reference for the tile's base polygon.
	///  @param[in]	tile		The tile.
	/// @return The polygon reference for the base polygon in the specified tile.
	GetPolyRefBase(tile *DtMeshTile) DtPolyRef
	/// Gets the endpoints for an off-mesh connection, ordered by "direction of travel".
	///  @param[in]		prevRef		The reference of the polygon before the connection.
	///  @param[in]		polyRef		The reference of the off-mesh connection polygon.
	///  @param[out]	startPos	The start position of the off-mesh connection. [(x, y, z)]
	///  @param[out]	endPos		The end position of the off-mesh connection. [(x, y, z)]
	/// @return The status flags for the operation.
	GetOffMeshConnectionPolyEndPoints(prevRef DtPolyRef, polyRef DtPolyRef, startPos []float32, endPos []float32) DtStatus
	/// Gets the specified off-mesh connection.
	///  @param[in]	ref		The polygon reference of the off-mesh connection.
	/// @return The specified off-mesh connection, or null if the polygon reference is not valid
	GetOffMeshConnectionByRef(ref DtPolyRef) *DtOffMeshConnection
	/// @}

	/// @{
	/// @name State Management
	/// These functions do not effect #DtTileRef or #DtPolyRef's.

	/// Sets the user defined flags for the specified polygon.
	///  @param[in]	ref		The polygon reference.
	///  @param[in]	flags	The new flags for the polygon.
	/// @return The status flags for the operation.
	SetPolyFlags(ref DtPolyRef, flags uint16) DtStatus

	/// Gets the user defined flags for the specified polygon.
	///  @param[in]		ref				The polygon reference.
	///  @param[out]	resultFlags		The polygon flags.
	/// @return The status flags for the operation.
	GetPolyFlags(ref DtPolyRef) (resultFlags uint16, status DtStatus)

	/// Sets the user defined area for the specified polygon.
	///  @param[in]	ref		The polygon reference.
	///  @param[in]	area	The new area id for the polygon. [Limit: < #DT_MAX_AREAS]
	/// @return The status flags for the operation.
	SetPolyArea(ref DtPolyRef, area uint8) DtStatus
	/// Gets the user defined area for the specified polygon.
	///  @param[in]		ref			The polygon reference.
	///  @param[out]	resultArea	The area id for the polygon.
	/// @return The status flags for the operation.
	GetPolyArea(ref DtPolyRef) (resultArea uint8, status DtStatus)

	/// Gets the size of the buffer required by #storeTileState to store the specified tile's state.
	///  @param[in]	tile	The tile.
	/// @return The size of the buffer required to store the state.

	/// Restores the state of the tile.
	///  @param[in]	tile			The tile.
	///  @param[in]	data			The new state. (Obtained from #storeTileState.)
	///  @param[in]	maxDataSize		The size of the state within the data buffer.
	/// @return The status flags for the operation.
	RestoreTileState(tile *DtMeshTile, tileState *dtTileState, polyStates []*dtPolyState, maxDataSize int32) DtStatus
	/// @name Encoding and Decoding
	/// These functions are generally meant for internal use only.

	/// Derives a standard polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	salt	The tile's salt value.
	///  @param[in]	it		The index of the tile.
	///  @param[in]	ip		The index of the polygon within the tile.
	EncodePolyId(salt, it, ip uint32) DtPolyRef
	/// Decodes a standard polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	ref   The polygon reference to decode.
	///  @param[out]	salt	The tile's salt value.
	///  @param[out]	it		The index of the tile.
	///  @param[out]	ip		The index of the polygon within the tile.
	///  @see #encodePolyId
	DecodePolyId(ref DtPolyRef) (salt, it, ip uint32)
	/// Extracts a tile's salt value from the specified polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	ref		The polygon reference.
	///  @see #encodePolyId
	DecodePolyIdPoly(ref DtPolyRef) uint32
	/// Extracts the tile's index from the specified polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	ref		The polygon reference.
	///  @see #encodePolyId
	DecodePolyIdTile(ref DtPolyRef) uint32
	/// Extracts a tile's salt value from the specified polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	ref		The polygon reference.
	///  @see #encodePolyId
	DecodePolyIdSalt(ref DtPolyRef) uint32

	GetPolyHeight(tile *DtMeshTile, poly *DtPoly, pos []float32) (height float32, ok bool)
	ClosestPointOnPoly(ref DtPolyRef, pos []float32, closest []float32, posOverPoly *bool)
}

type DtNavMesh struct {
	m_params                  *NavMeshParams ///< Current initialization params. TODO: do not store this info twice.
	m_orig                    [3]float32     ///< Origin of the tile (0,0)
	m_tileWidth, m_tileHeight float32        ///< Dimensions of each tile.
	m_maxTiles                int32          ///< Max number of tiles.
	m_tileLutSize             int32          ///< Tile hash lookup size (must be pot).
	m_tileLutMask             int32          ///< Tile hash lookup mask.
	m_posLookup               []*DtMeshTile  ///< Tile hash lookup.
	m_nextFree                *DtMeshTile    ///< Freelist of tiles.
	m_tiles                   []*DtMeshTile  ///< List of tiles.

	m_saltBits uint32 ///< Number of salt bits in the tile ID.
	m_tileBits uint32 ///< Number of tile bits in the tile ID.
	m_polyBits uint32 ///< Number of poly bits in the tile ID.
}

func (mesh *DtNavMesh) GetTile(i int) *DtMeshTile {
	return mesh.m_tiles[i]
}
func (mesh *DtNavMesh) GetMaxTiles() int32 {
	return mesh.m_maxTiles
}

/// @{
/// @name Encoding and Decoding
/// These functions are generally meant for internal use only.

// / Derives a standard polygon reference.
// /  @note This function is generally meant for internal use only.
// /  @param[in]	salt	The tile's salt value.
// /  @param[in]	it		The index of the tile.
// /  @param[in]	ip		The index of the polygon within the tile.
func (mesh *DtNavMesh) EncodePolyId(salt, it, ip uint32) DtPolyRef {
	if DT_POLYREF64 == 1 {
		return DtPolyRef(int64(salt)<<(DT_POLY_BITS+DT_TILE_BITS) | int64(it<<DT_POLY_BITS) | int64(ip))
	} else {
		return DtPolyRef(salt<<(mesh.m_polyBits+mesh.m_tileBits) | (it << mesh.m_polyBits) | ip)
	}
}

func (mesh *DtNavMesh) DecodePolyIdTile(ref DtPolyRef) uint32 {
	if DT_POLYREF64 == 1 {
		tileMask := (1 << DT_TILE_BITS) - 1
		return uint32((int(ref) >> DT_POLY_BITS) & tileMask)
	} else {
	}
	tileMask := (1 << mesh.m_tileBits) - 1
	return uint32((int(ref) >> mesh.m_polyBits) & tileMask)
}

// / Decodes a standard polygon reference.
// /  @note This function is generally meant for internal use only.
// /  @param[in]	ref   The polygon reference to decode.
// /  @param[out]	salt	The tile's salt value.
// /  @param[out]	it		The index of the tile.
// /  @param[out]	ip		The index of the polygon within the tile.
// /  @see #encodePolyId
func (mesh *DtNavMesh) DecodePolyId(r DtPolyRef) (salt, it, ip uint32) {
	ref := int(r)
	if DT_POLYREF64 == 1 {
		saltMask := (1 << DT_SALT_BITS) - 1
		tileMask := (1 << DT_TILE_BITS) - 1
		polyMask := (1 << DT_POLY_BITS) - 1
		salt = uint32((ref >> (DT_POLY_BITS + DT_TILE_BITS)) & saltMask)
		it = uint32((ref >> DT_POLY_BITS) & tileMask)
		ip = uint32(ref & polyMask)
	} else {
		saltMask := (1 << mesh.m_saltBits) - 1
		tileMask := (1 << mesh.m_tileBits) - 1
		polyMask := (1 << mesh.m_polyBits) - 1
		salt = uint32((ref >> (mesh.m_polyBits + mesh.m_tileBits)) & saltMask)
		it = uint32((ref >> mesh.m_polyBits) & tileMask)
		ip = uint32(ref & polyMask)
	}
	return
}

// / Get flags for edge in detail triangle.
// / @param[in]	triFlags		The flags for the triangle (last component of detail vertices above).
// / @param[in]	edgeIndex		The index of the first vertex of the edge. For instance, if 0,
// /								returns flags for edge AB.
func DtGetDetailTriEdgeFlags(triFlags uint8, edgeIndex int32) int32 {
	return int32(triFlags>>(edgeIndex*2)) & 0x3
}

func (mesh *DtNavMesh) DecodePolyIdPoly(ref DtPolyRef) uint32 {
	if DT_POLYREF64 == 1 {
		polyMask := (1 << DT_POLY_BITS) - 1
		return uint32(int(ref) & polyMask)
	} else {
		polyMask := (1 << mesh.m_polyBits) - 1
		return uint32(int(ref) & polyMask)
	}
}

// / Extracts a tile's salt value from the specified polygon reference.
// /  @note This function is generally meant for internal use only.
// /  @param[in]	ref		The polygon reference.
// /  @see #encodePolyId
func (mesh *DtNavMesh) DecodePolyIdSalt(ref DtPolyRef) uint32 {
	if DT_POLYREF64 == 1 {
		saltMask := (1 << DT_SALT_BITS) - 1
		return uint32((int(ref) >> (DT_POLY_BITS + DT_TILE_BITS)) & saltMask)
	}
	saltMask := (1 << mesh.m_saltBits) - 1
	return uint32((int(ref) >> (mesh.m_polyBits + mesh.m_tileBits)) & saltMask)

}

func (mesh *DtNavMesh) CalcTileLoc(pos []float32) (tx, ty int32) {
	tx = int32(math.Floor(float64((pos[0] - mesh.m_orig[0]) / mesh.m_tileWidth)))
	ty = int32(math.Floor(float64((pos[2] - mesh.m_orig[2]) / mesh.m_tileHeight)))
	return tx, ty
}

func (mesh *DtNavMesh) GetTileByRef(ref DtTileRef) *DtMeshTile {
	if ref == 0 {
		return nil
	}

	tileIndex := mesh.DecodePolyIdTile(DtPolyRef(ref))
	tileSalt := mesh.DecodePolyIdSalt(DtPolyRef(ref))
	if int(tileIndex) >= int(mesh.m_maxTiles) {
		return nil
	}

	tile := mesh.m_tiles[tileIndex]
	if tile.salt != tileSalt {
		return nil
	}

	return tile
}
