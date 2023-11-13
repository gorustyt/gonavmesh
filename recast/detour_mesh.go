package recast

import "math"

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
	/// For an example, see dtNavMesh::addTile().

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

type dtPolyRef int
type dtTileRef int

// / Defines a polygon within a dtMeshTile object.
// / @ingroup detour
type dtPoly struct {
	/// Index to first link in linked list. (Or #DT_NULL_LINK if there is no link.)
	firstLink int

	/// The indices of the polygon's vertices.
	/// The actual vertices are located in dtMeshTile::verts.
	verts [DT_VERTS_PER_POLYGON]int

	/// Packed data representing neighbor polygons references and flags for each edge.
	neis [DT_VERTS_PER_POLYGON]int

	/// The user defined polygon flags.
	flags int

	/// The number of vertices in the polygon.
	vertCount int

	/// The bit packed area id and polygon type.
	/// @note Use the structure's set and get methods to access this value.
	areaAndtype int
}

// / Sets the user defined area id. [Limit: < #DT_MAX_AREAS]
func (p *dtPoly) setArea(a int) { p.areaAndtype = (p.areaAndtype & 0xc0) | (a & 0x3f) }

// / Sets the polygon type. (See: #dtPolyTypes.)
func (p *dtPoly) setType(t int) { p.areaAndtype = (p.areaAndtype & 0x3f) | (t << 6) }

// / Gets the user defined area id.
func (p *dtPoly) getArea() int { return p.areaAndtype & 0x3f }

// / Gets the polygon type. (See: #dtPolyTypes)
func (p *dtPoly) getType() int { return p.areaAndtype >> 6 }

// / Defines the location of detail sub-mesh data within a dtMeshTile.
type dtPolyDetail struct {
	vertBase  int ///< The offset of the vertices in the dtMeshTile::detailVerts array.
	triBase   int ///< The offset of the triangles in the dtMeshTile::detailTris array.
	vertCount int ///< The number of vertices in the sub-mesh.
	triCount  int ///< The number of triangles in the sub-mesh.
}

// Defines a link between polygons.
// / @note This structure is rarely if ever used by the end user.
// / @see dtMeshTile
type dtLink struct {
	ref  dtPolyRef ///< Neighbour reference. (The neighbor that is linked to.)
	next int       ///< Index of the next link.
	edge int       ///< Index of the polygon edge that owns this link.
	side int       ///< If a boundary link, defines on which side the link is.
	bmin int       ///< If a boundary link, defines the minimum sub-edge area.
	bmax int       ///< If a boundary link, defines the maximum sub-edge area.
}

// / Bounding volume node.
// / @note This structure is rarely if ever used by the end user.
// / @see dtMeshTile
type dtBVNode struct {
	bmin [3]int ///< Minimum bounds of the node's AABB. [(x, y, z)]
	bmax [3]int ///< Maximum bounds of the node's AABB. [(x, y, z)]
	i    int    ///< The node's index. (Negative for escape sequence.)
}

// / Defines a navigation mesh tile.
// / @ingroup detour
type dtMeshTile struct {
	salt int ///< Counter describing modifications to the tile.

	linksFreeList int             ///< Index to the next free link.
	header        *dtMeshHeader   ///< The tile header.
	polys         []*dtPoly       ///< The tile polygons. [Size: dtMeshHeader::polyCount]
	verts         []float64       ///< The tile vertices. [(x, y, z) * dtMeshHeader::vertCount]
	links         []*dtLink       ///< The tile links. [Size: dtMeshHeader::maxLinkCount]
	detailMeshes  []*dtPolyDetail ///< The tile's detail sub-meshes. [Size: dtMeshHeader::detailMeshCount]

	/// The detail mesh's unique vertices. [(x, y, z) * dtMeshHeader::detailVertCount]
	detailVerts []float64

	/// The detail mesh's triangles. [(vertA, vertB, vertC, triFlags) * dtMeshHeader::detailTriCount].
	/// See dtDetailTriEdgeFlags and dtGetDetailTriEdgeFlags.
	detailTris []int

	/// The tile bounding volume nodes. [Size: dtMeshHeader::bvNodeCount]
	/// (Will be null if bounding volumes are disabled.)
	bvTree []*dtBVNode

	offMeshCons []*dtOffMeshConnection ///< The tile off-mesh connections. [Size: dtMeshHeader::offMeshConCount]

	data     []int       ///< The tile data. (Not directly accessed under normal situations.)
	dataSize int         ///< Size of the tile data.
	flags    int         ///< Tile flags. (See: #dtTileFlags)
	next     *dtMeshTile ///< The next free tile, or the next tile in the spatial grid.
}

// / Provides high level information related to a dtMeshTile object.
// / @ingroup detour
type dtMeshHeader struct {
	magic           int ///< Tile magic number. (Used to identify the data format.)
	version         int ///< Tile data format version number.
	x               int ///< The x-position of the tile within the dtNavMesh tile grid. (x, y, layer)
	y               int ///< The y-position of the tile within the dtNavMesh tile grid. (x, y, layer)
	layer           int ///< The layer of the tile within the dtNavMesh tile grid. (x, y, layer)
	userId          int ///< The user defined id of the tile.
	polyCount       int ///< The number of polygons in the tile.
	vertCount       int ///< The number of vertices in the tile.
	maxLinkCount    int ///< The number of allocated links.
	detailMeshCount int ///< The number of sub-meshes in the detail mesh.

	/// The number of unique vertices in the detail mesh. (In addition to the polygon vertices.)
	detailVertCount int

	detailTriCount  int        ///< The number of triangles in the detail mesh.
	bvNodeCount     int        ///< The number of bounding volume nodes. (Zero if bounding volumes are disabled.)
	offMeshConCount int        ///< The number of off-mesh connections.
	offMeshBase     int        ///< The index of the first polygon which is an off-mesh connection.
	walkableHeight  float64    ///< The height of the agents using the tile.
	walkableRadius  float64    ///< The radius of the agents using the tile.
	walkableClimb   float64    ///< The maximum climb height of the agents using the tile.
	bmin            [3]float64 ///< The minimum bounds of the tile's AABB. [(x, y, z)]
	bmax            [3]float64 ///< The maximum bounds of the tile's AABB. [(x, y, z)]

	/// The bounding volume quantization factor.
	bvQuantFactor float64
}

// / Defines an navigation mesh off-mesh connection within a dtMeshTile object.
// / An off-mesh connection is a user defined traversable connection made up to two vertices.
type dtOffMeshConnection struct {
	/// The endpoints of the connection. [(ax, ay, az, bx, by, bz)]
	pos [6]float64

	/// The radius of the endpoints. [Limit: >= 0]
	rad float64

	/// The polygon reference of the connection within the tile.
	poly int

	/// Link flags.
	/// @note These are not the connection's user defined flags. Those are assigned via the
	/// connection's dtPoly definition. These are link flags used for internal purposes.
	flags int

	/// End point side.
	side int

	/// The id of the offmesh connection. (User assigned when the navigation mesh is built.)
	userId int
}

// / Configuration parameters used to define multi-tile navigation meshes.
// / The values are used to allocate space during the initialization of a navigation mesh.
// / @see dtNavMesh::init()
// / @ingroup detour
type dtNavMeshParams struct {
	orig       [3]float64 ///< The world space origin of the navigation mesh's tile space. [(x, y, z)]
	tileWidth  float64    ///< The width of each tile. (Along the x-axis.)
	tileHeight float64    ///< The height of each tile. (Along the z-axis.)
	maxTiles   int        ///< The maximum number of tiles the navigation mesh can contain. This and maxPolys are used to calculate how many bits are needed to identify tiles and polygons uniquely.
	maxPolys   int        ///< The maximum number of polygons each tile can contain. This and maxTiles are used to calculate how many bits are needed to identify tiles and polygons uniquely.
}

type DtNavMesh interface {
	/// Adds a tile to the navigation mesh.
	///  @param[in]		data		Data for the new tile mesh. (See: #dtCreateNavMeshData)
	///  @param[in]		dataSize	Data size of the new tile mesh.
	///  @param[in]		flags		Tile flags. (See: #dtTileFlags)
	///  @param[in]		lastRef		The desired reference for the tile. (When reloading a tile.) [opt] [Default: 0]
	///  @param[out]	result		The tile reference. (If the tile was succesfully added.) [opt]
	/// @return The status flags for the operation.
	addTile(header *dtMeshHeader, titleData *dtMeshTile, dataSize int, flags int, lastRef dtTileRef) (result dtTileRef, status dtStatus)
	/// Removes the specified tile from the navigation mesh.
	///  @param[in]		ref			The reference of the tile to remove.
	///  @param[out]	data		Data associated with deleted tile.
	///  @param[out]	dataSize	Size of the data associated with deleted tile.
	/// @return The status flags for the operation.
	removeTile(ref dtTileRef) (data []int, dataSize int, status dtStatus)
	/// @name Query Functions

	/// Calculates the tile grid location for the specified world position.
	///  @param[in]	pos  The world position for the query. [(x, y, z)]
	///  @param[out]	tx		The tile's x-location. (x, y)
	///  @param[out]	ty		The tile's y-location. (x, y)
	calcTileLoc(pos []float64) (tx, ty int)
	/// Gets the tile at the specified grid location.
	///  @param[in]	x		The tile's x-location. (x, y, layer)
	///  @param[in]	y		The tile's y-location. (x, y, layer)
	///  @param[in]	layer	The tile's layer. (x, y, layer)
	/// @return The tile, or null if the tile does not exist.
	getTileAt(x, y, layer int) *dtMeshTile
	/// Gets all tiles at the specified grid location. (All layers.)
	///  @param[in]		x			The tile's x-location. (x, y)
	///  @param[in]		y			The tile's y-location. (x, y)
	///  @param[out]	tiles		A pointer to an array of tiles that will hold the result.
	///  @param[in]		maxTiles	The maximum tiles the tiles parameter can hold.
	/// @return The number of tiles returned in the tiles array.
	getTilesAt(x, y int, maxTiles int) ([]*dtMeshTile, int)
	/// Gets the tile reference for the tile at specified grid location.
	///  @param[in]	x		The tile's x-location. (x, y, layer)
	///  @param[in]	y		The tile's y-location. (x, y, layer)
	///  @param[in]	layer	The tile's layer. (x, y, layer)
	/// @return The tile reference of the tile, or 0 if there is none.
	getTileRefAt(x, y, layer int) dtTileRef
	/// Gets the tile reference for the specified tile.
	///  @param[in]	tile	The tile.
	/// @return The tile reference of the tile.
	getTileRef(tile *dtMeshTile) dtTileRef
	/// Gets the tile for the specified tile reference.
	///  @param[in]	ref		The tile reference of the tile to retrieve.
	/// @return The tile for the specified reference, or null if the
	///		reference is invalid.
	getTileByRef(ref dtTileRef) *dtMeshTile
	/// The maximum number of tiles supported by the navigation mesh.
	/// @return The maximum number of tiles supported by the navigation mesh.
	getMaxTiles() int
	/// Gets the tile at the specified index.
	///  @param[in]	i		The tile index. [Limit: 0 >= index < #getMaxTiles()]
	/// @return The tile at the specified index.
	getTile(i int) *dtMeshTile
	/// Gets the tile and polygon for the specified polygon reference.
	///  @param[in]		ref		The reference for the a polygon.
	///  @param[out]	tile	The tile containing the polygon.
	///  @param[out]	poly	The polygon.
	/// @return The status flags for the operation.
	getTileAndPolyByRef(ref dtPolyRef) (tile *dtMeshTile, poly *dtPoly, status dtStatus)
	/// Returns the tile and polygon for the specified polygon reference.
	///  @param[in]		ref		A known valid reference for a polygon.
	///  @param[out]	tile	The tile containing the polygon.
	///  @param[out]	poly	The polygon.
	getTileAndPolyByRefUnsafe(ref dtPolyRef) (tile *dtMeshTile, poly *dtPoly)
	/// Checks the validity of a polygon reference.
	///  @param[in]	ref		The polygon reference to check.
	/// @return True if polygon reference is valid for the navigation mesh.
	isValidPolyRef(ref dtPolyRef) bool
	/// Gets the polygon reference for the tile's base polygon.
	///  @param[in]	tile		The tile.
	/// @return The polygon reference for the base polygon in the specified tile.
	getPolyRefBase(tile *dtMeshTile) dtPolyRef
	/// Gets the endpoints for an off-mesh connection, ordered by "direction of travel".
	///  @param[in]		prevRef		The reference of the polygon before the connection.
	///  @param[in]		polyRef		The reference of the off-mesh connection polygon.
	///  @param[out]	startPos	The start position of the off-mesh connection. [(x, y, z)]
	///  @param[out]	endPos		The end position of the off-mesh connection. [(x, y, z)]
	/// @return The status flags for the operation.
	getOffMeshConnectionPolyEndPoints(prevRef dtPolyRef, polyRef dtPolyRef, startPos []float64, endPos []float64) dtStatus
	/// Gets the specified off-mesh connection.
	///  @param[in]	ref		The polygon reference of the off-mesh connection.
	/// @return The specified off-mesh connection, or null if the polygon reference is not valid
	getOffMeshConnectionByRef(ref dtPolyRef) *dtOffMeshConnection
	/// @}

	/// @{
	/// @name State Management
	/// These functions do not effect #dtTileRef or #dtPolyRef's.

	/// Sets the user defined flags for the specified polygon.
	///  @param[in]	ref		The polygon reference.
	///  @param[in]	flags	The new flags for the polygon.
	/// @return The status flags for the operation.
	setPolyFlags(ref dtPolyRef, flags int) dtStatus

	/// Gets the user defined flags for the specified polygon.
	///  @param[in]		ref				The polygon reference.
	///  @param[out]	resultFlags		The polygon flags.
	/// @return The status flags for the operation.
	getPolyFlags(ref dtPolyRef) (resultFlags int, status dtStatus)

	/// Sets the user defined area for the specified polygon.
	///  @param[in]	ref		The polygon reference.
	///  @param[in]	area	The new area id for the polygon. [Limit: < #DT_MAX_AREAS]
	/// @return The status flags for the operation.
	setPolyArea(ref dtPolyRef, area int) dtStatus
	/// Gets the user defined area for the specified polygon.
	///  @param[in]		ref			The polygon reference.
	///  @param[out]	resultArea	The area id for the polygon.
	/// @return The status flags for the operation.
	getPolyArea(ref dtPolyRef) (resultArea int, status dtStatus)

	/// Gets the size of the buffer required by #storeTileState to store the specified tile's state.
	///  @param[in]	tile	The tile.
	/// @return The size of the buffer required to store the state.

	/// Restores the state of the tile.
	///  @param[in]	tile			The tile.
	///  @param[in]	data			The new state. (Obtained from #storeTileState.)
	///  @param[in]	maxDataSize		The size of the state within the data buffer.
	/// @return The status flags for the operation.
	restoreTileState(tile *dtMeshTile, tileState *dtTileState, polyStates []*dtPolyState, maxDataSize int) dtStatus
	/// @name Encoding and Decoding
	/// These functions are generally meant for internal use only.

	/// Derives a standard polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	salt	The tile's salt value.
	///  @param[in]	it		The index of the tile.
	///  @param[in]	ip		The index of the polygon within the tile.
	encodePolyId(salt, it, ip int) dtPolyRef
	/// Decodes a standard polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	ref   The polygon reference to decode.
	///  @param[out]	salt	The tile's salt value.
	///  @param[out]	it		The index of the tile.
	///  @param[out]	ip		The index of the polygon within the tile.
	///  @see #encodePolyId
	decodePolyId(ref dtPolyRef) (salt, it, ip int)
	/// Extracts a tile's salt value from the specified polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	ref		The polygon reference.
	///  @see #encodePolyId
	decodePolyIdPoly(ref dtPolyRef) int
	/// Extracts the tile's index from the specified polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	ref		The polygon reference.
	///  @see #encodePolyId
	decodePolyIdTile(ref dtPolyRef) int
	/// Extracts a tile's salt value from the specified polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	ref		The polygon reference.
	///  @see #encodePolyId
	decodePolyIdSalt(ref dtPolyRef) int
}

type dtNavMesh struct {
	m_params                  *dtNavMeshParams ///< Current initialization params. TODO: do not store this info twice.
	m_orig                    [3]float64       ///< Origin of the tile (0,0)
	m_tileWidth, m_tileHeight float64          ///< Dimensions of each tile.
	m_maxTiles                int              ///< Max number of tiles.
	m_tileLutSize             int              ///< Tile hash lookup size (must be pot).
	m_tileLutMask             int              ///< Tile hash lookup mask.
	m_posLookup               []*dtMeshTile    ///< Tile hash lookup.
	m_nextFree                *dtMeshTile      ///< Freelist of tiles.
	m_tiles                   []*dtMeshTile    ///< List of tiles.

	m_saltBits int ///< Number of salt bits in the tile ID.
	m_tileBits int ///< Number of tile bits in the tile ID.
	m_polyBits int ///< Number of poly bits in the tile ID.
}

func (mesh *dtNavMesh) getTile(i int) *dtMeshTile {
	return mesh.m_tiles[i]
}
func (mesh *dtNavMesh) getMaxTiles() int {
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
func (mesh *dtNavMesh) encodePolyId(salt, it, ip int) dtPolyRef {
	if DT_POLYREF64 == 1 {
		return salt<<(DT_POLY_BITS+DT_TILE_BITS) | (it << DT_POLY_BITS) | ip
	} else {
		return salt<<(mesh.m_polyBits+mesh.m_tileBits) | (it << mesh.m_polyBits) | ip
	}
}

func (mesh *dtNavMesh) decodePolyIdTile(ref dtPolyRef) int {
	if DT_POLYREF64 == 1 {
		tileMask := (1 << DT_TILE_BITS) - 1
		return int((ref >> DT_POLY_BITS) & tileMask)
	} else {
	}
	tileMask := (1 << mesh.m_tileBits) - 1
	return (ref >> mesh.m_polyBits) & tileMask
}

// / Decodes a standard polygon reference.
// /  @note This function is generally meant for internal use only.
// /  @param[in]	ref   The polygon reference to decode.
// /  @param[out]	salt	The tile's salt value.
// /  @param[out]	it		The index of the tile.
// /  @param[out]	ip		The index of the polygon within the tile.
// /  @see #encodePolyId
func (mesh *dtNavMesh) decodePolyId(ref dtPolyRef) (salt, it, ip int) {
	if DT_POLYREF64 == 1 {
		saltMask := (1 << DT_SALT_BITS) - 1
		tileMask := (1 << DT_TILE_BITS) - 1
		polyMask := (1 << DT_POLY_BITS) - 1
		salt = int((ref >> (DT_POLY_BITS + DT_TILE_BITS)) & saltMask)
		it = int((ref >> DT_POLY_BITS) & tileMask)
		ip = int(ref & polyMask)
	} else {
		saltMask := (1 << mesh.m_saltBits) - 1
		tileMask := (1 << mesh.m_tileBits) - 1
		polyMask := (1 << mesh.m_polyBits) - 1
		salt = int((ref >> (mesh.m_polyBits + mesh.m_tileBits)) & saltMask)
		it = int((ref >> mesh.m_polyBits) & tileMask)
		ip = int(ref & polyMask)
	}
	return
}

// / Get flags for edge in detail triangle.
// / @param[in]	triFlags		The flags for the triangle (last component of detail vertices above).
// / @param[in]	edgeIndex		The index of the first vertex of the edge. For instance, if 0,
// /								returns flags for edge AB.
func dtGetDetailTriEdgeFlags(triFlags int, edgeIndex int) int {
	return (triFlags >> (edgeIndex * 2)) & 0x3
}

func (mesh *dtNavMesh) decodePolyIdPoly(ref dtPolyRef) int {
	if DT_POLYREF64 == 1 {
		polyMask := (1 << DT_POLY_BITS) - 1
		return int(ref & polyMask)
	} else {
		polyMask := (1 << mesh.m_polyBits) - 1
		return int(ref & polyMask)
	}

}

// / Extracts a tile's salt value from the specified polygon reference.
// /  @note This function is generally meant for internal use only.
// /  @param[in]	ref		The polygon reference.
// /  @see #encodePolyId
func (mesh *dtNavMesh) decodePolyIdSalt(ref dtPolyRef) int {
	if DT_POLYREF64 == 1 {
		saltMask := (1 << DT_SALT_BITS) - 1
		return ((ref >> (DT_POLY_BITS + DT_TILE_BITS)) & saltMask)
	}
	saltMask := (1 << mesh.m_saltBits) - 1
	return ((ref >> (mesh.m_polyBits + mesh.m_tileBits)) & saltMask)

}

func (mesh *dtNavMesh) calcTileLoc(pos []float64) (tx, ty int) {
	tx = int(math.Floor((pos[0] - mesh.m_orig[0]) / mesh.m_tileWidth))
	ty = int(math.Floor((pos[2] - mesh.m_orig[2]) / mesh.m_tileHeight))
	return tx, ty
}

func (mesh *dtNavMesh) getTileByRef(ref dtTileRef) *dtMeshTile {
	if ref == 0 {
		return nil
	}

	tileIndex := mesh.decodePolyIdTile(dtPolyRef(ref))
	tileSalt := mesh.decodePolyIdSalt(dtPolyRef(ref))
	if tileIndex >= mesh.m_maxTiles {
		return nil
	}

	tile := mesh.m_tiles[tileIndex]
	if tile.salt != tileSalt {
		return nil
	}

	return tile
}
