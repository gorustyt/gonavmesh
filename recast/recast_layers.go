package recast

const RC_MAX_LAYERS_DEF = 63
const RC_MAX_NEIS_DEF = 16

// Keep type checking.
const RC_MAX_LAYERS = RC_MAX_LAYERS_DEF

const RC_MAX_NEIS = RC_MAX_NEIS_DEF

type rcLayerRegion struct {
	layers     [RC_MAX_LAYERS]int
	neis       [RC_MAX_NEIS]int
	ymin, ymax int
	layerId    int // Layer ID
	nlayers    int // Layer count
	nneis      int // Neighbour count
	base       int // Flag indicating if the region is the base of merged regions.
}

func contains(a []int, an int, v int) bool {
	n := an
	for i := 0; i < n; i++ {
		if a[i] == v {
			return true
		}

	}
	return false
}

func addUnique(a []int, an *int, anMax int, v int) bool {
	if contains(a, *an, v) {
		return true
	}

	if *an >= anMax {
		return false
	}

	a[*an] = v
	*an++
	return true
}

func LayersOverlapRange(amin, amax, bmin, bmax int) bool {
	if amin > bmax || amax < bmin {
		return false
	}
	return true
}

type rcLayerSweepSpan struct {
	ns  int // number samples
	id  int // region id
	nei int // neighbour id
}
/// Represents a set of heightfield layers.
/// @ingroup recast
/// @see rcAllocHeightfieldLayerSet, rcFreeHeightfieldLayerSet
 type  rcHeightfieldLayerSet struct {
layers *rcHeightfieldLayer ;			///< The layers in the set. [Size: #nlayers]
 nlayers int					///< The number of layers in the set.
};

/// Represents a heightfield layer within a layer set.
/// @see rcHeightfieldLayerSet
 type  rcHeightfieldLayer struct {
bmin[3]float64;				///< The minimum bounds in world space. [(x, y, z)]
 bmax[3]float64;				///< The maximum bounds in world space. [(x, y, z)]
cs float64;					///< The size of each cell. (On the xz-plane.)
 ch float64;					///< The height of each cell. (The minimum increment along the y-axis.)
 width int				///< The width of the heightfield. (Along the x-axis in cell units.)
height int 			///< The height of the heightfield. (Along the z-axis in cell units.)
minx int 			///< The minimum x-bounds of usable data.
 maxx int		///< The maximum x-bounds of usable data.
 miny int				///< The minimum y-bounds of usable data. (Along the z-axis.)
 maxy int				///< The maximum y-bounds of usable data. (Along the z-axis.)
 hmin int			///< The minimum height bounds of usable data. (Along the y-axis.)
 hmax int				///< The maximum height bounds of usable data. (Along the y-axis.)
 heights []int;		///< The heightfield. [Size: width * height]
areas []int ;		///< Area ids. [Size: Same as #heights]
 cons []int ;		///< Packed neighbor connection information. [Size: Same as #heights]
};


/// @par
///
/// See the #rcConfig documentation for more information on the configuration parameters.
///
/// @see rcAllocHeightfieldLayerSet, rcCompactHeightfield, rcHeightfieldLayerSet, rcConfig
 func  rcBuildHeightfieldLayers( chf* rcCompactHeightfield , borderSize,  walkableHeight int , lset *rcHeightfieldLayerSet )bool {


 w := chf.width;
 h := chf.height;
	 srcReg:=make([]int,chf.spanCount)
for i:=range srcReg{
	srcReg[i]=0xff
}

 nsweeps := chf.width;
 sweeps:=make([]*rcLayerSweepSpan,nsweeps)


// Partition walkable area into monotone regions.
prevCount:=make([]int,256);
regId := 0;

for  y := borderSize; y < h-borderSize; y++{
sweepId := 0;

for  x := borderSize; x < w-borderSize; x++{
 c := chf.cells[x+y*w];

		 i := c.index
		 ni := (c.index+c.count)
for ; i < ni; i++{
 s := chf.spans[i];
if (chf.areas[i] == RC_NULL_AREA) {continue;}
 sid := 0xff;

// -x
if (rcGetCon(s, 0) != RC_NOT_CONNECTED) {
 ax := x + rcGetDirOffsetX(0);
 ay := y + rcGetDirOffsetY(0);
 ai := chf.cells[ax+ay*w].index + rcGetCon(s, 0);
if (chf.areas[ai] != RC_NULL_AREA && srcReg[ai] != 0xff){sid = srcReg[ai];}

}

if (sid == 0xff) {
sid = sweepId;
	sweepId++
sweeps[sid].nei = 0xff;
sweeps[sid].ns = 0;
}

// -y
if (rcGetCon(s,3) != RC_NOT_CONNECTED) {
 ax := x + rcGetDirOffsetX(3);
 ay := y + rcGetDirOffsetY(3);
 ai := chf.cells[ax+ay*w].index + rcGetCon(s, 3);
nr := srcReg[ai];
if (nr != 0xff) {
// Set neighbour when first valid neighbour is encoutered.
if (sweeps[sid].ns == 0){sweeps[sid].nei = nr;}


if (sweeps[sid].nei == nr) {
// Update existing neighbour
sweeps[sid].ns++;
prevCount[nr]++;
} else {
// This is hit if there is nore than one neighbour.
// Invalidate the neighbour.
sweeps[sid].nei = 0xff;
}
}
}

srcReg[i] = sid;
}
}

// Create unique ID.
for  i := 0; i < sweepId; i++{
// If the neighbour is set and there is only one continuous connection to it,
// the sweep will be merged with the previous one, else new region is created.
if (sweeps[i].nei != 0xff && prevCount[sweeps[i].nei] == (int)sweeps[i].ns)
{
sweeps[i].id = sweeps[i].nei;
}
else
{
if (regId == 255)
{
ctx->log(RC_LOG_ERROR, "rcBuildHeightfieldLayers: Region ID overflow.");
return false;
}
sweeps[i].id = regId++;
}
}

// Remap local sweep ids to region ids.
for (int x = borderSize; x < w-borderSize; ++x)
{
const rcCompactCell& c = chf.cells[x+y*w];
for (int i = (int)c.index, ni = (int)(c.index+c.count); i < ni; ++i)
{
if (srcReg[i] != 0xff)
srcReg[i] = sweeps[srcReg[i]].id;
}
}
}

// Allocate and init layer regions.
const int nregs = (int)regId;
rcScopedDelete<rcLayerRegion> regs((rcLayerRegion*)rcAlloc(sizeof(rcLayerRegion)*nregs, RC_ALLOC_TEMP));
if (!regs)
{
ctx->log(RC_LOG_ERROR, "rcBuildHeightfieldLayers: Out of memory 'regs' (%d).", nregs);
return false;
}
memset(regs, 0, sizeof(rcLayerRegion)*nregs);
for (int i = 0; i < nregs; ++i)
{
regs[i].layerId = 0xff;
regs[i].ymin = 0xffff;
regs[i].ymax = 0;
}

// Find region neighbours and overlapping regions.
for (int y = 0; y < h; ++y)
{
for (int x = 0; x < w; ++x)
{
const rcCompactCell& c = chf.cells[x+y*w];

unsigned char lregs[RC_MAX_LAYERS];
int nlregs = 0;

for (int i = (int)c.index, ni = (int)(c.index+c.count); i < ni; ++i)
{
const rcCompactSpan& s = chf.spans[i];
const unsigned char ri = srcReg[i];
if (ri == 0xff) continue;

regs[ri].ymin = rcMin(regs[ri].ymin, s.y);
regs[ri].ymax = rcMax(regs[ri].ymax, s.y);

// Collect all region layers.
if (nlregs < RC_MAX_LAYERS)
lregs[nlregs++] = ri;

// Update neighbours
for (int dir = 0; dir < 4; ++dir)
{
if (rcGetCon(s, dir) != RC_NOT_CONNECTED)
{
const int ax = x + rcGetDirOffsetX(dir);
const int ay = y + rcGetDirOffsetY(dir);
const int ai = (int)chf.cells[ax+ay*w].index + rcGetCon(s, dir);
const unsigned char rai = srcReg[ai];
if (rai != 0xff && rai != ri)
{
// Don't check return value -- if we cannot add the neighbor
// it will just cause a few more regions to be created, which
// is fine.
addUnique(regs[ri].neis, regs[ri].nneis, RC_MAX_NEIS, rai);
}
}
}

}

// Update overlapping regions.
for (int i = 0; i < nlregs-1; ++i)
{
for (int j = i+1; j < nlregs; ++j)
{
if (lregs[i] != lregs[j])
{
rcLayerRegion& ri = regs[lregs[i]];
rcLayerRegion& rj = regs[lregs[j]];

if (!addUnique(ri.layers, ri.nlayers, RC_MAX_LAYERS, lregs[j]) ||
!addUnique(rj.layers, rj.nlayers, RC_MAX_LAYERS, lregs[i]))
{
ctx->log(RC_LOG_ERROR, "rcBuildHeightfieldLayers: layer overflow (too many overlapping walkable platforms). Try increasing RC_MAX_LAYERS.");
return false;
}
}
}
}

}
}

// Create 2D layers from regions.
unsigned char layerId = 0;

static const int MAX_STACK = 64;
unsigned char stack[MAX_STACK];
int nstack = 0;

for (int i = 0; i < nregs; ++i)
{
rcLayerRegion& root = regs[i];
// Skip already visited.
if (root.layerId != 0xff)
continue;

// Start search.
root.layerId = layerId;
root.base = 1;

nstack = 0;
stack[nstack++] = (unsigned char)i;

while (nstack)
{
// Pop front
rcLayerRegion& reg = regs[stack[0]];
nstack--;
for (int j = 0; j < nstack; ++j)
stack[j] = stack[j+1];

const int nneis = (int)reg.nneis;
for (int j = 0; j < nneis; ++j)
{
const unsigned char nei = reg.neis[j];
rcLayerRegion& regn = regs[nei];
// Skip already visited.
if (regn.layerId != 0xff)
continue;
// Skip if the neighbour is overlapping root region.
if (contains(root.layers, root.nlayers, nei))
continue;
// Skip if the height range would become too large.
const int ymin = rcMin(root.ymin, regn.ymin);
const int ymax = rcMax(root.ymax, regn.ymax);
if ((ymax - ymin) >= 255)
continue;

if (nstack < MAX_STACK)
{
// Deepen
stack[nstack++] = (unsigned char)nei;

// Mark layer id
regn.layerId = layerId;
// Merge current layers to root.
for (int k = 0; k < regn.nlayers; ++k)
{
if (!addUnique(root.layers, root.nlayers, RC_MAX_LAYERS, regn.layers[k]))
{
ctx->log(RC_LOG_ERROR, "rcBuildHeightfieldLayers: layer overflow (too many overlapping walkable platforms). Try increasing RC_MAX_LAYERS.");
return false;
}
}
root.ymin = rcMin(root.ymin, regn.ymin);
root.ymax = rcMax(root.ymax, regn.ymax);
}
}
}

layerId++;
}

// Merge non-overlapping regions that are close in height.
const unsigned short mergeHeight = (unsigned short)walkableHeight * 4;

for (int i = 0; i < nregs; ++i)
{
rcLayerRegion& ri = regs[i];
if (!ri.base) continue;

unsigned char newId = ri.layerId;

for (;;)
{
unsigned char oldId = 0xff;

for (int j = 0; j < nregs; ++j)
{
if (i == j) continue;
rcLayerRegion& rj = regs[j];
if (!rj.base) continue;

// Skip if the regions are not close to each other.
if (!overlapRange(ri.ymin,ri.ymax+mergeHeight, rj.ymin,rj.ymax+mergeHeight))
continue;
// Skip if the height range would become too large.
const int ymin = rcMin(ri.ymin, rj.ymin);
const int ymax = rcMax(ri.ymax, rj.ymax);
if ((ymax - ymin) >= 255)
continue;

// Make sure that there is no overlap when merging 'ri' and 'rj'.
bool overlap = false;
// Iterate over all regions which have the same layerId as 'rj'
for (int k = 0; k < nregs; ++k)
{
if (regs[k].layerId != rj.layerId)
continue;
// Check if region 'k' is overlapping region 'ri'
// Index to 'regs' is the same as region id.
if (contains(ri.layers,ri.nlayers, (unsigned char)k))
{
overlap = true;
break;
}
}
// Cannot merge of regions overlap.
if (overlap)
continue;

// Can merge i and j.
oldId = rj.layerId;
break;
}

// Could not find anything to merge with, stop.
if (oldId == 0xff)
break;

// Merge
for (int j = 0; j < nregs; ++j)
{
rcLayerRegion& rj = regs[j];
if (rj.layerId == oldId)
{
rj.base = 0;
// Remap layerIds.
rj.layerId = newId;
// Add overlaid layers from 'rj' to 'ri'.
for (int k = 0; k < rj.nlayers; ++k)
{
if (!addUnique(ri.layers, ri.nlayers, RC_MAX_LAYERS, rj.layers[k]))
{
ctx->log(RC_LOG_ERROR, "rcBuildHeightfieldLayers: layer overflow (too many overlapping walkable platforms). Try increasing RC_MAX_LAYERS.");
return false;
}
}

// Update height bounds.
ri.ymin = rcMin(ri.ymin, rj.ymin);
ri.ymax = rcMax(ri.ymax, rj.ymax);
}
}
}
}

// Compact layerIds
unsigned char remap[256];
memset(remap, 0, 256);

// Find number of unique layers.
layerId = 0;
for (int i = 0; i < nregs; ++i)
remap[regs[i].layerId] = 1;
for (int i = 0; i < 256; ++i)
{
if (remap[i])
remap[i] = layerId++;
else
remap[i] = 0xff;
}
// Remap ids.
for (int i = 0; i < nregs; ++i)
regs[i].layerId = remap[regs[i].layerId];

// No layers, return empty.
if (layerId == 0)
return true;

// Create layers.
rcAssert(lset.layers == 0);

const int lw = w - borderSize*2;
const int lh = h - borderSize*2;

// Build contracted bbox for layers.
float bmin[3], bmax[3];
rcVcopy(bmin, chf.bmin);
rcVcopy(bmax, chf.bmax);
bmin[0] += borderSize*chf.cs;
bmin[2] += borderSize*chf.cs;
bmax[0] -= borderSize*chf.cs;
bmax[2] -= borderSize*chf.cs;

lset.nlayers = (int)layerId;

lset.layers = (rcHeightfieldLayer*)rcAlloc(sizeof(rcHeightfieldLayer)*lset.nlayers, RC_ALLOC_PERM);
if (!lset.layers)
{
ctx->log(RC_LOG_ERROR, "rcBuildHeightfieldLayers: Out of memory 'layers' (%d).", lset.nlayers);
return false;
}
memset(lset.layers, 0, sizeof(rcHeightfieldLayer)*lset.nlayers);


// Store layers.
for (int i = 0; i < lset.nlayers; ++i)
{
unsigned char curId = (unsigned char)i;

rcHeightfieldLayer* layer = &lset.layers[i];

const int gridSize = sizeof(unsigned char)*lw*lh;

layer->heights = (unsigned char*)rcAlloc(gridSize, RC_ALLOC_PERM);
if (!layer->heights)
{
ctx->log(RC_LOG_ERROR, "rcBuildHeightfieldLayers: Out of memory 'heights' (%d).", gridSize);
return false;
}
memset(layer->heights, 0xff, gridSize);

layer->areas = (unsigned char*)rcAlloc(gridSize, RC_ALLOC_PERM);
if (!layer->areas)
{
ctx->log(RC_LOG_ERROR, "rcBuildHeightfieldLayers: Out of memory 'areas' (%d).", gridSize);
return false;
}
memset(layer->areas, 0, gridSize);

layer->cons = (unsigned char*)rcAlloc(gridSize, RC_ALLOC_PERM);
if (!layer->cons)
{
ctx->log(RC_LOG_ERROR, "rcBuildHeightfieldLayers: Out of memory 'cons' (%d).", gridSize);
return false;
}
memset(layer->cons, 0, gridSize);

// Find layer height bounds.
int hmin = 0, hmax = 0;
for (int j = 0; j < nregs; ++j)
{
if (regs[j].base && regs[j].layerId == curId)
{
hmin = (int)regs[j].ymin;
hmax = (int)regs[j].ymax;
}
}

layer->width = lw;
layer->height = lh;
layer->cs = chf.cs;
layer->ch = chf.ch;

// Adjust the bbox to fit the heightfield.
rcVcopy(layer->bmin, bmin);
rcVcopy(layer->bmax, bmax);
layer->bmin[1] = bmin[1] + hmin*chf.ch;
layer->bmax[1] = bmin[1] + hmax*chf.ch;
layer->hmin = hmin;
layer->hmax = hmax;

// Update usable data region.
layer->minx = layer->width;
layer->maxx = 0;
layer->miny = layer->height;
layer->maxy = 0;

// Copy height and area from compact heightfield.
for (int y = 0; y < lh; ++y)
{
for (int x = 0; x < lw; ++x)
{
const int cx = borderSize+x;
const int cy = borderSize+y;
const rcCompactCell& c = chf.cells[cx+cy*w];
for (int j = (int)c.index, nj = (int)(c.index+c.count); j < nj; ++j)
{
const rcCompactSpan& s = chf.spans[j];
// Skip unassigned regions.
if (srcReg[j] == 0xff)
continue;
// Skip of does nto belong to current layer.
unsigned char lid = regs[srcReg[j]].layerId;
if (lid != curId)
continue;

// Update data bounds.
layer->minx = rcMin(layer->minx, x);
layer->maxx = rcMax(layer->maxx, x);
layer->miny = rcMin(layer->miny, y);
layer->maxy = rcMax(layer->maxy, y);

// Store height and area type.
const int idx = x+y*lw;
layer->heights[idx] = (unsigned char)(s.y - hmin);
layer->areas[idx] = chf.areas[j];

// Check connection.
unsigned char portal = 0;
unsigned char con = 0;
for (int dir = 0; dir < 4; ++dir)
{
if (rcGetCon(s, dir) != RC_NOT_CONNECTED)
{
const int ax = cx + rcGetDirOffsetX(dir);
const int ay = cy + rcGetDirOffsetY(dir);
const int ai = (int)chf.cells[ax+ay*w].index + rcGetCon(s, dir);
unsigned char alid = srcReg[ai] != 0xff ? regs[srcReg[ai]].layerId : 0xff;
// Portal mask
if (chf.areas[ai] != RC_NULL_AREA && lid != alid)
{
portal |= (unsigned char)(1<<dir);
// Update height so that it matches on both sides of the portal.
const rcCompactSpan& as = chf.spans[ai];
if (as.y > hmin)
layer->heights[idx] = rcMax(layer->heights[idx], (unsigned char)(as.y - hmin));
}
// Valid connection mask
if (chf.areas[ai] != RC_NULL_AREA && lid == alid)
{
const int nx = ax - borderSize;
const int ny = ay - borderSize;
if (nx >= 0 && ny >= 0 && nx < lw && ny < lh)
con |= (unsigned char)(1<<dir);
}
}
}

layer->cons[idx] = (portal << 4) | con;
}
}
}

if (layer->minx > layer->maxx)
layer->minx = layer->maxx = 0;
if (layer->miny > layer->maxy)
layer->miny = layer->maxy = 0;
}

return true;
}
