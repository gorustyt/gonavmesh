package recast

type DetourTitleCacheData struct {
}

type DetourTitleCacheLayerData struct {
	Header   *DtTileCacheLayerHeader
	Comp     DtTileCacheCompressor //压缩器
	Heights  []int
	Areas    []int
	Cons     []int
	Regs     []int
	RegCount int
}
